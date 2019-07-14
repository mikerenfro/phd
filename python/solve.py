# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 10:48:58 2017

@author: Renfro
"""

import datetime
import numpy as np
import os
import scipy
import scipy.interpolate
import scipy.optimize
import scipy.stats
import preprocess as pre
import postprocess as post
import section
import subprocess
import util


def calculate_rphi(phi, crack_geometry, plate_geometry):
    """Calculate the r_phi values for a given crack and plate geometry at
    a particular angle phi.
    """
    t = plate_geometry['t']
    W = plate_geometry['W']
    a = crack_geometry['a']
    c = crack_geometry['c']
    # Standard ellipse function: A*x^2 + B*y^2 + C*x + D*y + E*x*y + F = 0
    # In terms of crack geometry: (x^2)/(c^2)+(y^2)/(a^2) = 1
    # This, A = 1/(c^2), B=1/(a^2), C = D = E = 0, F = -1
    ellipse_A = 1/(c*c)
    ellipse_B = 1/(a*a)
    (ellipse_C, ellipse_D, ellipse_E, ellipse_F) = (0, 0, 0, -1) # analysis:ignore
    # Calculate point on circle inscribed in ellipse
    x_c = a*np.cos(phi) # analysis:ignore
    y_c = a*np.sin(phi)
    # Calculate point on ellipse
    x_e = c*np.cos(phi)
    y_e = y_c
    # print("x_c, y_c, x_e, y_e:", x_c, y_c, x_e, y_e)
    # Find slope of normal to point on ellipse at angle phi
    slope = (2*ellipse_B*y_e+ellipse_D*ellipse_E*x_e
             )/(
            2*ellipse_A*x_e+ellipse_E*y_e)
    # Find equation of line with given slope passing through (x_e, y_e)
    # (y-y0) = m*(x-x0)
    # (x-x0) = (y-y0)/m
    # x = (y-y0)/m + x0
    if not (phi == 0):
        y_range = np.array([0, y_e, t])
        x_range = (y_range-y_e)/slope + x_e
    else:
        x_range = np.array([0, c, W/2])
        y_range = slope*x_range
    r_phia = np.sqrt(np.power((x_range[1]-x_range[0]), 2) +
                     np.power((y_range[1]-y_range[0]), 2))
    r_phib = np.sqrt(np.power((x_range[2]-x_range[1]), 2) +
                     np.power((y_range[2]-y_range[1]), 2))
    return (r_phia, r_phib)


def calculate_J_final(J_table, phi_table, phi):
    """Return a J value for the last load step at a given angle phi.
    """
    J_final_int = scipy.interpolate.interp1d(phi_table, J_table[-1, :])
    J_value = J_final_int(phi)
    return J_value


def calculate_M(J_table, phi_table, crack_geometry, plate_geometry, material):
    """Calculate M = (r_phi*st_far)/J at phi=30 and phi=90.
    Tension: st_far = reaction_forces/area
    Bending: st_far = (bending_moment*thickness)/(2*area_moment_of_inertia)
    """
    Sys = material['Sys']
    (r_a_30, r_b_30) = calculate_rphi(np.pi/6, crack_geometry, plate_geometry)
    (r_a_90, r_b_90) = calculate_rphi(np.pi/2, crack_geometry, plate_geometry)
    M_a_30 = (r_a_30*Sys)/calculate_J_final(J_table, phi_table, np.pi/6)
    M_b_30 = (r_b_30*Sys)/calculate_J_final(J_table, phi_table, np.pi/6)
    M_a_90 = (r_a_90*Sys)/J_table[-1, -1]
    M_b_90 = (r_b_90*Sys)/J_table[-1, -1]
    return (M_a_30, M_b_30, M_a_90, M_b_90)


def M_deviation(M):
    """By default, returns the difference between dimensionless M value and
    the middle of the range of [20, 25]. Returns 0 if in the range of [20, 25].
    """
    # return deviation of M from range of [20, 25].
    M_range = np.array([20, 25])
    if M > np.min(M_range) and M < np.max(M_range):
        return 0
    else:
        return np.abs(M-np.median(M_range))


def slope_M_objective(cmod, J, phi_table, crack_geometry, plate_geometry,
                      material):
    """Objective function to determine when a plate has been deformed
    enough to extrapolate from. Allen, Wells (2014) used M=rphi*Sys/J,
    and indicated that when M<20 or so, models were pulled far enough.
    About half the models in the TASC database had M>25, and about 12%
    had M>50.

    This attempts to improve the linear extrapolation metric.
    Qualitatively, if the slope of the last few segments of the J-CMOD
    curve are similar enough, extrapolation should be accurate.
    But if we've only pulled into the linear regime, we're not predicting
    elastic-plastic behavior. We want the J-CMOD curve to be relatively
    steep when we terminate.

    The new objective function is M*(slope percent change)/(avg slope)
    where the slopes are taken from the J-CMOD curve (80-90% range of
    CMOD values, and 90-100$ range of CMOD values.), M is taken at the
    last time step at phi=30 and phi=90.

    This drives the objective toward zero if slopes are high, if slopes
    are consistent, and if M is low.

    If the slopes differ by <5% and M<50, returns 0.
    """
    cmod = cmod.values
    breakpoints = (0.6, 0.8)
    range2 = np.where((cmod >= breakpoints[0]*np.max(cmod)) &
                      (cmod < breakpoints[1]*np.max(cmod)))
    range3 = np.where((cmod >= breakpoints[1]*np.max(cmod)))

    (slope2, intercept2, r_value2,
     p_value2, std_err2) = scipy.stats.linregress(cmod[range2],
                                                  J[range2, -1][0])
    (slope3, intercept3, r_value3,
     p_value3, std_err3) = scipy.stats.linregress(cmod[range3],
                                                  J[range3, -1][0])
    pct_change = 100*(1-slope2/slope3)

    (M_a_30, M_b_30, M_a_90, M_b_90) = calculate_M(J,
                                                   phi_table,
                                                   crack_geometry,
                                                   plate_geometry,
                                                   material)
    min_M = np.min(np.array([M_a_30, M_b_30, M_a_90, M_b_90]))
    objective = min_M*pct_change/(0.5*(slope2+slope3))
    print(("min(M) = {0:.4e}, slope2 = {1:.4e}, "
           "slope3 = {2:.4e}, ({3:.1f} % change), ").format(min_M, slope2,
                                                            slope3,
                                                            pct_change),
          end='')
    if np.abs(pct_change) < 5 and min_M < 50:
        return 0
    else:
        return objective


def feacrack_increment_calc_scalar(step, num_steps):
    """Calculate a scalar to apply to a boundary condition at a given
    load step. Matches the FEACrack default behavior for Increment Calc Mode
    in the General tab of the Boundary Conditions dialog.
    """
    x = np.array([0.00, 0.30, 0.40, 0.50, 0.70, 0.90, 1.00])
    y = np.array([0.00, 0.75, 0.85, 0.90, 0.95, 0.99, 1.00])
    interpolant = scipy.interpolate.interp1d(x, y, kind='linear')
    scalar = interpolant(float(step)/float(num_steps))
    # don't need float() on Python 3.x, but doesn't hurt
    return scalar


def modify_bc(bc, input_file, model_type):
    """Adjust the boundary condition applied to a tension or bending model.
    Returns nothing, raises a ValueError if model_type is not 'tension'
    or 'bending'.
    """

    # bc = bc[0]  # only needed for minimize, not for newton
    print("BC = {0:.3e}: ".format(bc), end='', flush=True)
    with open(input_file) as f:
        input_lines = f.readlines()

    if model_type == 'tension':
        # If we're a tension case, want to adjust non-zero w constraints to
        # some value.
        # - Find line starting with "constraints"
        # - Find next line starting with "c *echo off"
        # - In that range, adjust lines like
        #   "  10294  u 0.000000E+00  w 1.000000E+00" to a new w value
        start_index = util.find_first_line_matching(
                input_lines, "constraints")+1
        end_index = util.find_first_line_matching(
                input_lines, "c *echo off", start_index)
        for i in range(start_index, end_index):
            line = input_lines[i]
            if 'w' in line:
                (rest, w_value) = line.rsplit(maxsplit=1)
                w_value = float(w_value)
                if not (w_value == 0):
                    input_lines[i] = "{0} {1:.6e}\n".format(rest, bc)
    elif model_type == 'bending':
        # If we're a bending case, want to adjust element tractions to some
        # value.
        # - Find line starting with "loading set1"
        # - Find next line starting with "c"
        # - In that range, adjust lines like
        #   "   2571 face  2  traction ty    3.52000000E-01" to a new ty value
        start_index = util.find_first_line_matching(
                input_lines, "loading set1")+1
        end_index = util.find_first_line_matching(
                input_lines, "c", start_index)
        for i in range(start_index, end_index):
            line = input_lines[i]
            if 'traction ty' in line:
                (rest, ty_value) = line.rsplit(maxsplit=1)
                ty_value = float(ty_value)
                if not (ty_value == 0):
                    input_lines[i] = "{0} {1:.8e}\n".format(rest, bc)
    else:
        raise(ValueError, "model_type must be either 'tension' or 'bending'")

    with open(input_file, 'w') as f:
        f.writelines(input_lines)
    # print("Using BC of {0:.8e}".format(bc))


def run_model(bc, input_file, crack_geometry, plate_geometry, material,
              model_type):
    """Run a model with a given overall boundary condition value. Return
    a value corresponding to its J-CMOD fitness (using J_CMOD_objective()).
    """

    if bc <= 0:
        # Negative BCs not allowed. Drive objective up enough to get BC back
        # to positive values. +0.1 to keep from driving BC back to 0.
        objective = -bc+0.1
    else:
        # Adjust boundary condition in the WARP3D input file
        modify_bc(bc, input_file, model_type)

        input_basename = os.path.splitext(input_file)[0]
        output_file = "{0}.out".format(input_basename)
        # bpf_filename = input_basename[:-4]+".bpf"
        bpf_filename = pre.get_bpf_filename(input_file)

        # Remove BPF file (if it exists) before running
        if os.path.exists(bpf_filename):
            os.remove(bpf_filename)

        # Run WARP3D
        subprocess.run('warp3d.omp < {0} > {1} 2>&1'.format(
                input_file, output_file), shell=True)

        with open(output_file, 'r', errors='ignore') as f:
            # On a major crash, had bad bytes in the output file.
            # https://stackoverflow.com/a/41652865
            output_lines = f.readlines()
        if util.find_first_line_containing(output_lines,
                                           '>>>>> analysis terminated.') > -1:
            # Exceeded strain limits, most likely. Drive objective way up to
            # get BC back down to normal levels. -0.01 to keep from driving
            # BC negative.
            objective = bc-0.01
        else:
            with open(input_file) as f:
                input_lines = f.readlines()

            (phi, J_table) = post.find_J(input_lines, output_lines)

            output_filenames = post.extract_results(bpf_filename)
            msh_filename = pre.get_msh_filename(input_file)
            coordinates = post.get_node_coordinates(msh_filename)
            displacement_name = output_filenames['displacements']
            cmod = post.get_cmod(coordinates, displacement_name, x=0, y=0, z=0)
            J_30 = post.J_interpolated(phi, J_table, np.pi/6)
            objective = J_CMOD_objective(J_30, cmod)/bc

    print("objective = {0:.4e}, ".format(objective), end='')
    print(datetime.datetime.now().isoformat(' ', timespec='seconds'))
    return objective


def J_CMOD_objective(J, cmod, num_points=None):
    """Want consistent slope for the later part of the J-CMOD curve,
    and a much shallower slope for the earlier (linear part). Also don't
    want to pull so far as to exceed the stress-strain curve data.

    So return a value proportional to change in slope for
    last two segments, and ratio of slope from first segment to last.

    If the slope changes by <10%, and the slope ratio is <1:20, return 0.
    """
    # Examine slope at [60%, 80%) and [80%, 100%] ranges
    breakpoints = (None, 0.6, 0.8)

    # Examine first several points, or all points
    if not (num_points is None):
        J = J[:num_points]
        cmod = cmod[:num_points]
    # Normalize J and CMOD so plot is bounded by (0, 0) and (1, 1)
    J = np.insert(J, 0, 0)
    cmod = np.insert(cmod.values, 0, 0)
    J = J/np.max(J)
    cmod = cmod/np.max(cmod)
#    print("J")
#    print(J)
#    print("CMOD")
#    print(cmod)

    # Calculate slope of initial linear region of J-CMOD curve
    slope1 = (J[1]-J[0])/(cmod[1]-cmod[0])

    # Make finer interpolation of CMOD and J to prevent too few points
    # from being used in linear regressions
    J_int_function = scipy.interpolate.interp1d(cmod, J)
    cmod_int = np.linspace(0, np.max(cmod))  # 50 interpolation points
    J_int = J_int_function(cmod_int)

    # Define ranges of points for in mostly-plastic regime
    range2 = np.where((cmod_int >= breakpoints[1]) &
                      (cmod_int < breakpoints[2]))
    range3 = np.where((cmod_int >= breakpoints[2]))

    # Find slopes in mostly-plastic ranges
    (slope2, intercept2, r_value2,
     p_value2, std_err2) = scipy.stats.linregress(cmod_int[range2],
                                                  J_int[range2])
    (slope3, intercept3, r_value3,
     p_value3, std_err3) = scipy.stats.linregress(cmod_int[range3],
                                                  J_int[range3])
    # Calculate quality of curve: change in slope in mostly-plastic regime
    # and ratio of last range to initial linear region
    slope_change23 = np.abs(1-slope2/slope3)
    slope_ratio31 = np.abs(slope3/slope1)
    print(("slope1 = {0:.4e}, "
           "slope2 = {1:.4e}, "
           "slope3 = {2:.4e} "
           "({3:.1f}%, {4:.1f}:1), ").format(slope1,
                                             slope2,
                                             slope3,
                                             100*slope_change23,
                                             slope_ratio31),
          end='')
    if (slope_change23 < 0.10) and (slope_ratio31 > 20):
        objective = 0
    else:
        objective = slope_change23/slope_ratio31

    return objective


def St_CMOD_objective(St_far, cmod):
    """Want consistent slope for the later part of the Stress-CMOD curve,
    and a much steeper slope for the earlier (linear part). Also don't
    want to pull so far as to exceed the stress-strain curve data.

    So return a value proportional to max(CMOD), change in slope for
    last two segments, and ratio of slope from last segment to first.
    """
    cmod = np.insert(cmod.values, 0, 0)
    St_far = np.insert(St_far, 0, 0)
    breakpoints = (0.1, 0.6, 0.8)
    range1 = np.where((cmod <= breakpoints[0]*np.max(cmod)))
    range2 = np.where((cmod >= breakpoints[1]*np.max(cmod)) &
                      (cmod < breakpoints[2]*np.max(cmod)))
    range3 = np.where((cmod >= breakpoints[2]*np.max(cmod)))

    (slope1, intercept1, r_value1,
     p_value1, std_err1) = scipy.stats.linregress(cmod[range1],
                                                  St_far[range1])
    (slope2, intercept2, r_value2,
     p_value2, std_err2) = scipy.stats.linregress(cmod[range2],
                                                  St_far[range2])
    (slope3, intercept3, r_value3,
     p_value3, std_err3) = scipy.stats.linregress(cmod[range3],
                                                  St_far[range3])
    print(("slope1 = {0:.4e}, slope2 = {1:.4e}, "
           "slope3 = {2:.4e}, max(CMOD) = {3:.4e}, ").format(slope1,
                                                             slope2,
                                                             slope3,
                                                             np.max(cmod)),
          end='')
#    slope_change23 = np.abs(1-slope2/slope3)
#    slope_ratio31 = np.abs(slope3/slope1)
#    objective = np.max(cmod)*slope_change23*slope_ratio31
    objective = np.abs(slope3/slope2)*np.abs(slope3/slope1)

    return objective


def slope_objective(min_M, slope1, slope2):
    percent_change = np.abs((1-slope2/slope1)*100)
    if percent_change < 5 and min_M < 50:
        return 0
    else:
        return (min_M*percent_change)/(0.5*(slope1+slope2))


def find_pin_locations(input_file):
    """Fnid the locations of roller pins for 4-point bending models.
    returns the z coordinate values of the inner and outer pin locations.
    """
    with open(input_file) as f:
        input_lines = f.readlines()

    # Bending model input files will have pins
    # e.g.,
    #
    # c    top load pin mesh line Reaction Node
    # c    number of nodes in list =       33
    # c        3452    0.000000E+00    0.000000E+00   -1.250000E+00
    # ...
    start_index = util.find_first_line_matching(
            input_lines, "c    top load pin mesh line Reaction Node")
    first_node_line = input_lines[start_index+2]
    (rest, pin_location) = first_node_line.rsplit(maxsplit=1)
    outer_pin_location = float(pin_location)
    #
    # and
    #
    # c    bottom load pin mesh line
    # c    number of nodes in list =       25
    # c        3299    1.000000E+00   -1.000000E+00   -6.250000E-01
    # ...
    start_index = util.find_first_line_matching(
            input_lines, "c    bottom load pin mesh line")
    first_node_line = input_lines[start_index+2]
    (rest, pin_location) = first_node_line.rsplit(maxsplit=1)
    inner_pin_location = float(pin_location)
    return (inner_pin_location, outer_pin_location)


def get_initial_bc(input_file, crack_geometry, plate_geometry, material,
                   model_type):
    """Set the initial boundary condition for optimization.
    Detect if the model is in tension or bending, and set the appropriate
    BC to X% of level causing net section yielding.
    """
    E = material['E']
    Sys = material['Sys']
    L = plate_geometry['L']
    Sys_ratio = 0.9
    if model_type == 'tension':
        # axial deflection required to cause remote stress of 0.5*Sys
        # u = (P*L)/(A*E) = (P/A)*(L/E) = sigma*(L/E)
        init_bc = Sys_ratio*float(Sys)*float(L)/float(E)
    else:
        # transverse traction along pin line required to cause remote stress
        # of Sys_ratio*Sys (Sys_ratio can be 0.5, 0.9, etc.)
        # Fundamental equations:
        # sigma = -(M*ybar)/I
        # M(x) = P*(z_outer-z_inner)
        # P = traction*W*t
        #
        # become:
        # M = -(I*sigma)/ybar
        # P*(z_outer-z_inner) = -(I*sigma)/ybar
        # traction*W*t*(z_outer-z_inner) = -(I*sigma)/ybar
        # traction = -(I*sigma)/(ybar*t*W*(z_outer-z_inner))
        W = plate_geometry['W']
        t = plate_geometry['t']
        z_inner = np.abs(plate_geometry['z_inner'])
        z_outer = np.abs(plate_geometry['z_outer'])
        a = crack_geometry['a']
        c = crack_geometry['c']
        shape = [(0, -t),
                 (W, -t),
                 (W, 0), ]

        for phi in np.pi*np.arange(0, 0.55, 0.05):
            x = c*np.cos(phi)
            y = -a*np.sin(phi)
            shape.append((x, y))
        I, _, _ = section.inertia(shape)
        xbar, ybar = section.centroid(shape)
        init_bc = -(I*Sys_ratio*Sys)/(ybar*t*W*(z_outer-z_inner))

    # modify_bc(init_bc, input_file, model_type)
    return init_bc


def get_upper_bd(input_file, crack_geometry, plate_geometry, material,
                 model_type):
    """Set the upper bound of boundary condition for optimization.
    Detect if the model is in tension or bending, and set the appropriate
    BC to level causing maximum defined strain.
    """
    E = material['E']
    Sys = material['Sys']
    n = material['n']
    L = plate_geometry['L']
    (strain, stress) = pre.lppl_unformatted(E, Sys, n)
    strain_ys = Sys/E
    strain_ys_ratio = (strain_ys+(np.max(strain)))/2
    if model_type == 'tension':
        # axial deflection required to cause remote strain of
        # strain_ys_ratio*Sys/E
        # u = (P*L)/(A*E) = (P/A)*(L/E) = sigma*(L/E)
        # u/L = strain = P/(A*E) = sigma/E
        upper_bd = strain_ys_ratio*L
    else:
        # TODO: solve this equation for bending strain
        # transverse traction along pin line required to cause remote stress
        # of Sys_ratio*Sys (Sys_ratio can be 0.5, 0.9, etc.)
        # Fundamental equations:
        # sigma = (M*t)/(2*I)
        # M(x) = P*(z_outer-z_inner)
        # P = traction*W*t
        #
        # become:
        # M = (2*I*sigma)/t
        # P*(z_outer-z_inner) = (2*I*sigma)/t
        # traction*W*t*(z_outer-z_inner) = (2*I*sigma)/t
        # traction = (2*I*sigma)/(t^2*W*(z_outer-z_inner))
        W = plate_geometry['W']
        t = plate_geometry['t']
        z_inner = np.abs(plate_geometry['z_inner'])
        z_outer = np.abs(plate_geometry['z_outer'])
        Iyy = W*np.power(t, 3)/12.0
        upper_bd = (2*Iyy*Sys)/(t*t*W*(z_outer-z_inner))

    # modify_bc(init_bc, input_file, model_type)
    return upper_bd


def get_model_type(input_file):
    """Detect if a WARP3D input file is for a tension model or a bending
    model. Returns 'tension' for tension models, 'bending' for bending
    models, and raises a ValueError otherwise.
    """
    with open(input_file) as f:
        input_lines = f.readlines()

    # If we're a tension case, will have non-zero w constraints on nodes.
    # - Find line starting with "constraints"
    # - Find next line starting with "c *echo off"
    # - In that range, find lines like
    #   "  10294  u 0.000000E+00  w 1.000000E+00"
    start_index = util.find_first_line_matching(
            input_lines, "constraints")+1
    end_index = util.find_first_line_matching(
            input_lines, "c *echo off", start_index)
    for i in range(start_index, end_index):
        line = input_lines[i]
        if 'w' in line:
            (rest, w_value) = line.rsplit(maxsplit=1)
            w_value = float(w_value)
            if not (w_value == 0):
                return 'tension'
    # If we're a bending case, will have ty tractions on elements.
    # - Find line starting with "loading set1"
    # - Find next line starting with "c"
    # - In that range, find lines like
    #   "   2571 face  2  traction ty    3.52000000E-01"
    start_index = util.find_first_line_matching(
            input_lines, "loading set1")+1
    end_index = util.find_first_line_matching(
            input_lines, "c", start_index)
    for i in range(start_index, end_index):
        line = input_lines[i]
        if 'ty' in line:
            (rest, ty_value) = line.rsplit(maxsplit=1)
            ty_value = float(ty_value)
            if not (ty_value == 0):
                return 'bending'

    raise ValueError("Could not detect model type")


def optimize_bc(input_file, crack_geometry, plate_geometry, material):
    """Optimize boundary condition for run_model until dimensionless M
    value is within range of [20, 25]. Returns final boundary condition
    value.
    """
    # Modify values for boundary conditions at each load step to ensure
    # model deforms well into the plastic regime, but not so far as to become
    # numerically unstable. Requires figuring out dimensionless value M
    # at phi=30 and phi=90 at the final load step for a given boundary
    # condition magnitude, then adjusting the magnitude until 20 < M < 25.
    # Can be run on any platform where the WARP3D input file exists, but in
    # a practical sense, will be run on an HPC system to avoid excessive
    # file transfers.

    model_type = get_model_type(input_file)
    initial_bc = get_initial_bc(input_file, crack_geometry, plate_geometry,
                                material, model_type)
#    upper_bd = get_upper_bd(input_file, crack_geometry, plate_geometry,
#                            material, model_type)
    if model_type == 'tension':
        bc_optimal = scipy.optimize.newton(
                run_model, x0=initial_bc, tol=9e-5,
                args=(input_file, crack_geometry, plate_geometry, material,
                      model_type))
    else:
        bc = initial_bc
        while run_model(bc, input_file, crack_geometry, plate_geometry,
                        material, model_type) > 0:
            bc = bc+0.1*initial_bc
        bc_optimal = bc
#    bc_optimal = scipy.optimize.minimize(
#            run_model, x0=np.array(initial_bc), bounds=((0, upper_bd),),
#            args=(input_file, crack_geometry, plate_geometry, material,
#                  model_type))
#    bc_optimal = scipy.optimize.minimize(
#            run_model, x0=np.array(initial_bc), bounds=((0.9*initial_bc,
#                                                         0.7),),
#            args=(input_file, crack_geometry, plate_geometry, material,
#                  model_type),
#            method='SLSQP', options={'eps': 0.01})
    return bc_optimal
