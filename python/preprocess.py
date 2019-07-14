# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 10:39:12 2017

@author: Renfro
"""

import numpy as np
import os
import shutil
import subprocess
import util


def calculate_geometry(template_filename, depth_ratio, aspect_ratio,
                       thickness):
    """Calculate plate and crack geometry values from base parameters.
    Returns a tuple of (a, c, W, t, L, S_inner, S_outer), where S_inner
    and S_outer are None for tension cases.
    """
    a = depth_ratio*thickness
    c = a/aspect_ratio
    if (5*c > 5*thickness):
        W = 5*c
    else:
        W = 5*thickness
    t = thickness
    model_type = get_elt_type(template_filename)
    if (model_type == 'bending'):
        # Need to define S_inner and S_outer
        S_inner = W  # full 2*W in Excel, but .elt uses symmetry
        S_outer = 2*S_inner
        L = np.max(np.array([2*W, 1.1*S_outer]))
    else:
        # Tension model. No S_inner, S_outer
        S_inner = None
        S_outer = None
        L = 2*W
    return (a, c, W, t, L, S_inner, S_outer)


def get_elt_type(elt_file):
    with open(elt_file) as f:
        elt_lines = f.readlines()

    # If we're a bending case, will have specified ty tractions on elements.
    # - Find line starting with "Notes:" (including quotes)
    # - On next line, look for "*use bottom load pin plate ty "
    #   (without quotes)
    start_index = util.find_first_line_matching(
            elt_lines, '"Notes:"')+1
    if not(start_index == -1):
        # Found the notes line
        notes_content = elt_lines[start_index]
        if "*use bottom load pin plate ty " in notes_content:
            elt_type = 'bending'
        else:
            elt_type = 'tension'
    else:
        elt_type = 'tension'

    # Extra sanity check: bending models should have two pin definitions.
    # May or may not be possible in the FEACrack GUI, but they have
    # to be present to locate the force and reaction roller locations.
    # First an outer roller definition:
    #  "RigidSurfaceData_Radius",.75
    #  "RigidSurfaceData_PinLocation",10
    # and then an inner roller definition:
    #  "RigidSurfaceData_Radius",.75
    #  "RigidSurfaceData_PinLocation",5
    if elt_type == 'bending':
        outer_roller_radius_index = util.find_first_line_containing(
                elt_lines, '"RigidSurfaceData_Radius",')
        outer_roller_location_index = util.find_first_line_containing(
                elt_lines, '"RigidSurfaceData_PinLocation",')
        if (outer_roller_radius_index ==
                -1) or (outer_roller_location_index == -1):
            elt_type = 'invalid'
            extra_info = "No rollers fully defined"
        else:
            inner_roller_radius_index = util.find_first_line_containing(
                    elt_lines, '"RigidSurfaceData_Radius",',
                    outer_roller_radius_index+1)
            inner_roller_location_index = util.find_first_line_containing(
                    elt_lines, '"RigidSurfaceData_PinLocation",',
                    outer_roller_location_index+1)
            if (inner_roller_radius_index ==
                    -1) or (inner_roller_location_index == -1):
                elt_type = 'invalid'
                extra_info = "Only one roller fully defined"

    if (elt_type == 'invalid'):
        raise ValueError("Malformed .elt file for bending: {0}".format(
                extra_info))
    else:
        return elt_type


def adjust_roller_locations(elt_filename, S_inner, S_outer):
    """Modifies the position of inner and outer rollers for FEACrack
    bending models. Must be done before building the WARP3D input file.
    """
    # FIXME: is the roller order (outer, inner) or (inner, outer)?
    # Comments say one thing, code says another
    with open(elt_filename) as f:
        elt_lines = f.readlines()

    outer_roller_location_index = util.find_first_line_containing(
            elt_lines, '"RigidSurfaceData_PinLocation",')
    if not (outer_roller_location_index == -1):
        # Adjust roller location to S_inner
        roller_line = elt_lines[outer_roller_location_index]
        (prefix, location) = roller_line.split(',')
        roller_line = "{0},{1}\n".format(prefix, S_outer)
        elt_lines[outer_roller_location_index] = roller_line
    else:
        raise ValueError("Couldn't find line for first roller location")

    inner_roller_location_index = util.find_first_line_containing(
            elt_lines, '"RigidSurfaceData_PinLocation",',
            outer_roller_location_index+1)
    if not (inner_roller_location_index == -1):
        # Adjust roller location to S_outer
        roller_line = elt_lines[inner_roller_location_index]
        (prefix, location) = roller_line.split(',')
        roller_line = "{0},{1}\n".format(prefix, S_inner)
        elt_lines[inner_roller_location_index] = roller_line
    else:
        raise ValueError("Couldn't find line for second roller location")

    with open(elt_filename, 'w') as f:
        f.writelines(elt_lines)


def get_elt_filename(global_elt_filename, a, c, L, W, t):
    """Given a global template filename like 'tens.elt' and a set of
    geometry parameters (a, c, L, W, t), return a filename like:
    'tens_ac1.0_at0.8_L10.00_W05.00.elt'
    """
    basename = os.path.splitext(global_elt_filename)[0]
    elt_filename = ("{0}_ac{1:.1f}_at{2:.1f}_"
                    "L{3:05.2f}_W{4:05.2f}.elt").format(
            basename, a/c, a/t, L, W)
    return elt_filename


def get_msh_filename(inp_filename):
    """Given a input filename like
    'tens_ac1.0_at0.8_L10.00_W05.00_E0500_n06_wrp.inp',
    return a filename like 'tens_ac1.0_at0.8_L10.00_W05.00_msh.out'
    """
    basename = os.path.splitext(inp_filename)[0]
    basename = basename[:-4]  # remove _wrp suffix
    (base, _, _) = basename.rsplit(sep='_', maxsplit=2)
    msh_filename = "{0}_msh.out".format(base)
    return msh_filename


def get_generic_inp_filename(elt_filename):
    """Given a FEACrack filename like
    'tens_ac1.0_at0.8_L10.00_W05.00.elt', return a filename like:
    'tens_ac1.0_at0.8_L10.00_W05.00_wrp.inp'
    """
    basename = os.path.splitext(elt_filename)[0]
    elt_filename = ("{0}_wrp.inp").format(basename)
    return elt_filename


def get_specific_inp_filename(generic_inp_filename, E, n):
    """Given a WARP3D input filename like
    'tens_ac1.0_at0.8_L10.00_W05.00_wrp.inp' and a set of material parameters
    (E, n), return a filename like:
    'tens_ac1.0_at0.8_L10.00_W05.00_E0500_n06_wrp.inp'
    """
    basename = os.path.splitext(generic_inp_filename)[0]
    basename = basename[:-4]  # remove _wrp suffix
    inp_filename = ("{0}_E{1:04d}_n{2:02d}_wrp.inp").format(
            basename, int(E), int(n))
    return inp_filename


def get_bpf_filename(inp_filename):
    """Given a WARP3D input filename like
    'tens_ac1.0_at0.8_L10.00_W05.00_E0500_n06_wrp.inp', return a filename
    like:
    'tens_ac1.0_at0.8_L10.00_W05.00_E0500_n06.bpf'
    """
    basename = os.path.splitext(inp_filename)[0]
    basename = basename[:-4]  # remove _wrp suffix
    bpf_filename = ("{0}.bpf").format(basename)
    return bpf_filename


def get_out_filename(inp_filename):
    """Given a input filename like
    'tens_ac1.0_at0.8_L10.00_W05.00_E0500_n06_wrp.inp',
    return a filename like 'tens_ac1.0_at0.8_L10.00_W05.00_E0500_n06_wrp.out'
    """
    basename = os.path.splitext(inp_filename)[0]
#    (base, _, _) = basename.rsplit(sep='_', maxsplit=2)
    out_filename = "{0}.out".format(basename)
    return out_filename


def build_mesh(global_elt_filename, a, c, L, W, t, S_inner, S_outer):
    """Modify values for a, c, L, W, t, S_inner, S_outer in the .elt file.
    a, c, L, W, t can be passed as command-line parameters to feacrack.exe.
    S_inner and S_outer have to be modified in the .elt file before
    running feacrack.exe with command-line parameters to build a new
    WARP3D input file. Has to be run on Windows, due to dependency on
    FEACrack executable.

    Returns the name of the specific WARP3D input file created.
    """
    feacrack_path = r'C:\Program Files (x86)\Quest Integrity Group\FEACrack'
    feacrack_exe = feacrack_path+r'\feacrack.exe'
    # Make copy of global .elt file to a new name, as that will control the
    # name of the WARP3D input file
    elt_filename = get_elt_filename(global_elt_filename, a, c, L, W, t)
    shutil.copy(global_elt_filename, elt_filename)
    # IF a bending model (S_inner and S_outer defined), adjust roller
    # positions in new .elt file to new values before creating mesh.
    if (not (S_inner is None)) and (not (S_outer is None)):
        adjust_roller_locations(elt_filename, S_inner, S_outer)
    # Run feacrack.exe with desired geometry parameters to build WARP3D file
    run_feacrack = True
    feacrack_cmd = ("feacrack.exe {0} /ModelWidth {1} /ModelThickness {2} "
                    "/ModelLength {3} /CrackLength {4} /CrackDepth {5} "
                    "/BuildMesh")
    print(feacrack_cmd.format(elt_filename, W, t, L, 2.0*c, a))
    if run_feacrack:
        subprocess.check_output([feacrack_exe,
                                 os.path.join(os.getcwd(), elt_filename),
                                 "/ModelWidth", str(W),
                                 "/ModelThickness", str(t),
                                 "/ModelLength", str(L),
                                 "/CrackLength", str(2.0*c),
                                 "/CrackDepth", str(a),
                                 "/BuildMesh",
                                 ],
                                stderr=subprocess.STDOUT)
    # inp_filename = elt_filename.rsplit(sep='.', maxsplit=1)[0]+"_wrp.inp"
    elt_filename = get_elt_filename(global_elt_filename, a, c, L, W, t)
    inp_filename = get_generic_inp_filename(elt_filename)
    return inp_filename


def lppl_unformatted(E, Sys, n):
    """Return arrays of stress strain data for a given elastic modulus,
    yield strength, and hardening exponent.
    """
    plastic_strain = 0.001*np.arange(1, 9)
    plastic_strain = np.append(plastic_strain,
                               plastic_strain[-1]+0.005*np.arange(1, 5))
    plastic_strain = np.append(plastic_strain,
                               plastic_strain[-1]+0.01*np.arange(1, 8))
    strain_ys = Sys/E
    strain_powerlaw = strain_ys+plastic_strain
    stress_powerlaw = Sys*np.power(strain_powerlaw/strain_ys, 1.0/n)
    return (strain_powerlaw, stress_powerlaw)


def lppl_formatted(E, Sys, n):
    """Return a formatted list of lines of stress strain data for a given
    elastic modulus, yield strength, and hardening exponent. Lines are
    in a format ready to replace an existing stress-strain curve inside
    a WARP3D input file.

    List has columns of total strain (elastic+plastic) and stress level.
    """
    (strain_powerlaw, stress_powerlaw) = lppl_unformatted(E, Sys, n)
    strain_ys = Sys/E
    table = np.column_stack((strain_powerlaw, stress_powerlaw))
    formatted_table_lines = ["{0:.6e} {1:.6e}\n".format(strain_ys, Sys), ]
    for row in table:
        formatted_table_lines.append("{0:.6e} {1:.6e},\n".format(
                row[0], row[1]))
    last_line = formatted_table_lines[-1]
    formatted_table_lines[-1] = last_line[:-2]+'\n'
    return formatted_table_lines


def stress_strain_formatted(stress_strain_table):
    pass


def change_material_properties(inp_filename,
                               E=None, Sys=None, n=None,
                               stress_strain_table=None):
    """Write a new WARP3D input file with specific material properties.
    Returns the name of the final WARP3D input file.
    """
    # Sanity check parameters
    if inp_filename is None:
        raise ValueError(('inp_filename must be the filename '
                          'of the .inp file for this specific geometry'))
    elif (E is None) or (Sys is None):
        raise ValueError('E and Sys must both be values')
    elif not((n is None) ^ (stress_strain_table is None)):
        raise ValueError(('1 and only 1 of n and stress_strain_table '
                          'must have a value'))
#    inp_filename = get_generic_inp_filename(global_elt_filename,
#                                            a, c, L, W, t)
    model_filename = get_specific_inp_filename(inp_filename, E, n)
    bpf_filename = get_bpf_filename(model_filename)
    shutil.copy(inp_filename, model_filename)
    with open(model_filename) as f:
        model_lines = f.readlines()
    # Find lines in the .inp file corresponding to material properties
    # Adjust material properties to match either LPPL or stress-strain table
    # values

    # Look for "stress-strain curve      1" and "c" to find stress-strain data
    start_index = util.find_first_line_matching(
            model_lines, "stress-strain curve      1")+1
    if start_index == 0:
        raise ValueError("Didn't find start of stress-strain curve")
    end_index = util.find_first_line_matching(
            model_lines, "c", start_index)
    if end_index == -1:
        raise ValueError("Didn't find end of stress-strain curve")
    # Decide on contents of replacement stress-strain data
    if n is None:
        # stress-strain table in use, just plug in those values instead
        stress_strain_table_formatted = stress_strain_formatted(
                stress_strain_table)
    else:
        # LPPL in use, calculate 20 entries for stress-strain table
        stress_strain_table_formatted = lppl_formatted(E, Sys, n)

    model_lines[start_index:end_index] = stress_strain_table_formatted

    # Look for "material mat1" to find elastic properties 3 lines below
    E_line_index = util.find_first_line_matching(model_lines,
                                                 "material mat1")+3
    (E_label, E_orig, nu_label, nu) = model_lines[E_line_index].split()
    new_e_line = "{0} {1:6e} {2} {3}\n".format(E_label, E, nu_label, nu)
    model_lines[E_line_index] = new_e_line

    # Fix comments that describe material properties, too.
    start_index = util.find_first_line_matching(model_lines,
                                                "c  Mesh Material Data")+1
    end_index = util.find_first_line_matching(model_lines,
                                              "c  Mesh Material Data",
                                              start_index)

    mat_notes = """c      number of material groups =         1
c                 material group =         1
c                  material type = total-strain vs stress table
c          modulus of elasticity =   {0}
c                  poisson ratio =   3.0000001E-01
c           number of table rows =        20
c                total-strain        stress
""".format(E)
    for line in stress_strain_table_formatted:
        mat_notes = mat_notes+"c {0}".format(line)
    model_lines[start_index:end_index] = mat_notes

    # Fix line setting which binary packet file (bpf) is written
    bpf_index = util.find_first_line_containing(model_lines,
                                                'binary packets on file')
    bpf_line = "  binary packets on file '{0}'\n".format(bpf_filename)
    model_lines[bpf_index] = bpf_line

    with open(model_filename, 'w') as f:
        f.writelines(model_lines)

    return model_filename


def add_rigid_plane(model_filename=None, stiffness=None, origin=None,
                    size=None):
    if model_filename is None:
        raise ValueError("Missing model_filename")
    if stiffness is None:
        raise ValueError("Missing rigid plane stiffness")
    if origin is None:
        raise ValueError("Missing rigid plane origin")
    if size is None:
        raise ValueError("Missing rigid plane size")
    plane_string = """    contact surface
      surface 1 plane
        point {0} {1} 0
        point {0} {2} 0
        point {3} {1} 0 $ relative normal vector ends up as (0, 0, -1)
        stiffness {4}
""".format(origin[0], origin[1],
           origin[1]+size[1],
           origin[0]+size[0],
           stiffness)
    # print(plane_string)
    with open(model_filename) as f:
        model_lines = f.readlines()
    start_index = util.find_first_line_matching(model_lines,
                                                "c          ANALYSIS PARAMETERS")-2
    model_lines.insert(start_index, plane_string)
    with open(model_filename, 'w') as f:
        f.writelines(model_lines)


def enable_t_stresses(model_filename=None):
    if model_filename is None:
        raise ValueError("Missing model_filename")
    with open(model_filename) as f:
        model = f.read()
    model = model.replace('output packets J', 'output packets I J')
    model = model.replace('compute domain integral',
                          'compute domain interaction integral')
    with open(model_filename, 'w') as f:
        model = f.write(model)


def get_filenames(global_elt_filename, crack_geometry, plate_geometry,
                  material):
    """Given a global template filename like 'bend.elt' and a set of
    geometry and material parameters, return a dictionary of filenames for
    this model.
    """
    filenames = {}
    a = float(crack_geometry['a'])
    c = float(crack_geometry['c'])
    L = float(plate_geometry['L'])
    W = float(plate_geometry['W'])
    t = float(plate_geometry['t'])
    E = int(material['E'])
    n = int(material['n'])

    filenames['global_elt'] = global_elt_filename  # bend.elt

    basename = os.path.splitext(global_elt_filename)[0]  # bend
    geometry_format = "{0}_ac{1:.1f}_at{2:.1f}_L{3:05.2f}_W{4:05.2f}"
    geometry_material_format = geometry_format+"_E{5:04d}_n{6:02d}"

    elt = (geometry_format+".elt").format(basename, a/c, a/t, L, W)
    filenames['elt'] = elt  # bend_ac0.8_at0.2_L11.00_W05.00.elt

    msh = (geometry_format+"_msh.out").format(basename, a/c, a/t, L, W)
    filenames['msh'] = msh  # bend_ac0.8_at0.2_L11.00_W05.00_msh.out

    generic_inp = (geometry_format+"_wrp.inp").format(basename, a/c, a/t, L, W)
    # bend_ac0.8_at0.2_L11.00_W05.00_wrp.inp
    filenames['generic_inp'] = generic_inp

    specific_inp = (geometry_material_format+"_wrp.inp").format(basename, a/c,
                                                                a/t, L, W, E,
                                                                n)
    # bend_ac0.8_at0.2_L11.00_W05.00_E0100_n10_wrp.inp
    filenames['specific_inp'] = specific_inp

    bpf = (geometry_material_format+".bpf").format(basename, a/c, a/t, L, W,
                                                   E, n)
    filenames['bpf'] = bpf  # bend_ac0.8_at0.2_L11.00_W05.00_E0100_n10.bpf

    out = (geometry_material_format+"_wrp.out").format(basename, a/c, a/t, L,
                                                       W, E, n)
    filenames['out'] = out  # bend_ac0.8_at0.2_L11.00_W05.00_E0100_n10_wrp.out

    return filenames
