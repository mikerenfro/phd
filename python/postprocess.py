# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 10:49:17 2017

@author: Renfro
"""

import numpy as np
import scipy.interpolate
import os
from subprocess import Popen, PIPE
import util
import pandas as pd


# Extract J values from WARP3D model output. Used in calculating M,
# also used in TASC itself. Can be run on any platform where the WARP3D
# output exists (but most likely on an HPC system, since M is a function
# of J).
def find_J(input_lines, output_lines):
    """Calculate J values from the contents of a WARP3D input and output file.
    J values are calculated at each position around the crack front for each
    load step.
    """
    start_index = util.find_first_line_matching(input_lines,
                                                'c  Analysis Load Step Data')
    (_, steps) = input_lines[start_index+1].rsplit(maxsplit=1)
    # print("Steps: {0}".format(steps))
    steps = int(steps)

    start_index = util.find_first_line_matching(input_lines,
                                                'c  Crack Node Data')
    (_, crack_nodes) = input_lines[start_index+14].rsplit(maxsplit=1)
    # print("Crack nodes: {0}".format(crack_nodes))
    crack_nodes = int(crack_nodes)

    # Create list of phi and J values for the corner nodes and
    # number of load steps
    phi_table = np.zeros((int((crack_nodes+1)/2),))
    J_table = np.zeros((steps, int((crack_nodes+1)/2)))

    j = 0
    label_list = []
    for i in range(start_index+16, start_index+16+crack_nodes, 2):
        line = input_lines[i]
        # print(line)
        (c, node_id, phi, label, _) = line.split(maxsplit=4)
        phi_table[j] = float(phi)
        label_list.append(label)
        j = j+1

    start_index = 0
    for step in range(steps):
        start_index = util.find_first_line_containing(output_lines,
                  "loading: history1     step:", # analysis:ignore
                  start=start_index) # analysis:ignore
        for label in label_list:
            J_column = label_list.index(label)
            J_row = step
            # Need this line after adding interaction integrals (T stress)q
            start_index = util.find_first_line_containing(output_lines,
                      "J-integral components", start=start_index)
            start_index = util.find_first_line_containing(output_lines,
                      "average      minimum      maximum", # analysis:ignore
                      start=start_index+1) # analysis:ignore
            (_, _, J) = output_lines[start_index+1].split(maxsplit=2)
            J_table[J_row, J_column] = float(J)

    return (phi_table, J_table)


def find_J_list(input_lines, output_lines):
    """Calculate J values from the contents of a WARP3D input and output file.
    J values are calculated at each position around the crack front for each
    load step.
    """
    start_index = util.find_first_line_matching(input_lines,
                                                'c  Analysis Load Step Data')
    (_, steps) = input_lines[start_index+1].rsplit(maxsplit=1)
    # print("Steps: {0}".format(steps))
    steps = int(steps)

    start_index = util.find_first_line_matching(input_lines,
                                                'c  Crack Node Data')
    (_, crack_nodes) = input_lines[start_index+14].rsplit(maxsplit=1)
    # print("Crack nodes: {0}".format(crack_nodes))
    crack_nodes = int(crack_nodes)

    # Create list of phi and J values for the corner nodes and
    # number of load steps
    phi_table = np.zeros((int((crack_nodes+1)/2),))
    J_table = np.zeros((steps, int((crack_nodes+1)/2), 10))

    j = 0
    label_list = []
    for i in range(start_index+16, start_index+16+crack_nodes, 2):
        line = input_lines[i]
        # print(line)
        (c, node_id, phi, label, _) = line.split(maxsplit=4)
        phi_table[j] = float(phi)
        label_list.append(label)
        j = j+1

    start_index = 0
    for step in range(steps):
        # print("step: "+str(step))
        start_index = util.find_first_line_containing(output_lines,
                  "loading: history1     step:", # analysis:ignore
                  start=start_index) # analysis:ignore
        for label in label_list:
            start_index = util.find_first_line_containing(output_lines,
                                                          "domain id: "+label,
                                                          start=start_index)
            # print("label: "+label)
            J_column = label_list.index(label)
            J_row = step
            for i in range(10):
                start_index = util.find_first_line_containing(output_lines,
                          "killed ele", # analysis:ignore
                          start=start_index+1) # analysis:ignore
                # print("line: "+output_lines[start_index+1])
                J = output_lines[start_index+1].split()[10]
                J_table[J_row, J_column, i] = float(J)

    return (phi_table, J_table)


def J_interpolated(phi, J_table, phi_0):
    """Return an array of interpolated J-CMOD values at a specified angle.
    """
    J_new = []
    for index in range(J_table.shape[0]):
        J_interpolant = scipy.interpolate.interp1d(phi, J_table[index, :])
        J_new.append(J_interpolant(phi_0))
    return np.array(J_new)


def extract_results(bpf_file, run_packet_reader=True):
    """Use the packet_reader utility included with WARP3D to extract the
    following results from a binary packet file (.bpf) for all time steps:
    - nodal displacements
    - nodal reactions (numerically zero everywhere except at BCs)

    Each result goes into its own text file, named according to the basename
    of the .bpf file and the result type (basename-type.out).

    Returns a dict of output filenames, keyed by the result type
    (displacements, reactions).
    """
    # Extract reactions at top roller to calculate P
    # Extract positions of traction surface and roller to calculate
    # a length.
    # P=0 at rollers, so Moment = P*Distance at traction surface
    # and Szz = (M*t)/(2*I)
    results_map = [
            ('displacements', 1),
            ('reactions', 4),
            # ('stresses', 15),
            # ('j', 17),
            ]
    bpf_basename = os.path.splitext(bpf_file)[0]
    input_template = "{0}\n{1}\nn\n{2}\n"  # BPF file, code, output file

    filenames = {}
    for (kind, code) in results_map:
        output_filename = "{0}-{1}.out".format(
                bpf_basename,
                kind
                )
        input_bytes = bytearray(input_template.format(
                bpf_file,
                code,
                output_filename), 'utf8')
#        print("Writing {0} to {1}".format(
#                kind,
#                output_filename))
        if run_packet_reader:
            p = Popen(['packet_reader'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            (p_stdout, p_stderr) = p.communicate(input=input_bytes)
        filenames[kind] = output_filename
    return filenames


def get_node_coordinates(mesh_name):
    """Returns a Dataframe of node numbers and (x, y, z) coordinates from
    a FEACrack _msh.out file.
    """
    # Lines 4-5 of _msh.out contain:
    #
    # 3-D Mesh: #nodes, #elements, #nodes/element, #face elements, #face nodes
    #    43330     9324  20   3538  11457
    mesh_stats = pd.read_table(mesh_name,
                               sep='\s+',
                               skiprows=4,
                               nrows=1,
                               names=['nodes',
                                      'elements',
                                      'nodes_per_element',
                                      'face_elements',
                                      'face_nodes'])
    n_nodes = mesh_stats.iloc[0]['nodes']
    # Actual coordinates start on line 11, and now we know how many nodes
    # are in the mesh
    coordinates = pd.read_table(mesh_name,
                                sep='\s+',
                                skiprows=10,
                                nrows=n_nodes,
                                names=['node', 'x', 'y', 'z'],
                                )
    return coordinates


# Extract reaction forces and nodal displacements from WARP3D model
# output. Used in calculating P-CMOD curves for TASC. Can be run on any
# platform where the WARP3D output exists, (but most likely on an HPC
# system, since M is a function of P).
def get_displacements(displacement_name):
    d = pd.read_table(displacement_name, sep='\s+')
    d_final = None
    for step in np.unique(d['step']):
        # Extract displacements for this step: copying to get around
        # SettingWithCopyWarning, from "How to deal with SettingWithCopyWarning
        # in Pandas?" https://stackoverflow.com/a/20644369
        d_step = d[d['step'] == step].copy()
        # node column can't store values above 99999, so renumber correctly
        d_step['node'] = d_step.index.values-np.min(d_step.index.values)+1
        d_final = pd.concat((d_final, d_step))
    return d_final


def get_cmod(coordinates, displacement_name, x=0, y=0, z=0):
    """Returns a Dataframe of CMOD displacements, by doubling the nodal
    displacement of a specified node (by default, the node at
    (x, y, z) = (0, 0, 0)). Requires a Dataframe of node numbers and
    coordinates(typically returned by get_node_coordinates()), and a file
    containing the displacement output (typically provided by
    extract_results()).
    """
    # CMOD node is assumed to be at (0, 0, 0), and will normally be node 1.
    cmod_node = coordinates[(coordinates['x'] == x) &
                            (coordinates['y'] == y) &
                            (coordinates['z'] == z)]
    # https://stackoverflow.com/questions/33282119/
    displacements = get_displacements(displacement_name)
    keys = ['node']
    i1 = displacements.set_index(keys).index
    i2 = cmod_node.set_index(keys).index
    cmod_displacement = displacements[i1.isin(i2)]

    # cmod_number = cmod_node.iloc[0]['node']
    # cmod_data = displacements.loc[displacements['node'] == cmod_number]
    return np.abs(2*cmod_displacement['z-disp'])


def get_reaction(coordinates, reactions_name, z=0):
    """Return a list of total reaction force on the plate for all time steps.
    Requires a Dataframe of node numbers and coordinates (typically returned
    by get_node_coordinates()), and a file containing the reaction output
    (typically provided by extract_results()). By default, all nodes along
    the plane z=0 will be used, but other z values can be used instead.
    """
    reactions = pd.read_table(reactions_name,
                              sep='\s+',
                              skiprows=1,
                              names=['step', 'node', 'Rx', 'Ry', 'Rz'],
                              dtype={'step': np.int32,
                                     'node': np.float64,
                                     'Rx': np.float64,
                                     'Ry': np.float64,
                                     'Rz': np.float64},
                              na_values='TOT')
    # Node column isn't really a float, but can't treat 'TOT' as 'N/A' on an
    # integer column

    # Find all nodes on designated plane along length of plate
    reaction_nodes = coordinates[coordinates['z'] == z]
    # https://stackoverflow.com/questions/33282119/
    keys = ['node']
    i1 = reactions.set_index(keys).index
    i2 = reaction_nodes.set_index(keys).index
    actual_reactions = reactions[i1.isin(i2)]
    P_list = []
    for step in actual_reactions.step.unique():
        P = -np.sum(actual_reactions[actual_reactions['step'] == step]['Rz'])
        P_list.append(P)

    return P_list


def get_nodal_results(filename):
    table = pd.read_table(filename, sep='\s+', na_values='TOT')
    table = table.dropna()
    table_final = None
    for step in np.unique(table['step']):
        # Extract displacements for this step: copying to get around
        # SettingWithCopyWarning, from "How to deal with
        # SettingWithCopyWarning in Pandas?"
        # https://stackoverflow.com/a/20644369
        t_step = table[table['step'] == step].copy()
        # node column can't store values above 99999, so renumber correctly
        t_step['node'] = t_step.index.values-np.min(t_step.index.values)+1
        table_final = pd.concat((table_final, t_step))
    return table_final


def get_roller_reactions(coordinates, reactions, y=0, z=0):
    roller_nodes = coordinates[(coordinates['y'] == y) &
                               (coordinates['z'] == z)]
    keys = ['node']
    i1 = reactions.set_index(keys).index
    i2 = roller_nodes.set_index(keys).index
    roller_reactions = reactions[i1.isin(i2)]
    return roller_reactions


def get_line_displacements(coordinates, displacements, y=0, z=0):
    roller_nodes = coordinates[(coordinates['y'] == y) &
                               (coordinates['z'] == z)]
    keys = ['node']
    i1 = displacements.set_index(keys).index
    i2 = roller_nodes.set_index(keys).index
    line_displacements = displacements[i1.isin(i2)]
    return line_displacements
