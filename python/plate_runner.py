#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plate_runner.py

Created on Fri Dec 22 17:07:36 2017

@author: renfro
"""

import postprocess as post
import preprocess as pre
import solve
import numpy as np
import pandas as pd
import os

# from v2014_config import aspect_ratio_list, depth_ratio_list
# from v2014_config import plate_thickness, E_Sys_ratio_list
# from v2014_config import Sys, n_list
# from v2014_config import elt_global_template_filename
from minimal_config import aspect_ratio_list, depth_ratio_list
from minimal_config import plate_thickness, E_Sys_ratio_list
from minimal_config import Sys, n_list
from minimal_config import elt_global_template_filename
# from bend_config_all import aspect_ratio_list, depth_ratio_list
# from bend_config_all import plate_thickness, E_Sys_ratio_list
# from bend_config_all import Sys, n_list
# from bend_config_all import elt_global_template_filename
# from validate_experiment_config import aspect_ratio_list, depth_ratio_list
# from validate_experiment_config import plate_thickness, E_Sys_ratio_list
# from validate_experiment_config import Sys, n_list
# from validate_experiment_config import elt_global_template_filename

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def do_postprocessing(input_file):
    bpf_file = pre.get_bpf_filename(input_file)
    print("Postprocessing: {0}".format(bpf_file))
    results_files = post.extract_results(bpf_file)
    nfile = pre.get_msh_filename(input_file)
#    dtable = post.get_nodal_results(dfile)
    rtable = post.get_nodal_results(results_files['reactions'])
    ntable = post.get_node_coordinates(nfile)
    output_file = pre.get_out_filename(input_file)
    msh_file = pre.get_msh_filename(input_file)
    with open(output_file) as f:
        output_lines = f.readlines()
    with open(input_file) as f:
        input_lines = f.readlines()
    coordinates = post.get_node_coordinates(msh_file)
    cmod = post.get_cmod(coordinates, results_files['displacements'],
                         x=0, y=0, z=0)
    (phi, J_table) = post.find_J(input_lines, output_lines)
    J_30 = post.J_interpolated(phi, J_table, np.pi/6)
    basename = os.path.splitext(input_file)[0]

    plt.figure()
    plt.plot(cmod, J_table[:, 0], ':', label=r'$\phi=0$')
    plt.plot(cmod, J_30, 'x-', label=r'$\phi=30$')
    plt.plot(cmod, J_table[:, -1], '--', label=r'$\phi=90$')
    plt.xlabel('CMOD')
    plt.ylabel(r'$J(\phi)$')
    plt.legend()
    plt.title(basename)
    plt.savefig("J_CMOD_"+basename+".pdf")
    plt.close()

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # for step in range(J_table.shape[0]):
    #     ax.plot(180*phi/np.pi, J_table[step, :], cmod[step+1])
    # ax.set_xlabel(r'$\phi$ (degrees)')
    # ax.set_ylabel('J')
    # plt.yticks(np.arange(0, 0.025, 0.005))
    # ax.set_zlabel('CMOD', rotation=90)
    # plt.title(basename)
    # plt.xticks(np.arange(0, 105, 15))
    # plt.savefig("J_phi_CMOD_"+basename+".pdf")
    # plt.close()

    if solve.get_model_type(input_file) == 'bending':
        z_outer = -10.0
        roller_reactions = post.get_roller_reactions(ntable, rtable,
                                                     z=z_outer)
        for step in roller_reactions['step'].unique():
            reactions = roller_reactions[
                    roller_reactions['step'] == step]
            nodes = pd.merge(ntable, reactions, how='inner', on=['node'])
            nodes.columns = ['node', 'x', 'y', 'z', 'step', 'Rx', 'Ry', 'Rz']
            if (nodes[nodes['Ry'] > 0].shape[0]) > 0:
                print(("On step {0}, at least one node "
                       "lifted off").format(step))

if __name__ == "__main__":
    for depth_ratio in depth_ratio_list:
        for aspect_ratio in aspect_ratio_list:
            (a, c, W, t, L,
             S_inner, S_outer) = pre.calculate_geometry(
                 elt_global_template_filename,
                 depth_ratio, aspect_ratio,
                 plate_thickness)

            for n in n_list:
                for E_Sys_ratio in E_Sys_ratio_list:
                    E = Sys*E_Sys_ratio
                    elt_file = pre.get_elt_filename(elt_global_template_filename,
                                                    a, c, L, W, t)
                    generic_file = pre.get_generic_inp_filename(elt_file)
                    input_file = pre.get_specific_inp_filename(generic_file,
                                                               E, n)
                    crack_geometry = {'a': a, 'c': c}
                    plate_geometry = {'L': L, 'W': W, 't': t}
                    if solve.get_model_type(input_file) == 'bending':
                        # Determine z_inner and z_outer from pin locations
                        (z_inner, z_outer) = solve.find_pin_locations(input_file)
                        plate_geometry['z_inner'] = z_inner
                        plate_geometry['z_outer'] = z_outer
                        material = {'E': E, 'Sys': Sys, 'n': n}

                        print("Running {0}".format(input_file))
                        solve.optimize_bc(input_file, crack_geometry, plate_geometry,
                                          material)

                        do_postprocessing(input_file)
