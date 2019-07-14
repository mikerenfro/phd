#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plate_creator.py

Created on Fri Dec 22 17:07:36 2017

@author: renfro
"""

# import numpy as np
import preprocess as pre

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

# modelIndex = 1
for depth_ratio in depth_ratio_list:
    for aspect_ratio in aspect_ratio_list:
        (a, c, W, t, L,
         S_inner, S_outer) = pre.calculate_geometry(
                     elt_global_template_filename,
                     depth_ratio, aspect_ratio,
                     plate_thickness)

        # Have all geometry parameters now. Call FEACrack with arguments
        # to write new WARP3D template file with this geometry
        geom_filename = pre.build_mesh(
                elt_global_template_filename, a, c, L, W, t,
                S_inner, S_outer)
        print("Wrote generic mesh {0}".format(geom_filename))
        for n in n_list:
            for E_Sys_ratio in E_Sys_ratio_list:
                E = Sys*E_Sys_ratio
                # Have all material parameters now. Write new .inp file with
                # those parameters
                model_filename = pre.change_material_properties(
                        inp_filename=geom_filename,
                        E=E, Sys=Sys, n=n)
                # Have all model parameters now.
                print("Wrote {0}".format(model_filename))
                # modelIndex = modelIndex+1
                pre.add_rigid_plane(model_filename=model_filename,
                                    stiffness=100.0*E, origin=(0, -t),
                                    size=(W, 11*t))
                print("Added rigid plane to {0}".format(model_filename))
                pre.enable_t_stresses(model_filename=model_filename)
                print("Added T stresses to {0}".format(model_filename))
