# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 13:36:26 2017

@author: Renfro
"""

# Full set of model parameters
plate_thickness = .374
depth_ratio_list = [.215/plate_thickness]
aspect_ratio_list = [(plate_thickness*depth_ratio_list[0])/(.490/2)]
Sys = 56.0e3
E_Sys_ratio_list = [10.8e6/Sys, ]
n_list = [9, ]
elt_global_template_filename = 'bend.elt'
