#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:16:44 2018

@author: renfro
"""

# Ripped from Dr. Drang
# (http://leancrew.com/all-this/2018/01/python-module-for-section-properties/)

#import matplotlib.pyplot as plt
from math import atan2, sin, cos, sqrt, pi, degrees


def area(pts):
    'Area of cross-section.'

    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [c[0] for c in pts]
    y = [c[1] for c in pts]
    s = 0
    for i in range(len(pts) - 1):
        s += x[i]*y[i+1] - x[i+1]*y[i]
    return s/2


def centroid(pts):
    'Location of centroid.'

    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [c[0] for c in pts]
    y = [c[1] for c in pts]
    sx = sy = 0
    a = area(pts)
    for i in range(len(pts) - 1):
        sx += (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
        sy += (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
    return sx/(6*a), sy/(6*a)


def inertia(pts):
    'Moments and product of inertia about centroid.'

    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [c[0] for c in pts]
    y = [c[1] for c in pts]
    sxx = syy = sxy = 0
    a = area(pts)
    cx, cy = centroid(pts)
    for i in range(len(pts) - 1):
        sxx += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        syy += (x[i]**2 + x[i]*x[i+1] + x[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        sxy += (x[i]*y[i+1] + 2*x[i]*y[i] +
                2*x[i+1]*y[i+1] + x[i+1]*y[i])*(x[i]*y[i+1] - x[i+1]*y[i])
    return sxx/12 - a*cy**2, syy/12 - a*cx**2, sxy/24 - a*cx*cy


def principal(Ixx, Iyy, Ixy):
    'Principal moments of inertia and orientation.'

    avg = (Ixx + Iyy)/2
    diff = (Ixx - Iyy)/2      # signed
    I1 = avg + sqrt(diff**2 + Ixy**2)
    I2 = avg - sqrt(diff**2 + Ixy**2)
    theta = atan2(-Ixy, diff)/2
    return I1, I2, theta


def summary(pts):
    'Text summary of cross-sectional properties.'

    a = area(pts)
    cx, cy = centroid(pts)
    Ixx, Iyy, Ixy = inertia(pts)
    I1, I2, theta = principal(Ixx, Iyy, Ixy)
    summ = """Area
  A = {}
Centroid
  cx = {}
  cy = {}
Moments and product of inertia
  Ixx = {}
  Iyy = {}
  Ixy = {}
Principal moments of inertia and direction
  I1 = {}
  I2 = {}
  θ︎ = {}°""".format(a, cx, cy, Ixx, Iyy, Ixy, I1, I2, degrees(theta))
    return summ


# def outline(pts, basename='section', format='pdf', size=(8, 8), dpi=100):
#     'Draw an outline of the cross-section with centroid and principal axes.'
#
#     if pts[0] != pts[-1]:
#         pts = pts + pts[:1]
#     x = [c[0] for c in pts]
#     y = [c[1] for c in pts]
#
#     # Get the bounds of the cross-section
#     minx = min(x)
#     maxx = max(x)
#     miny = min(y)
#     maxy = max(y)
#
#     # Whitespace border is 5% of the larger dimension
#     b = .05*max(maxx - minx, maxy - miny)
#
#     # Get the properties needed for the centroid and principal axes
#     cx, cy = centroid(pts)
#     i = inertia(pts)
#     p = principal(*i)
#
#     # Principal axes extend 10% of the minimum dimension from the centroid
#     length = min(maxx-minx, maxy-miny)/10
#     a1x = [cx - length*cos(p[2]), cx + length*cos(p[2])]
#     a1y = [cy - length*sin(p[2]), cy + length*sin(p[2])]
#     a2x = [cx - length*cos(p[2] + pi/2), cx + length*cos(p[2] + pi/2)]
#     a2y = [cy - length*sin(p[2] + pi/2), cy + length*sin(p[2] + pi/2)]
#
#     # Plot and save
#     # Axis colors chosen from http://mkweb.bcgsc.ca/colorblind/
#     fig, ax = plt.subplots(figsize=size)
#     ax.plot(x, y, 'k*-', lw=2)
#     ax.plot(a1x, a1y, '-', color='#0072B2', lw=2)     # blue
#     ax.plot(a2x, a2y, '-', color='#D55E00')           # vermillion
#     ax.plot(cx, cy, 'ko', mec='k')
#     ax.set_aspect('equal')
#     plt.xlim(xmin=minx-b, xmax=maxx+b)
#     plt.ylim(ymin=miny-b, ymax=maxy+b)
#     filename = basename + '.' + format
#     plt.savefig(filename, format=format, dpi=dpi)
#     plt.close()
