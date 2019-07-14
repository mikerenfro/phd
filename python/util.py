# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 10:36:33 2017

@author: Renfro
"""


def find_first_line_matching(lines, text, start=0):
    """Find first line matching some text"""
    for i in range(start, len(lines)):
        line = lines[i].strip()
        if line == text:
            return i
    return -1


def find_last_line_matching(lines, text, end):
    """Find last line matching some text"""
    for i in range(end, 0, -1):
        line = lines[i].strip()
        if line == text:
            return i
    return -1


def find_first_line_containing(lines, text, start=0):
    """Find first line containing some text"""
    for i in range(start, len(lines)):
        line = lines[i].strip()
        if text in line:
            return i
    return -1


def find_last_line_containing(lines, text, end):
    """Find last line containing some text"""
    for i in range(end, 0, -1):
        line = lines[i].strip()
        if text in line:
            return i
    return -1
