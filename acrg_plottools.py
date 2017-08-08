#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 12:13:25 2017

@author: as13988
"""

# Generates N HSV or RGB tuples spread over colour space 
class generatecolours:
    def __init__(self, N):
        
        import colorsys        
        
        HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
        RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
        
        self.RGB = RGB_tuples
        self.HSV = HSV_tuples