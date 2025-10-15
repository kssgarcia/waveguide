# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 16:21:52 2025

@author: agarz
"""


def fun_parent(x, a=1):
    
    def fun_child(x):
        return a + x
    
    # a = 1
    return fun_child(x)

a = 10
b = fun_parent(1, 3)

print(b)
    
    
