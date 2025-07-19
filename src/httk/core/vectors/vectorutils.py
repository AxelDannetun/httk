#!/usr/bin/env python

import math
from httk.core.vectors.fracvector import FracVector

#TODO: Might want to generalize to FracVectors of any size
#TODO: Might want to change "Cube" to sup-norm 
def norm(v, ord="L2"):
    tmp = v.to_floats()
    match ord:
        case "Cube":
            return max(abs(tmp[0]), abs(tmp[1]), abs(tmp[2]))
        case "L1":
            return abs(tmp[0]) + abs(tmp[1]) + abs(tmp[2])
        case "L2":
            return math.sqrt(tmp[0]**2 + tmp[1]**2 + tmp[2]**2) 
        case _:
            raise Exception('Norm "' + ord + '" is not supported')


#TODO: Might want L1 dist to be a cube instead of diamond shape
def dist(a, b, ord = "L2"):
    diff = a - b
    return norm(diff, ord=ord)    