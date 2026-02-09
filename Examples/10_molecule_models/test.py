#!/usr/bin/env python

import httk.httkio
import httk.atomistic.vis
from httk.atomistic import Structure
import httk.atomistic.structure
from httk.core.vectors.fracvector import FracVector

from httk.core.vectors.vectorutils import *

diamond_host_struct = httk.load("input/Diamond.cif")
diamond_defect_struct = httk.load("input/Diamond_NV_center.cif")

CaO_host_struct = httk.load("input/Calcium_oxide.cif")
CaO_defect_struct = httk.load("input/BiV_in_CaO.cif")

print(diamond_host_struct.assignments)
print(diamond_defect_struct.assignments)

print(diamond_defect_struct.uc)

diamond_defect_coords = diamond_defect_struct.uc_reduced_coordgroups[1][0]

defect, host = diamond_defect_struct.to_molecule(host_struct = diamond_host_struct, layers = 3, defect_type=  "NV", defect_coords = diamond_defect_coords, termination_group = "H")
defect.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
#host.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
#host.vis.wait()

print(defect.assignments)
print(defect.uc_reduced_coordgroups.to_floats()[0])
print(defect.uc_reduced_coordgroups.to_floats()[1])
print(defect.uc_reduced_coordgroups.to_floats()[2])
print(defect.uc_reduced_coordgroups.to_floats()[3])

CaO_host_struct = httk.load("input/Calcium_oxide.cif")
CaO_defect_struct = httk.load("input/BiV_in_CaO.cif")

