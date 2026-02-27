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

#diamond_host_struct.vis.show()

defect, host = diamond_defect_struct.to_molecule(host_struct = diamond_host_struct, layers = 3, defect_type=  "NV", defect_coords = diamond_defect_coords, termination_group = "H")
#defect.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
#host.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
#host.vis.wait()

print(defect.nested_layers['assignments'])
print('Core layer')
print(defect.get_assignments_layer('core_layer'))
print('Host layer')
print(defect.get_assignments_layer('host_layer'))
print('Termination layer')
print(defect.get_assignments_layer('termination_layer'))

CaO_host_struct = httk.load("input/Calcium_oxide.cif")
CaO_defect_struct = httk.load("input/BiV_in_CaO.cif")

