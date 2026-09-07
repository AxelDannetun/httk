
#!/usr/bin/env python
from __future__ import print_function
from unittest import case

import numpy as np

import httk
import httk.atomistic.vis
from httk.atomistic import Structure, StructureTag
import httk.task
from httk.core.vectors.fracvector import FracVector


from pprint import pprint


setting = "cao"

match setting:
    case "diamond":
        defect_type = 'NV'
        host = httk.load('input/Diamond.cif')
    case "cao":
        defect_type = 'BiV'
        host = httk.load('input/Calcium_oxide.cif')
    case "tricky":
        defect_type = 'Tricky'
        host = httk.load('input/Silicon_carbide.cif')
    case _:
        raise ValueError(f"Unknown setting: {setting}")


print(host.assignments.symbols)
molecule_structure = host.to_molecule_no_defect(layers=2, structure_type=defect_type, termination_group='H')
molecule_structure.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})

print(molecule_structure.nested_layers['assignments'])
print('Core layer')
print(molecule_structure.get_assignments_layer('core_layer'))
print('Termination layer')
print(molecule_structure.get_assignments_layer('termination_layer'))