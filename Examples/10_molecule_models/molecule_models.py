#!/usr/bin/env python

import httk.httkio
import httk.atomistic.vis
from httk.atomistic import Structure
from httk.core.vectors.fracvector import FracVector



diamond_host_struct = httk.load("input/Diamond.cif")
diamond_defect_struct = httk.load("input/Diamond_NV_center.cif")
diamond_defect_coords = [FracVector((15, 15, 19), 32), FracVector((17, 17, 21), 32)]

CaO_host_struct = httk.load("input/Calcium_oxide.cif")
CaO_defect_struct = httk.load("input/BiV_in_CaO.cif")
CaO_defect_coords = [FracVector((3, 3, 4), 8), FracVector((4, 3, 4), 8)]


molecule = diamond_host_struct.to_molecule_no_defect(layers=4, structure_type = 'Diamond', ref=FracVector((1, 1, 1), 2), termination_group='H')
molecule.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
molecule.vis.wait()

defect, host = diamond_defect_struct.to_molecule(host_struct = diamond_host_struct, layers = 3, defect_type=  "NV", defect_coords = diamond_defect_coords, termination_group = "H")
defect.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
host.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
host.vis.wait()

defect, host = CaO_defect_struct.to_molecule(host_struct = CaO_host_struct, layers = 2, defect_type=  "BiV", defect_coords = CaO_defect_coords, termination_group = "H")
defect.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
host.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
host.vis.wait()
