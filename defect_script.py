#!/usr/bin/env python
from __future__ import print_function
from unittest import case

import numpy as np

import httk
import httk.db
import httk.atomistic.vis
from httk.atomistic import Structure, StructureTag
import httk.task
from httk.atomistic.moleculeutils import *


from pprint import pprint

from classes import *


setting = "tricky"  # "diamond", "cao", "tricky"

match setting:
    case "diamond":
        backend = httk.db.backend.Sqlite('diamond/defects.sqlite')
        defect_type = 'NV'
        #host = httk.load('input/Diamond.cif')
    case "cao":
        backend = httk.db.backend.Sqlite('cao/defects.sqlite')
        defect_type = 'BiV'
        #host = httk.load('input/Calcium_oxide.cif')
    case "tricky":
        backend = httk.db.backend.Sqlite('sic/defects.sqlite')
        defect_type = 'Tricky'
        #host = httk.load('input/Silicon_carbide.cif')
    case _:
        raise ValueError(f"Unknown setting: {setting}")

#backend = httk.db.backend.Sqlite('httk/sic/defects.sqlite')
store = httk.db.store.SqlStore(backend)

# search
search = store.searcher()

search_host = search.variable(HostSuperCell)

search.output(search_host, 'hostsupercell')

for match in search:
    host = match[0][0].host_supercell
    

store = httk.db.store.SqlStore(backend)

# search
search = store.searcher()

search_defectinfo = search.variable(DefectInfo)
search_defectcell = search.variable(DefectCell)
search_screencell = search.variable(ScreenCell)

search.add(search_defectinfo.key == search_defectcell.key)
search.add(search_defectinfo.key == search_defectcell.key)
search.add(search_defectinfo.key == search_screencell.defect_key)

match setting:
    case "diamond":
        # filter out NV center in diamond
        search.add(search_defectinfo.key == -4871140043124584231)
        search.add(search_screencell.charge == -1)
        search.add(search_screencell.spin == 1.0)
    case "cao":
# filter out BiV center in CaO
        search.add(search_defectinfo.key == 6848080561178893677)
        search.add(search_screencell.charge == -1)
        search.add(search_screencell.spin == 1.0)
    case "tricky":
        # filter out tricky defect in SiC
        search.add(search_defectinfo.key == 7226619195624951017)
        search.add(search_screencell.charge == 0)
        search.add(search_screencell.spin == -1.0)

search.output(search_defectinfo, 'defectinfo')
search.output(search_defectcell, 'defectcell')
search.output(search_screencell, 'screencell')

for match in search:
    defect_info = match[0][0]
    defect_cell = match[0][1]
    screen_cell = match[0][2]

    print(defect_cell.defect_types)
    # sommerjobb
    print_atom_statistics(screen_cell.structure)
    molecule_structure, host_structure = screen_cell.structure.to_molecule(host_struct=host, layers=4, defect_type=defect_type, defect_cell=defect_cell, center='Replace', termination_group='H')
    molecule_structure.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
    host_structure.vis.show({'bonds': True, 'extbonds': False, 'polyhedra': False})
    print_atom_statistics(molecule_structure)



    print(molecule_structure.nested_layers['assignments'])
    print('Core layer')
    print(molecule_structure.get_assignments_layer('core_layer'))
    print('Host layer')
    print(molecule_structure.get_assignments_layer('host_layer'))
    print('Termination layer')
    print(molecule_structure.get_assignments_layer('termination_layer'))