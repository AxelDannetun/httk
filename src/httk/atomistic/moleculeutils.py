#!/usr/bin/env python

from httk.core.vectors.fracvector import FracVector

from httk.core.vectors.vectorutils import *

#0.10825317547305482
#0.1083

def empty(coordgroups):
    for group in coordgroups:
        if group:
            return False
    return True

#NV: 013
#BiV: 0.21/0.19
def append_layer(new_layer, coordgroups, ref_coord, radius, ord):
    to_be_removed = [[] for i in range(len(coordgroups))]
    for group in range(len(coordgroups)):
        for coord in coordgroups[group]:
            if dist(coord, ref_coord, ord=ord) < radius:
                new_layer[group].append(coord)
                to_be_removed[group].append(coord)

    for group in range(len(to_be_removed)):
        for coord in to_be_removed[group]:
            coordgroups[group].remove(coord)

def atomlayer(coordgroups, ref=FracVector((1, 1, 1), 2), radius=0.13, ord="L2"):

    number_of_elements = len(coordgroups)

    for i in range(number_of_elements):
        for coord in coordgroups[i]:
            closest = coord
            group = i
            break
        break


    for i in range(number_of_elements):
        for coord in coordgroups[i]:
            if dist(coord, ref) < dist(closest, ref):
                closest = coord
                group = i

    layer = [[] for element in range(number_of_elements)]
    layer[group].append(closest)
    layers = [layer]
    coordgroups[group].remove(closest)
    while not empty(coordgroups):
        new_layer = [[] for element in range(number_of_elements)]
        for group in layer:
            for coord in group:
                append_layer(new_layer, coordgroups, coord, radius=radius, ord=ord)
        layer = new_layer
        if empty(layer):
            break
        layers.append(new_layer)

    return layers    


def closest_coord(coordgroups, defect_coord = FracVector((1, 1, 1), 2)):
    best = coordgroups[0][0]

    for group in coordgroups:
        for coord in group:
            if dist(coord, defect_coord, "L2") < dist(best, defect_coord, "L2"):
                    best = coord
    return best

def closest_coord2(coords, ref):
    best = coords[0]
    for coord in coords:
        if dist(coord, ref, "L2") < dist(best, ref, "L2"):
            best = coord

    return best

def center_coordgroups(coordgroups, ref):
    cellcenter = FracVector((1, 1, 1), 2)
    diff = ref - cellcenter
    
    centered_coordgroups = [[coord - diff for coord in coords] for coords in coordgroups]
    return centered_coordgroups

def terminationlayer(outer_layer, termination_coords):
    termination_layer = []
    for group in outer_layer:
        for atom in group:
            for group2 in termination_coords:
                for pos in group2:
                    if dist(atom, pos, "L2") < 0.13:
                        v = pos - atom
                        v = v* (1/norm(v, "L2"))
                        termination_layer.append(atom + 0.08*v)

    return termination_layer

def determine_norm_and_radius(defect_type):
    match defect_type:
        case 'NV' | 'Diamond':
            ord = 'L2'
            radius = 0.13
            max_layers = 5
        case 'BiV' | 'SiC':
            ord = 'Cube' 
            radius = 0.19
            max_layers = 3
        case 'Tricky':
            ord = 'L2'
            radius = 0.135
            max_layers = 5
        case _:
            raise Exception('Defect ' + defect_type + 'is not supported')
    return ord, radius, max_layers


def parse_defect_cell(defect_cell, element_to_idx):
    parsed_defects = []
    
    for defect, coord in zip(defect_cell.defect_types, defect_cell.defect_positions):
        parts = defect.split('_')
        
        info = {
            "original_string": defect,
            "index": parts[0],
            "type": None,
            "host_element": None,
            "defect_element": None,
            'coord' : coord
        }

        if "Vac" in defect:
            info["type"] = "Vac"  
            info["host_element"] = element_to_idx[ parts[2].split('(')[0] ]
            
        elif "Int" in defect:
            info["type"] = "Int" # Remove the extra atom
            info["defect_element"] = element_to_idx[ parts[2].split('(')[0] ]
            #TODO: This case is more complex and has additional information that needs to be parsed.

        else:
            info["type"] = "Replace"
            info["defect_element"] = element_to_idx[ parts[1] ]
            info["host_element"] = element_to_idx[ parts[2].split('(')[0] ]
            
        parsed_defects.append(info)
    
    return parsed_defects

def pre_process(defect_coordgroups, defect_infos):
    for defect_info in defect_infos:
        match defect_info['type']:
            case 'Vac':
                defect_coordgroups[defect_info['host_element']].append(defect_info['coord'])
            case 'Replace':
                defect_info['new_coord'] = closest_coord2(defect_coordgroups[defect_info['defect_element']], defect_info['coord'])
                defect_coordgroups[defect_info['host_element']].append(defect_info['coord'])
                #TODO: The coordinate of defect is the coordinate for a perfect material. This following codes assumes that the defect element
                #      is unique, i.e. that no other element in the structure is the same as the defect element.
                #       Might want to find the new coordinate by allowing an epsilon error.
                #defect_info['new_coord'] = defect_coordgroups[defect_info['defect_element']][0]            
                defect_coordgroups[defect_info['defect_element']].remove(defect_info['new_coord'])
            case 'Int':
                defect_info['coord'] = closest_coord2(defect_coordgroups[defect_info['defect_element']], defect_info['coord'])
                defect_coordgroups[defect_info['defect_element']].remove(defect_info['coord'])
            case _:
                raise Exception('The defect type' + defect_info['type'] + 'is not supported')
   


def post_process(new_coordgroups, defect_infos):
    for defect_info in defect_infos:
        match defect_info['type']:
            case 'Vac':
                new_coordgroups[defect_info['host_element']].remove(defect_info['coord'])
            case 'Replace':
                new_coordgroups[defect_info['host_element']].remove(defect_info['coord'])
                new_coordgroups[defect_info['defect_element']].append(defect_info['new_coord'])
            case 'Int':
                new_coordgroups[defect_info['defect_element']].append(defect_info['coord'])
            case _:
                raise Exception('The defect type' + defect_info['type'] + 'is not supported')
            
 

def merge_layers(atom_layers, layers):
    number_of_elements = len(atom_layers[0])

    new_coordgroups = [[] for group in range(number_of_elements)]

    for i in range(min(layers, len(atom_layers))):
        for group in range(len(new_coordgroups)):
            new_coordgroups[group] = new_coordgroups[group] + atom_layers[i][group]

    return new_coordgroups

def append_coordgroups(coordgroups, new_coordgroups):
    for group in range(len(new_coordgroups)):
            for coord in new_coordgroups[group]:
                coordgroups[group].append(coord)


def anchor(defect_infos, defect_type = 'Vac'):
    #TODO: Sophisticate calc
    for defect_info in defect_infos:
        if defect_info['type'] == defect_type:
            return defect_info['coord']
        
    return FracVector((1, 1, 1), 2)




