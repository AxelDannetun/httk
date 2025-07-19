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
        layers.append(new_layer)

    return layers    


def closest_coord(coordgroups, defect_coord = FracVector((1, 1, 1), 2)):
    best = coordgroups[0][0]

    for group in coordgroups:
        for coord in group:
            if dist(coord, defect_coord, "L2") < dist(best, defect_coord, "L2"):
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
            max_layers = 2
        #case 'Tricky':
        #    ord = 'L2'
        #    radius = 0.15
        case _:
            raise Exception('Defect ' + defect_type + 'is not supported')
    return ord, radius, max_layers

def pre_process(defect_coordgroups, defect_type, defect_coords):
    match defect_type:
        case 'NV' : 
            defect_coordgroups[0].append(defect_coords[0])
            defect_coordgroups[0].append(defect_coords[1])
        case 'BiV':
            defect_coordgroups[0].append(defect_coords[0])
            defect_coordgroups[1].append(defect_coords[1])
        #case 'Tricky':
        case _:
            raise Exception('Defect ' + defect_type + 'is not supported')


def post_process(new_coordgroups, defect_type, defect_coords):
    match defect_type:
        case 'NV': 
            if (defect_coords[0] in new_coordgroups[0]):
                new_coordgroups[0].remove(defect_coords[0])
            if (defect_coords[1] in new_coordgroups[0]):
                new_coordgroups[0].remove(defect_coords[1])
        case 'BiV':
            if (defect_coords[0] in new_coordgroups[0]):
                new_coordgroups[0].remove(defect_coords[0])
            if (defect_coords[1] in new_coordgroups[1]):
                new_coordgroups[1].remove(defect_coords[1])
        #case 'Tricky':
        case _:
            raise Exception('Defect ' + defect_type + 'is not supported')

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




