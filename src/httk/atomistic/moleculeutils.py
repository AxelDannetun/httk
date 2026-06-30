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

def terminationlayer(outer_layer, termination_coords, ref = FracVector((1, 1, 1), 2)):
    termination_layer = []
    #for group in outer_layer:
    #    for atom in group:
    #        for group2 in termination_coords:
    #            for pos in group2:
    #                if dist(atom, pos, "L2") < 0.13:
    #                    v = pos - atom
    #                    v = v* (1/norm(v, "L2"))
    #                    termination_layer.append(atom + 0.08*v)
    for group in termination_coords:
        for atom in group:
            v = ref - atom
            v = v* (0.02/norm(v, "L2"))
            termination_layer.append(atom + v)

    return termination_layer

def determine_norm_and_radius(defect_type):
    match defect_type:
        case 'NV' | 'Diamond':
            ord = 'L2'
            radius = 0.13#3.0
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
            info["host_element"] = element_to_idx[parts[2].split('(')[0]]
            
        elif "Int" in defect:
            info["type"] = "Int" # Remove the extra atom
            info["defect_element"] = element_to_idx[parts[2].split('(')[0]]
            #TODO: This case is more complex and has additional information that needs to be parsed.

        else:
            info["type"] = "Replace"
            info["defect_element"] = element_to_idx[parts[1]]
            info["host_element"] = element_to_idx[parts[2].split('(')[0]]
            
        parsed_defects.append(info)
    
    return parsed_defects

def pre_process(defect_coordgroups, defect_infos):
    #print(defect_coordgroups)
    #defect_coordgroups = copy_coordgroups(defect_coordgroups)
    for defect_info in defect_infos:
        match defect_info['type']:
            case 'Vac':
                host_group = find_closest_group(defect_coordgroups, defect_info['host_element'], defect_info['coord'])
                defect_info['host_group'] = host_group
                defect_coordgroups[host_group].append(defect_info['coord'])
            case 'Replace':
                defect_group = find_closest_group(defect_coordgroups, defect_info['defect_element'], defect_info['coord'])
                host_group = find_closest_group(defect_coordgroups, defect_info['host_element'], defect_info['coord'])
                defect_info['defect_group'] = defect_group
                defect_info['host_group'] = host_group
                defect_info['new_coord'] = closest_coord2(defect_coordgroups[defect_group], defect_info['coord'])
                defect_coordgroups[host_group].append(defect_info['coord'])
                defect_coordgroups[defect_group].remove(defect_info['new_coord'])
            case 'Int':
                defect_group = find_closest_group(defect_coordgroups, defect_info['defect_element'], defect_info['coord'])
                defect_info['defect_group'] = defect_group
                defect_info['coord'] = closest_coord2(defect_coordgroups[defect_group], defect_info['coord'])
                defect_coordgroups[defect_group].remove(defect_info['coord'])
            case _:
                raise Exception('The defect type' + defect_info['type'] + 'is not supported')
   


def post_process(new_coordgroups, defect_infos):
    for defect_info in defect_infos:
        match defect_info['type']:
            case 'Vac':
                new_coordgroups[defect_info['host_group']].remove(defect_info['coord'])
            case 'Replace':
                new_coordgroups[defect_info['host_group']].remove(defect_info['coord'])
                new_coordgroups[defect_info['defect_group']].append(defect_info['new_coord'])
            case 'Int':
                new_coordgroups[defect_info['defect_group']].append(defect_info['coord'])
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


def count_atoms_by_element(structure):
    """
    Count the number of atoms for each element in the structure.
    
    Args:
        structure: A Structure object with uc_reduced_coordgroups and assignments
    
    Returns:
        dict: A dictionary with element symbols as keys and atom counts as values.
    """
    atom_counts = {}
    coordgroups = structure.uc_reduced_coordgroups
    symbols = structure.assignments.symbols
    
    for element_idx, symbol in enumerate(symbols):
        if element_idx < len(coordgroups):
            count = len(coordgroups[element_idx])
            atom_counts[symbol] = atom_counts.get(symbol, 0) + count
    
    return atom_counts


def get_atom_statistics(structure):
    """
    Get detailed atom statistics including counts and percentages.
    
    Args:
        structure: A Structure object with uc_reduced_coordgroups and assignments
    
    Returns:
        dict: Dictionary containing:
            - 'counts': dict with element symbol -> count
            - 'total_atoms': total number of atoms
            - 'percentages': dict with element symbol -> percentage (%)
    """
    atom_counts = count_atoms_by_element(structure)
    total_atoms = getattr(structure, 'uc_nbr_atoms', sum(atom_counts.values()))
    
    percentages = {}
    if total_atoms > 0:
        for element, count in atom_counts.items():
            percentages[element] = (count / total_atoms) * 100
    
    return {
        'counts': atom_counts,
        'total_atoms': total_atoms,
        'percentages': percentages
    }


def print_atom_statistics(structure):
    """
    Print formatted atom statistics for the structure.
    Displays total atom count, count per element, and percentage distribution.
    
    Args:
        structure: A Structure object with uc_reduced_coordgroups and assignments
    """
    stats = get_atom_statistics(structure)
    
    print(f"Total atoms: {stats['total_atoms']}")
    print("Element counts:")
    for element in sorted(stats['counts'].keys()):
        count = stats['counts'][element]
        percentage = stats['percentages'][element]
        print(f"  {element}: {count} atoms ({percentage:.2f}%)")


def build_element_index(symbols):
    """
    Map each element symbol to all indices where it occurs.
    """
    element_to_idx = {}
    for i, element in enumerate(symbols):
        element_to_idx.setdefault(element, []).append(i)
    return element_to_idx


def find_closest_group(coordgroups, indices, target_coord, ord='L2'):
    """
    Find the coordgroup index for a given element that is closest to the target coordinate.
    """
    best_idx = None
    best_dist = float('inf')

    for idx in indices:
        for coord in coordgroups[idx]:
            d = dist(coord, target_coord, ord)
            if d < best_dist:
                best_dist = d
                best_idx = idx

    if best_idx is None:
        raise Exception('Unable to find matching element coordgroup for target coordinate')

    return best_idx


def copy_coordgroups(coordgroups, symbols=None, merge=False):
    """
    Return a copy of coordgroups.
    If merge is True, merge coordgroups for the same element symbol into one group.
    symbols must be provided if merge is True.
    """
    if merge:
        if symbols is None:
            raise ValueError("symbols must be provided when merge=True")
        from collections import defaultdict
        merged = defaultdict(list)
        for i, group in enumerate(coordgroups):
            symbol = symbols[i]
            merged[symbol].extend(group)
        # Return in the order of unique symbols as they first appear
        result = []
        seen = set()
        for symbol in symbols:
            if symbol not in seen:
                result.append(merged[symbol][:])  # copy the list
                seen.add(symbol)
        return result
    else:
        return [[coord for coord in group] for group in coordgroups]


def validate_host_layers(host_atom_layers, layers, require_termination=False):
    """
    Validate that the host atom layers contain the requested layer and optionally the termination layer.
    """
    min_layers = layers + 2 if require_termination else layers + 1
    if len(host_atom_layers) <= min_layers - 1:
        raise Exception(
            "Host structure does not contain enough layers for layers=%d. "
            "Need at least %d host layers but got %d." % (layers, min_layers, len(host_atom_layers))
        )


def append_termination_layer(new_defect_coordgroups, new_host_coordgroups,
                             host_atom_layers, layers, termination_group,
                             defect_core_end_idx, defect_host_end_idx, host_core_end_idx,
                             new_defect_assignments, new_host_assignments, ref):
    termination_groups = host_atom_layers[layers + 1]
    outer_layer_groups = host_atom_layers[layers]
    termination_layer = terminationlayer(outer_layer_groups, termination_groups, ref)

    new_defect_coordgroups.append(termination_layer)
    new_defect_assignments = new_defect_assignments + [termination_group]

    new_host_coordgroups.append(termination_layer)
    new_host_assignments = new_host_assignments + [termination_group]

    defect_tags = {'host_layer': defect_core_end_idx, 'termination_layer': defect_host_end_idx}
    host_tags = {'termination_layer': host_core_end_idx}

    return new_defect_coordgroups, new_defect_assignments, new_host_coordgroups, new_host_assignments, defect_tags, host_tags




