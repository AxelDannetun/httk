Molecule model
==============

This branch extends the structure class to support the conversion of model of defective or non-defective material structures, such as diamond, calcium oxide, etc, to smaller molecule models. This is mainly done by two functions within the structure class **to_molecule_no_defect()** and **to_molecule()**. The functions are very similair and the main difference is that the function **to_molecule_no_defect()** is for perfect structures, while **to_molecule** is meant for defect structures and requires the perfect material as an input parameter. Both functions are describe in more details below.

A lot of the functionality for the aforementioned function can be found in the file **moleculeutils.py**.

Currently not all structures and defects are supported.

The following structures are supported:
* Diamonds (C)
* Calcium oxide (CaO)
* Silicon carbide (SiC)

The following defects are supported:
* Vacancy (Vac)
* Interstitial (Int) **(Not guaranteed to work)**
* Substitution (Replace)

to_molecule_no_defect()
------
**Paramters**
* layers
* structure_type
* ref (Default = FracVector((1, 1, 1), 2))
* termination_group (Default ='H')

This is a member method of the structure class and it return a new structure which is a molecule model of itself. It works by creating layers of atoms around the atom in the structure that is closest to the coordinate of the input parameter ref. Layer one will contain the atom that is closest to the coordinate ref, the N:th layer will contains all atoms that are distance d or closer to the N - 1 layer. The distance d is determined by the input paramter structure_type (so Diamond, Calcium oxide, and Silicon carbide has different values of d). The input paramter layers is an integer and which layers should be a part of the resulting molecule model. If termination_group == None then a structure is created of these atoms and returned. If not another layer is added to the molecule model where the atom elements are set to the termination group, by default hydrogen.


to_molecule()
----
**Paramters**
* host_struct
* layers
* defect_cell (Default = None)
* defect_type (Default = None)
* center (Default = None)
* termination_group (Default ='H')

This function is similair to **to_molecule_no_defect()**, start by reading its description. The difference is that this function is meant for defected structures. The function creates layers for the structure itself, aswell as the perfect structure (input parameter host_structure). The layers are created around one of the structures defect (Vac, Int, or Replace) determined by the input paramter center (default to the first found defect if not specified). The layers are created in the same way as **to_molecule_no_defect** differing only in how the center of the layers are determined. The defect_cell parameter contains information about which defects the structure has and there coordinates. If the input paramter layer (same as for the function **to_molecule_no_defect()**) is the integer N, then the molecule model will be made up of the atoms in the first N - 1 layers in the defected structure and the N:th layer from the perfect structure (host_struct). If the termination_group is not None then the coordinates of the termination atoms will come from the N+1 layers in the perfect structure.

How to use
----
See files **defect_script.py** and **host_script.py** for example of how to use the functions.