import numpy as np
import MDAnalysis as mda
from typing import Optional, Tuple
from collections import Counter

#  Dependencies
#      - split_mem
#      - makeLipidNDX.py


def hasPhosGroup(atomgroup) -> np.ndarray | None:
    """
    Returns the indices of the FIRST phosphate P and its oxygens when successful

    Returns:
        None -> atomgroup does not contain a phosphate
        list -> atomgroup contains a phosphate (list contains indices of the phosphate and its oxygens)
            [<P_index> <O_index> <O_index> <O_index> <O_index>]
    """
    Ps = atomgroup.select_atoms('element P')

    if not Ps:
        return None
    else:
        for P in Ps:
            Os = []
            for bond in P.bonds:
                if bond.partner(P).element == "O":
                    Os.append(bond.partner(P).index)
            if len(Os) == 4:
                has_phos = []
                has_phos.append(P.index)
                has_phos.extend(Os)
                return np.array(has_phos)
    return None


def getEsterCarbonyl(bridgeO, alcoholC):
    """
    bridgeO = mda.Atom
    alcoholC = mda.Atom
    Returns:
        list of length 2 with indices of the carbon and oxygen of the carbonyl attached to the alcohol
    """
    carbonyl_indices = [] # [<Cindex>, <Oindex>]
    other_atoms = 0
    for bond in bridgeO.bonds:
        if bond.partner(bridgeO) == alcoholC:
            continue
        elif bond.partner(bridgeO).element == "C":
            carbonyl_indices.append(bond.partner(bridgeO).index)
            for bond2 in bond.partner(bridgeO).bonds:
                if bond2.partner(bond.partner(bridgeO)) == bridgeO:
                    continue
                elif bond2.partner(bond.partner(bridgeO)).element == "O":
                    carbonyl_indices.append(bond2.partner(bond.partner(bridgeO)).index)
                elif other_atoms > 1:
                    return None
                else:
                    other_atoms += 1
        else:
            return None
    return carbonyl_indices


def isGlycerolCarbon(carbontype, atom):
    """
    atom = mda.Atom (the potential glycerol carbon)
    carbontype = 2:
        atom should be terminal glycerol carbon
    carbontype = 3:
        atom should be interstitial glycerol carbon
    """
    CandO = Counter()
    for bond in atom.bonds:
        CandO.update([bond.partner(atom).element])
    if not ("C" in CandO and "O" in CandO):
        return False
    if carbontype == 2:
        if CandO["C"] == 1 and CandO["O"] == 1:
            return True
        else:
            return False
    elif carbontype == 3:
        if CandO["C"] == 2 and CandO["O"] == 1:
            return True
        else:
            return False


def hasGlycerolAttached(atomgroup, phos_indices) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Check if atomindex in molecule is an O-C-C or O-C-2C
    Returns:
        None, None
        np.array, np.array -> indices of the glycerol, indices of the carbonyl carbons attached to it
    """


    glycerol_found = False
    carbonylC_indices = []
    current_atom = None # the index of the atom the decider is at
    phosmol = atomgroup.select_atoms('index ' + ' '.join([str(num) for num in phos_indices]))
    for oxygen in phosmol:
        if glycerol_found:
            break
        glyc_indices = []
        if oxygen.element != "O" or len(oxygen.bonds) == 1:
            continue
        if oxygen.bonds[0].partner(oxygen).element == "P":
            current_atom = oxygen.bonds[1].partner(oxygen)
        else:
            current_atom = oxygen.bonds[0].partner(oxygen)
        glyc_indices.append(current_atom.index) # add first glycerol C
        if isGlycerolCarbon(2, current_atom):
            for bond in current_atom.bonds:
                if bond.partner(current_atom).element == "C":
                    prev_atom = current_atom
                    current_atom = bond.partner(current_atom)
                    glyc_indices.append(current_atom.index) # add second glycerol C
                    break
            if isGlycerolCarbon(3, current_atom):
                for bond in current_atom.bonds:
                    if prev_atom == bond.partner(current_atom) or bond.partner(current_atom).element == "H":
                        continue
                    elif bond.partner(current_atom).element == "C":
                        next_atom = bond.partner(current_atom)
                    else:
                        glyc_indices.append(bond.partner(current_atom).index) # add O bound to C2
                        CandO = getEsterCarbonyl(bond.partner(current_atom), current_atom)
                        carbonylC_indices.append(CandO[0])
                        glyc_indices.extend(CandO) # add C=O bound to C2O if present
                glyc_indices.append(next_atom.index) # add third glycerol C
                prev_atom = current_atom
                current_atom = next_atom
                if isGlycerolCarbon(2,current_atom):
                    glycerol_found = True
                    for bond in current_atom.bonds:
                        if bond.partner(current_atom).element == "O":
                            glyc_indices.append(bond.partner(current_atom).index) # add O bound to C3
                            CandO = getEsterCarbonyl(bond.partner(current_atom), current_atom)
                            carbonylC_indices.append(CandO[0])
                            glyc_indices.extend(CandO) # add C=O bound to C3O if present
    if not glycerol_found:
        glyc_indices = None
        carbonylC_indices = None
    return glyc_indices, carbonylC_indices


def isPhospholipid(atomgroup) -> bool:
    """
    Returns:
        True  -> atomgroup is a phospholipid
        False -> atomgroup is not a phospholipid (probably)

    Might not work when:
        - atomgroup contains multiple phosphates

    """
    phos_indices = hasPhosGroup(atomgroup)
    if phos_indices is not None:
        also_has_glycerol_attached, ester_carbonyl_Cs = hasGlycerolAttached(atomgroup, phos_indices)
        if also_has_glycerol_attached and ester_carbonyl_Cs:
            return True
    return False


def leafletResolution(atomgroup, ref_atomgroup=None, z_midpoint=None) -> list:
    """
    Takes a single atomgroup and returns two atomgroups in a list.
    The first is the upper leaflet, the second is the lower leaflet.
    """
    # Create two groups for every residue in the atomgroup
    if ref_atomgroup is None:
        ref_atomgroup = atomgroup
    if z_midpoint is None:
        z_midpoint = ref_atomgroup.center_of_geometry()[2]
    leaflet_p_ndx = []
    leaflet_n_ndx = []

    for residue in atomgroup.residues:
        if ref_atomgroup.select_atoms('resid ' + str(residue.resid)).center_of_geometry()[2] > z_midpoint:
            leaflet_p_ndx.extend(atomgroup.select_atoms('resid ' + str(residue.resid)).indices)
        else:
            leaflet_n_ndx.extend(atomgroup.select_atoms('resid ' + str(residue.resid)).indices)
    leaflet_p = atomgroup.select_atoms('index ' + ' '.join([str(num) for num in leaflet_p_ndx]))
    leaflet_n = atomgroup.select_atoms('index ' + ' '.join([str(num) for num in leaflet_n_ndx]))
    assert len(leaflet_p) + len(leaflet_n) == len(atomgroup), "Leaflet resolution failed."
    return [leaflet_p, leaflet_n]
