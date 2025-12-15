# from utils/mda.py
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
def hasPhosGroup(atomgroup) -> np.ndarray | None:
    """
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
