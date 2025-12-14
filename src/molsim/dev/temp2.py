# from utils/mda.py
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
