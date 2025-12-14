# from makeLipidNDX

def hasPhosGroup(u, resid):
    has_phos = []
    molmd = u.select_atoms(f'resid {resid}')
    Ps = molmd.select_atoms('element P')
    if not Ps:
        return []
    else:
        P = Ps[0]
    Os = []
    for bond in P.bonds:
        if bond.partner(P).element == "O":
            Os.append(bond.partner(P).index)
    if len(Os) == 4:
        has_phos.append(P.index)
        has_phos.extend(Os)
    return has_phos
