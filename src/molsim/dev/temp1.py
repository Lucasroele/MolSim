# from makeLipidNDX
# Check if atomindex in molmd is a O-C-C or O-C-2C
def isGlycerolCarbon(carbontype, atom):
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
