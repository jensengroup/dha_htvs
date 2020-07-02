

def read_rdkit_out(rdkit_conf, quantity='structure'):
    """Reads rdkit conformer object

    - quantity = 'structure' - structure from rdkit conf
    - quantity = 'atomic_numbers' 
    """

    if quantity == 'structure':
        return read_structure(rdkit_conf)
    
    if quantity == 'atomic_numbers':
        return read_atomic_numbers(rdkit_conf)


def read_structure(rdkit_conf):
    """Get structure from rdkit conformer"""
    return rdkit_conf.GetPositions().tolist()


def read_atomic_numbers(rdkit_conf):
    """Get atomic numbers from rdkit conformer"""

    rdkit_mol = rdkit_conf.GetOwningMol()
    
    atomic_numbers = list()
    for atom in rdkit_mol.GetAtoms():
        atomic_numbers.append(atom.GetAtomicNum())
    
    return atomic_numbers


