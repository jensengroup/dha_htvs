

from rdkit import Chem
from rdkit.Chem import AllChem


def read_smiles_out(content, quantity='structure'):
    """Read smiles string, and embed molecule to get qmconf.

    NB. If there are stereocenters in the molecule, you will only 
    explore one random stereoisomer. If the stereocenters are defined
    explicitly in the smiles string, the conformer will conform 
    to the given chirality. 
    """

    # This doesn't work. Does it have to?

    mol = Chem.MolFromSmiles(content)
    mol = Chem.AddHs(mol)

    if quantity == 'structure':
        AllChem.EmbedMolecule(mol, enforceChirality=True)
        AllChem.UFFOptimizeMolecule(mol)
        
        return mol.GetConformer().GetPositions().tolist()
    
    if quantity == 'atomic_numbers':
        
        atomic_numbers = list()
        for atom in mol.GetAtoms():
            atomic_numbers.append(atom.GetAtomicNum())

        return atomic_numbers

        




