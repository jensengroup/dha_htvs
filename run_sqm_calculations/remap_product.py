from rdkit import Chem
from rdkit.Chem import AllChem
import copy

def map_atoms(reactant, product):
    """Returns new atom ordering in product that matches the reactant """

    reac = copy.deepcopy(reactant)
    prod = copy.deepcopy(product)

    Chem.Kekulize(reac,clearAromaticFlags=True)
    Chem.Kekulize(prod,clearAromaticFlags=True)
     
    reac = change_mol(reac)
    prod = change_mol(prod)
    
    # Break Bond in reactant, in order to compare to product.
    smarts_bond = Chem.MolFromSmarts('[CX4;H0;R]-[CX4;H1;R]')
    atom_idx = list(reac.GetSubstructMatch(smarts_bond))
    
    if len(atom_idx) != 0:
        bond = reac.GetBondBetweenAtoms(atom_idx[0], atom_idx[1])
        broken_bond_reac = Chem.FragmentOnBonds(reac, [bond.GetIdx()], addDummies=False)
        
        # find new atom order for product
        prod_order = prod.GetSubstructMatch(broken_bond_reac)
    
    else:
        prod_order = prod.GetSubstructMatch(reac)
    
    return prod_order


def change_mol(mol):
    """ Change molecule to only compare conectivity """

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            bond.SetBondType(Chem.BondType.SINGLE)

    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() == 0:
            atom.SetFormalCharge(0)

    return mol


def reorder_product(reac, prod):
    """ change atom order of product, to match reactant  """
    new_product_order = map_atoms(reac, prod)
    reordered_product = Chem.RenumberAtoms(prod, new_product_order)

    return reordered_product
    

if __name__ == "__main__":
    
    reactant = Chem.MolFromSmiles('N#CC1(C#N)C(c2ccc(F)cc2)=CC2=CC=CC=CC21')
    AllChem.EmbedMolecule(reactant)
    
    product = Chem.MolFromSmiles('N#CC(C#N)=C(C=C1C=CC=CC=C1)c1ccc(F)cc1')
    AllChem.EmbedMolecule(product)
    
    product = reorder_product(reactant, product)

    Chem.SDWriter("reactant.sdf").write(reactant)
    Chem.SDWriter("reordered_product.sdf").write(product)


