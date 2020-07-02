import itertools
import copy

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem import rdMolTransforms

from calculator.rdkit import ff_clean

def RotatableBonds(mol):
    """ Find rotatable bonds """

    # Find indexes of atoms in dihedral angles.
    # Remember alcohols, thiols, and primery amines.
    rot_bonds_smarts = [
            "[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]",  # general rot  bonds
            "[*]~[*]-[O,S]-[#1]",                         # alchol, and thiols
            "[*]~[*]-[NX3;H2]-[#1]",                      # primery amines
            ]

    raw_dihedral_idx = list()
    for bond_smart in rot_bonds_smarts:
        raw_dihedral_idx += mol.GetSubstructMatches(Chem.MolFromSmarts(bond_smart))

    rot_bonds_idx = list()  # indexes of bonds to rotate around
    dihedrals = list()   # indexes of dihedral
    for k, i, j, l in raw_dihedral_idx:
        if (i,j) not in rot_bonds_idx:

            rot_bonds_idx.append((i,j))
            dihedrals.append((k,i,j,l))

    return dihedrals

def systematic_conformers(mol, init_rdkit_mol, theta=120.,
                          charged_fragments=True, ff_variant='UFF',
                          max_iters=200, threads=1, check_stero=False):
    """Systematic rotation of dihedral angles theta degrees"""

    # find dihedrals of parent molecule and create all combinations
    # where the angles are rotated theta degrees.
    dihedral_idx = RotatableBonds(mol)

    conf = mol.GetConformer()

    dihedrals = list()
    for k, i, j, l in dihedral_idx:
        mol_dihedral = rdMolTransforms.GetDihedralDeg(conf, k,i,j,l)

        new_dihedrals = [mol_dihedral + x*theta for x in range(int(360./theta))]
        dihedrals.append(new_dihedrals)

    dihedralCombs = list(itertools.product(*dihedrals))
    
    # Create the conformations according to angle combinations
    # in dihedralCombs.
    for idx, dihedrals in enumerate(dihedralCombs):

        for (k,i,j,l), angle in zip(dihedral_idx, dihedrals):
            rdMolTransforms.SetDihedralDeg(conf, k,i,j,l, angle )

        mol.AddConformer(conf,assignId=True)


    # Clean conformers. Test if structure is corrupted or the stereo
    # chemistry changed.
    ff_clean(mol, init_rdkit_mol, ff_variant, max_iters, threads,
        charged_fragments, check_stero)

    return mol


def random_conformers(m, numConfs, threads):
    """ Create numConfs new random conformers from m """

    AllChem.EmbedMultipleConfs(m, numConfs, numThreads=threads)
    AllChem.UFFOptimizeMoleculeConfs(m, numThreads=threads, maxIters=1000)
    return m




if __name__ == "__main__":

    smi = "C[C](O)CC[C@@H](C)OO"

    m = Chem.MolFromSmiles(smi)
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m)

    m.SetProp("_Name", "test")

    w = Chem.SDWriter("test.sdf")
    for x in systematic(m, theta=180.):
        w.write(x)
        #print(Chem.MolToSmiles(x, isomericSmiles=True))

    #Chem.DetectBondStereochemistry(m)
    #Chem.AssignStereochemistry(m, flagPossibleStereoCenters=True, force=True)
    #Chem.AssignAtomChiralTagsFromStructure(m)


    #print(Chem.MolToSmiles(m, isomericSmiles=True))
