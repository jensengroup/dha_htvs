import numpy as np
from multiprocessing.pool import ThreadPool

from rdkit import Chem
from rdkit.Chem import AllChem


import xyz2mol.xyz2mol as x2m


def ff_clean(rdkit_mol, init_rdkit_mol, ff_variant='UFF', max_iters=200,
             threads=1, charged_fragments=True, check_stero=True):
    """ """
    if ff_variant.lower() == 'uff':
        AllChem.UFFOptimizeMoleculeConfs(rdkit_mol, numThreads=threads,
                                         maxIters=max_iters)


    if ff_variant.lower() in ['mmff94', 'mmff94s']:
        AllChem.MMFFOptimizeMoleculeConfs(rdkit_mol, numThreads=threads,
                                          maxIters=max_iters,
                                          mmffVariant=ff_variant.upper())


    # Are the structures corrupted or did the chirality change?
    # compare structures to init_rdkit_conf.

    # 1) check if structures are corrupted.
    with ThreadPool(threads) as pool:
        good_confs = pool.map(check_structures, rdkit_mol.GetConformers())
    
    # 2) check if stero chemistry changed.
    if check_stero:
        inp = zip(good_confs, [charged_fragments]*len(good_confs))

        with ThreadPool(threads) as pool:
            conf_stero = pool.map(check_stereo_chemistry, inp) 
        
        # compare to init_rdkit_mol.
        init_stero = Chem.FindMolChiralCenters(init_rdkit_mol)
        
        for idx, stero in enumerate(conf_stero):
            if stero != init_stero:
                good_confs[idx] = None


    # collect conformers into one rdkit mol
    for idx, conf in enumerate(good_confs):
        rdkit_mol.RemoveAllConformers() # Remove unchecked confs.

        good_confs = [x for x in good_confs if x is not None] # remove None

        for conf in good_confs:
            rdkit_mol.AddConformer(conf)

    return rdkit_mol
        

def check_structures(conf):
    """ Compute 3D distance matrix, and if min distance between 
    atoms <0.85 AA then is an ivallid  structure."""
    
    dist_matrix = AllChem.Get3DDistanceMatrix(conf.GetOwningMol(), conf.GetId())
    np.fill_diagonal(dist_matrix, np.inf) # fill diagonal
    
    if dist_matrix.min() > 0.85:
        return conf


def check_stereo_chemistry(inp):
    """ Check if stereo changed during ff opt"""
    
    conf, charged = inp

    if conf != None:
        new_mol = x2m.xyz2mol(*xyz_helper(conf), charged_fragments=charged,quick=True)
        conf_chiral = Chem.FindMolChiralCenters(new_mol)
        
        return conf_chiral


def xyz_helper(conf):
    """helper function: I need to reset in order to change check stero."""
    
    xyz_coordinates = conf.GetPositions().tolist()
    charge = Chem.GetFormalCharge(conf.GetOwningMol())
    
    atomicNumList = list() 
    for atom in conf.GetOwningMol().GetAtoms():
        atomicNumList.append(atom.GetAtomicNum())

    return atomicNumList, charge, xyz_coordinates

