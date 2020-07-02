import copy

from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina

def butina_clustering_m(rdkit_mol, difference_matrix='tfd', threshold=0.001):
	""" Clustering conformers with RDKit's Butina algorithem """

	# calculate difference matrix
	if difference_matrix.lower() == 'tfd':
		diffmat = TorsionFingerprints.GetTFDMatrix(rdkit_mol)

	if difference_matrix.lower() == 'rms':
		diffmat = AllChem.GetConformerRMSMatrix(rdkit_mol, prealigned=False)

	# cluster conformers
	num_confs = rdkit_mol.GetNumConformers()
	clt = Butina.ClusterData(diffmat, num_confs, threshold, 
							 isDistData=True, reordering=True)

	# new conformers
	centroid_idx = [c[0] for c in clt] # centroid indexes.
	
	new_rdkit_mol = copy.deepcopy(rdkit_mol)
	new_rdkit_mol.RemoveAllConformers()

	for idx in centroid_idx:
		centroid_conf = rdkit_mol.GetConformer(idx)
		new_rdkit_mol.AddConformer(centroid_conf, assignId=True)

	del rdkit_mol # delete old mol, is this nessesary?

	return new_rdkit_mol
