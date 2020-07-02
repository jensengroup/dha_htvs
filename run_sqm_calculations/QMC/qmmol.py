import sys
import multiprocessing as mp
import copy

import numpy as np # not needed.

from qmconf import QMConf
from conformers.create_conformers import random_conformers
from conformers.create_conformers import systematic_conformers
from clustering.butina_clustering import butina_clustering_m

class QMMol:
    """ Molecule object.

    The molecule object contains qmconf conformers.
    """

    def __init__(self):

        self.conformers = []

        self.calc = None # calculator object, for all conformers
        self.rdkit_mol = None

        self.charge = None
        self.multiplicity = None
        self.charged_fragments = None

        self.initial_conformer = None
        self.label = None


    def add_conformer(self, input_mol=None, fmt='gaussian', label=None,
                      charge=0, multiplicity=1, read_freq=False,
                      charged_fragments=True, set_initial=False):
        """Read X input and add conformer. Furthermore updates
        charge and multiplicity """

        qmconf = QMConf(input_mol, fmt, label,
                        charge, multiplicity, read_freq,
                        charged_fragments)

        self.conformers.append(qmconf)

        # set initial conformer if True.
        if set_initial:
            self.initial_conformer = qmconf
            self.label = qmconf.label

            # set charge and multiplicity for QMMol
            self.charge = charge
            self.multiplicity = multiplicity
            self.charged_fragments = charged_fragments

    def remove_conformers(self):
        self.conformers = []

    def set_initial_conformer(self, qmconf):
        """Set initial conformer. Use this conformer to check 
        after molecular manipulations."""
        self.initial_conformer = qmconf

    def set_calculator(self, calc=None):
        """Attach calculator object to all confs"""

        if calc is not None:
            self._calc = calc

            for conf in self.conformers:
                calc = copy.copy(self._calc)
                calc.add_conf(conf)
                conf.calc = calc

    def get_calculator(self):
        """Return attached calculator object."""
        return self._calc

    calc = property(get_calculator, set_calculator)

    def optimize(self, num_procs=1, keep_files=False, 
                 quantities=['energy', 'structure']):
        """Run optimizations with external QM programs in parallel.
        This function is only here to make it clear that an
        optimization is performed.

        All it does is appending the 'opt' keyword to parameters,
        and structure to quantities. In case one forgot. 
        """
        
        for conf in self.conformers:
            if (self.calc.name in ['gaussian', 'xtb', 'orca'] 
                    and 'opt' not in conf.calc.parameters):
                # If both true add opt to parameters.
                conf.calc.parameters['opt'] = 'opt'

        if 'structure' not in quantities:
            quantities.append('structure')
        
        self.calculate(num_procs, keep_files, quantities)
        

    def calculate(self, num_procs=1, keep_files=False, quantities=['energy']):
        """This function runs the calculation defined with
        parameters in parallel."""
        
        inp = list(zip(self.conformers,
                       [keep_files]*len(self.conformers),
                       [quantities]*len(self.conformers)))

        with mp.Pool(num_procs) as p:
           updated_confs =  p.map(worker, inp)
        
        self.conformers = updated_confs


    def get_rdkit_mol(self, charged_fragments=True):
        """Create rdkit mol, with qmconf conformers as RDKit conformers."""

        for idx, conf in enumerate(self.conformers):
            
            try:
                if idx == 0:
                    rdkit_mol = conf.get_rdkit_mol(charged_fragments)
                else:
                    rdkit_conf = conf.get_rdkit_mol(charged_fragments).GetConformer()
                    rdkit_mol.AddConformer(rdkit_conf, assignId=True)
            except:
                sys.stderr.write("can't create RDKit conf")
                continue

        self.rdkit_mol = rdkit_mol

        return rdkit_mol

    def rdkit_to_qmconfs(self, rdkit, label='conf', remove_confs=True):
        """Append RDKit mol conformers to self.conformers.
        
        * rdkit - rdkit mol object with conformers.
        * label - prefix name for conformers (label - prefix-X).
        * remove_confs - Remove all old conformers. 
                         If False append conformers.
        """

        if remove_confs:
            self.remove_conformers()

        for idx, new_conf in enumerate(rdkit.GetConformers()):
            conf_name = str(label) +'-{}'.format(idx)

            self.add_conformer(input_mol=new_conf, fmt='rdkit', label=conf_name,
                               charge=self.charge, multiplicity=self.multiplicity,
                               read_freq=False)

    def create_random_conformers(self, num_confs=100, threads=1):
        """Create random conformers using RDKit embed function

        * num_confs - the number of random conformers created.
        * threads - the number of threads to utilize.
        """
        if self.initial_conformer == None:
            raise RuntimeError('No initial conformer set.')

        initial_rdkit_conf = self.initial_conformer.get_rdkit_mol()
        new_confs = random_conformers(initial_rdkit_conf, num_confs, threads)

        # if label is None, use label from init_conf
        self.rdkit_to_qmconfs(new_confs, self.label, remove_confs=True)


    def create_systematic_conformers(self, theta=120., threads=1,
                                     ff_variant='uff', max_iters=200,
                                     check_stero=True):
        """Create systematic conformers by changing dihedrals of
        rotatable bonds by theta degrees.

        * theta - number of degrees to rotate dihedrals.

        Varibels used to clean new conformers:
        * ff_variant - force-field to use (UFF or MMFF94).
        * max_iters - max steps in optimization.
        * check_stero - check if stereo chemistry changed with respect to
                        initial_conformer.
        """
        if self.initial_conformer == None:
            raise RuntimeError('No initial conformer set.')

        if self.charged_fragments == None:
            raise RuntimeError('charged_fragments for initial not conformer set.')

        # systematic conformers takes rdkit mol.
        init_conf_rdkit = self.initial_conformer.get_rdkit_mol()
        new_confs = systematic_conformers(init_conf_rdkit, init_conf_rdkit,
                                          theta, self.charged_fragments,
                                          ff_variant, max_iters, threads,
                                          check_stero)

        # if label is None, use label from init_conf.
        self.rdkit_to_qmconfs(new_confs, self.label, remove_confs=True)


    def cluster_conformers(self, method='butina', difference_matrix='tfd', 
                           threshold=0.001):
        """Cluster conformers, and only return unique conformers."""

        rdkit_mol = self.get_rdkit_mol(self.charged_fragments)
        new_rdkit_mol = butina_clustering_m(rdkit_mol, difference_matrix,
                                            threshold)

        
        self.rdkit_to_qmconfs(new_rdkit_mol, self.label, remove_confs=True)
    

    def nlowest(self, n, energy='energy'):
        """Return list of conformers with n lowest energies. """
        
        energy_list = list()
        for conf in self.conformers:
            energy_list.append( conf.results['energy'])

        sorted_energy_list = sorted(energy_list)
        
        nlowest_confs = list()
        for val in sorted_energy_list[:n]:
            conf_idx = energy_list.index(val)
            nlowest_confs.append(self.conformers[conf_idx])
        
        return nlowest_confs


def worker(inp):
    """Helper function to run external programs in parallel"""
    conf, keep, quant = inp
    return conf.conf_calculate(keep_files=keep, quantities=quant)
