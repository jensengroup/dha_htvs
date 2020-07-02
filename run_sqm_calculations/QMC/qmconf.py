import os
import copy

import xyz2mol.xyz2mol as x2m

implmented_input_formats = ['xyz',
                            'gaussian',
                            'xtb',
                            'rdkit_mol',
                            'smiles']

class QMConf:

    def __init__(self, input_mol=None, fmt='gaussian', label=None,
                 charge=0, multiplicity=1, read_freq=False,
                 charged_fragments=True):
        """Conformer object.

        The conformer object is an isolated conformer. Information
        about the atoms (atomic numbers and positions) is stored as
        lists.

        in order to calculate energies, frequencies etc. a calculator 
        object has to be attached to the qmconf object."""

        self.charge = charge
        self.multiplicity = multiplicity
        self.charged_fragments = charged_fragments

        self.structure = None
        self.atomic_symbols = None
        self.atomic_numbers = None
        self.results = {}
        
        # create QMConf from input file X.
        if input_mol is not None:
            
            reader = self.get_reader(fmt)
            
            # check if input mol is a file. If rdkit mol skip.
            if fmt != 'rdkit':
                if os.path.isfile(input_mol):
                    with open(input_mol) as f:
                        input_mol = f.read()

            # read quantities.
            quantities = ["structure", "atomic_numbers"]

            if read_freq:
                quantities += ['frequencies', 
                               'intensities',
                               'normal_coordinates']

            for quant in quantities:
                setattr(self, quant, reader(input_mol, quant))
                
            # set atomic_symbols from atomic_numbers
            self.number2symbols()

        self._calc = None
        self.rdkit_mol = None
        self.label = label


    def get_reader(self, fmt):
        """Return reader function from myio.X"""
        
        out_reader = 'read_' + fmt + '_out'
        reader = __import__('myio.io_' + fmt, {}, None, [out_reader])
        
        return getattr(reader, out_reader)
    

    def set_calculator(self, calc=None):
        """Attach calculator object."""
        calc_tmp = copy.copy(calc)
        
        calc_tmp.qmconf = self
        
        self._calc = calc_tmp
        
    def get_calculator(self):
        """Return attached calculator object."""
        return self._calc

    calc = property(get_calculator, set_calculator)


    def get_rdkit_mol(self, charged_fragments=True):
        """Return RDKit object. Uses xyz2mol."""

        quick = True # only option that works.

        rdkit_mol = x2m.xyz2mol(self.atomic_numbers, self.charge,
                                self.structure, charged_fragments, quick)

        self.rdkit_mol = rdkit_mol

        return rdkit_mol


    def number2symbols(self):
        """Convert atomic numbers to symbols, or sybols to numbers """

        atomic_list = x2m.__ATOM_LIST__

        # convert atomic numbers to symbols
        if self.atomic_symbols == None:

            atomic_symbols = list()

            for atom_num in self.atomic_numbers:
                atomic_symbols.append(atomic_list[atom_num - 1 ].title())

            self.atomic_symbols = atomic_symbols

        # convert atomic symbols to numbers
        if self.atomic_numbers == None:

            atomic_numbers = list()

            for atom_symbol in self.atomic_symbols:
                atom_symbol = atom_symbol.lower()
                atomic_numbers.append(atomic_list.index(atom_symbol) + 1)

            self.atomic_numbers = atomic_numbers
    

    def conf_calculate(self, keep_files=False, quantities=['energy']):
        """Run calculation on single conformation"""
    
        results = self.calc.calculate(keep_files=keep_files, 
                                      quantities=quantities)
        
        if 'structure' in quantities:
            self.structure = results.pop('structure')
        
        elif 'ts_guess' in quantities:
            self.structure = results.pop('ts_guess')
        
        self.results = results
        
        return self


    def write_xyz(self, to_file=True):
        """Write xyz string/file of qmconf"""

        xyz_string = str(len(self.atomic_numbers)) + '\n'
        xyz_string += str(self.label) + '\n'

        for symbol, pos in zip(self.atomic_symbols, self.structure):
            xyz_string += '{}  {:10.5f} {:10.5f} {:10.5f}\n'.format(symbol, *pos)

        if to_file:
            with open(self.label + '.xyz', 'w') as xyz:
                xyz.write(xyz_string)

        else:
            return xyz_string
