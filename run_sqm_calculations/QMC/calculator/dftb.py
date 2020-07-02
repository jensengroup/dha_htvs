import os
from itertools import product, repeat

from calculator.calculator import Calculator

class dftb(Calculator):
    """ dftb calculator as implemented in GAMESS

    This can only perform geometry optimizations.
    
    NB. It is worth considering changin the SCRDIR in rungms, to
        a scratch folder on the node itself.
    """

    implemented_properties = ['energy']

    program_path = '/home/koerstz/opt/gamess/gamess-sep18R3/rungms'
    #program_path = '/home/charnley/opt/gamess/github-2017-08-29/rungms'
    skf_path = '/home/andersx/dftb/3ob-3-1'
    scr_dir = '/home/koerstz/scr'

    default_parameters = {'method': 'rohf',
                          'opt': 'opt'}


    def __init__(self, qmconf=None, label='dftb', **kwargs):
        Calculator.__init__(self, qmconf, label, **kwargs)

    def write_input(self):
        """Write input file for dftb (GAMESS) optimization. """
        qmconf = self.qmconf # makes life easier

        charge = qmconf.charge
        multiplicity = qmconf.multiplicity
        method = self.parameters['method']


        header = """ $system
    mwords=20
    modio=32
 $end

 $scf
    npunch=0
 $end

 $contrl
     runtyp=optimize
     scftyp={}
     icharg={}
     mult={}
     nprint=-5
 $end

 $statpt
     nstep=500
     opttol=5.0e-4
 $end

 $basis
     gbasis=dftb
 $end

 $dftb
     ndftb=3
     dampxh=.t.
     dampex=4.0
 $end

 $dftbsk \n""".format(method, charge, multiplicity)

        # write SKF block
        for atom_pair in product(set(qmconf.atomic_symbols),repeat=2):
            header += '    {0} {1}  {2}/{0}-{1}.skf" \n'.format(*atom_pair,
                                                             self.skf_path)

        header += ' $end \n \n'

        # Write Mol block
        atom_data = zip(qmconf.atomic_symbols,
                        qmconf.atomic_numbers,
                        qmconf.structure)

        mol_block = ' $DATA \n \n C1 \n'
        for symbol, number, pos in atom_data:
            atom_line = symbol + '\t' + str(number) + '\t'
            atom_line += '{:10.5f} {:10.5f} {:10.5f}\n'.format(*pos)

            mol_block += atom_line

        mol_block += ' $END'

        with open(self.calc_dir + '/' + self.prefix + ".inp", 'w') as inp:
            inp.write(header + mol_block)


    def calculate(self, quantities=['energy'], keep_files=False):
        """ Run calculation """
        Calculator.write_input(self, self.qmconf)

        self.write_input() # write input

        command = self.program_path +" "+ self.label + '.inp' + ' 2>' + self.scr_dir + '/' + self.label +'.log'

        os.chdir(self.calc_dir) # create working dir
        output = os.popen(command).read() # run calculation
        os.chdir('..') # get out of working dir

        # extract results from quantities.
        results = Calculator.read_results(self, self.qmconf, output, quantities)

        if keep_files:
            with open(self.label + ".out", 'w') as f:
                f.write(output)

        self.clean(keep_files)

        return results


    def clean(self, keep_files):
        """ Cleans from previus run """

        if keep_files:
            os.rename(self.calc_dir+'/'+self.prefix + '.inp', self.label + '.inp')
        
        # remove working dir.
        for f in os.listdir(self.calc_dir):
            f = self.calc_dir+'/'+f

            if os.path.isfile(f):
                os.remove(f)

        os.rmdir(self.calc_dir)

        # Also remove .dat file from scratch dir
        os.remove(self.scr_dir + '/' + self.prefix + '.dat')
        os.remove(self.scr_dir + '/' + self.prefix + '.log')
