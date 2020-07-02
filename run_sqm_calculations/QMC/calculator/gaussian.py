import os
from calculator.calculator import Calculator

link0_keys = ['mem',
              'nproc',
              'nprocshared']

route_self_keys = ['opt', 'freq']


class Gaussian(Calculator):
    """Gaussian calculator """

    implemented_properties = ['energy', 'freqencies']
    #program_path = "/opt/gaussian/g16"
    program_path = "/opt/gaussian/g09"
    
    default_parameters = {'method': 'pm3',
                          'basis': ''}


    def __init__(self, qmconf=None, label='g16', **kwargs):

        Calculator.__init__(self, qmconf, label, **kwargs)


    def write_input(self, write_file=False):
        """Write input file"""

        # write Gaussian header
        link0 = str()
        route = '# {} {}'.format(self.parameters['method'],
                                 self.parameters['basis'])
        
        for key, val in self.parameters.items():
            if key.lower() in link0_keys:
                link0 += '%{}={}\n'.format(key, val)

            elif key.lower() in route_self_keys:
                if (val.lower() == key.lower()):
                    route += (' ' + val)

                else:
                    if ',' in val:
                        route += ' {}({})'.format(key, val)
                    else:
                        route += ' {}={}'.format(key,val)

        # write molecular block
        mol_block = str(self.qmconf.charge) + "  " + str(self.qmconf.multiplicity) + "\n"

        structure = self.qmconf.structure
        atomic_symbols = self.qmconf.atomic_symbols

        for atom in zip(atomic_symbols, structure):
            symbol, pos = atom
            mol_block += "{}  {:10.5f} {:10.5f} {:10.5f}\n".format(symbol, *pos)

        input_string = link0 + route + '\n \n' + 'input prepared by QMMol \n\n' + mol_block + "\n"
        
        with open(self.calc_dir + '/' + self.prefix + ".com", 'w') as inp:
            inp.write(input_string)


    def calculate(self, quantities=['energy'], keep_files=False):
        """ Run calculation """
        Calculator.write_input(self, self.qmconf)
        
        # set enviroment.
        os.environ['GAUSS_SCRDIR'] = os.path.abspath(self.calc_dir)
        os.environ['GAUSS_EXEDIR'] = self.program_path # needed for g09.
        
        # Run calculation.
        self.write_input() # write input

        #command = self.program_path + "/g16 " + self.label + '.com'
        command = self.program_path + "/g09 < " + self.label + '.com'
        
        os.chdir(self.calc_dir) # move to working dir
        output = os.popen(command).read() # run calculation
        os.chdir('..') # get out of working dir
        
        # extract results from quantities.  
        results = Calculator.read_results(self, self.qmconf, output, quantities)

        # write output if keep_files = True
        if keep_files:
            with open(self.label + ".out", 'w') as f:
                f.write(output)

        self.clean(keep_files)

        return results


    def clean(self, keep_files):
        """Cleans from previus run"""
        
        if keep_files:
            os.rename(self.calc_dir+'/'+self.prefix + '.com', self.label + '.com')

        for f in os.listdir(self.calc_dir):
            f = self.calc_dir+'/'+f

            if os.path.isfile(f):
                os.remove(f)

        os.rmdir(self.calc_dir)
