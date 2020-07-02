import os
from calculator.calculator import Calculator

calc_type = ['opt', 'ohess', 'hess', 'grad']

class xTB(Calculator):
    """ xTB calculator """

    implemented_properties = ['energy']

    program_path = '/opt/xtb/5.8.1'

    default_parameters = {'method': 'gfn2',
                          'opt': 'opt',
                          'cpus': 1}


    def __init__(self, qmconf=None, label='xtb', **kwargs):

        Calculator.__init__(self, qmconf, label, **kwargs)


    def write_input(self):
        """Write input file. For xtb input is xyz file """

        input_string = str(len(self.qmconf.atomic_numbers)) + '\n' # num atoms
        input_string += str(self.qmconf.label) + '\n' # title


        structure = self.qmconf.structure
        atomic_symbols = self.qmconf.atomic_symbols

        for atom in zip(atomic_symbols, structure):
            symbol, pos = atom

            input_string += "{}  {:10.5f} {:10.5f} {:10.5f}\n".format(symbol, *pos)

        with open(self.calc_dir + '/' + self.prefix + ".xyz", 'w') as inp:
            inp.write(input_string)


    def input_cmd(self):
        """Prepare command to run xTB, and set enviroment """

        # set xTB enviroment
        os.environ['XTBHOME'] = self.program_path

        if 'cpus' not in self.parameters.keys():
            cpus = 1
        else:
            cpus = self.parameters.pop('cpus')

        os.environ['OMP_NUM_THREADS'] = str(cpus)
        os.environ['MKL_NUM_THREADS'] = str(cpus)
        
        # Create xTB command. Set method (gfn, gfn2, gfn2d3, etc.)
        route = ' -{} '.format(self.parameters['method'])

        # Run uhf if not singlet.
        if self.qmconf.multiplicity != 1:
            route += '-uhf {} '.format(self.qmconf.multiplicity - 1)

        if self.qmconf.charge != 0:
            route += '-chrg {} '.format(self.qmconf.charge)
        
        # remove opt if ohess in self.parameters
        if 'ohess' in self.parameters.keys() and 'opt' in self.parameters.keys():
            del self.parameters['opt']

        # Set calculation type and solvent.
        for key, val in self.parameters.items():
            key, val = str(key), str(val)

            if key.lower() in calc_type:
                if (val.lower() == key.lower()):
                    route += '-{} '.format(val)
                else:
                    route += '-{} {}'.format(key,val)
        
        cmd = '{}/xtb {} {}'.format(self.program_path, self.label + '.xyz',  route)

        return cmd


    def calculate(self, quantities=['energy'], keep_files=False):
        """Run xTB calculation and return dict of results"""
        Calculator.write_input(self, self.qmconf)

        self.write_input() # write input file

        command = self.input_cmd() # command to run xTB
        
        os.chdir(self.calc_dir) # move into working dir.
        output = os.popen(command).read() # run calculation.
        os.chdir('..') # get out of working dir.
        
        # extract results from quantities.
        results = Calculator.read_results(self, self.qmconf, output, quantities)

        if keep_files:
            with open(self.label + ".out", 'w') as f:
                f.write(output)

        self.clean(keep_files)

        return results


    def clean(self, keep_files):
        """Cleans up from previus run"""

        if keep_files:
            os.rename(self.calc_dir+'/'+self.prefix + '.xyz', self.label + '.xyz')

        for f in os.listdir(self.calc_dir):
            f = self.calc_dir+'/'+f

            if os.path.isfile(f):
                os.remove(f)

        os.rmdir(self.calc_dir)


    #    for fil in os.listdir(self.directory):
    #        os.remove(self.directory + '/' + fil)
    #
    #    # TODO: keep files. Don't input or output file.
