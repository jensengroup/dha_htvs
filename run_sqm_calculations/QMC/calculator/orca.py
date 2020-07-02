import os
from calculator.calculator import Calculator


route_keys = ['opt', 'freq']

link_keys = ['geom scan']
link1_keys = ['']


class ORCA(Calculator):
    """ ORCA calculator """

    implemented_properties = ['energy', 'ts_guess']
    program_path = '/opt/orca/4.0.0/orca'

    default_parameters = {'method': 'pm3',
                          'basis': '',
                          'cpus': 1,
                          'mem': '4GB'
                          }

    def __ini__(self, qmconf=None, label='orca', **kwargs):

        Calculator.__init__(self, qmconf, label, **kwargs)


    def write_input(self, write_file=True):
        """Write input file """

        # Write ORCA header
        link = str()
        route = '! {} {}'.format(self.parameters.pop('method'),
                                 self.parameters.pop('basis'))

        # cpus
        if 'cpus' in self.parameters:
            cpus = self.parameters.pop('cpus')
        else:
            cpus = 1

        link += '%pal nprocs {} end \n'.format(cpus)

        # mem pr. core
        if 'mem' in self.parameters:
            total_mem = self.parameters.pop('mem')
            total_mem = total_mem[:-2] # remove GB

            mem_pr_core = int(round(float(total_mem + '000') / cpus, -3))

            link += '%maxcore {} \n'.format(mem_pr_core)


        for key, val in self.parameters.items():

            if key.lower() in link_keys:
                link += '%{} {} end end \n'.format(key,val)

            if key.lower() in route_keys:
                    route += (' ' + val)

        # write mol block
        mol_block = '*xyz {} {} \n'.format(self.qmconf.charge,
                                           self.qmconf.multiplicity)

        structure = self.qmconf.structure
        atomic_symbols = self.qmconf.atomic_symbols

        for atom in zip(atomic_symbols, structure):
            symbol, pos = atom
            mol_block += "{}  {:10.5f} {:10.5f} {:10.5f}\n".format(symbol, *pos)
        mol_block += '*'

        input_string = route + '\n' + link + '\n' + mol_block

        if write_file:
            with open(self.calc_dir + '/' + self.prefix + ".inp", 'w') as inp:
                inp.write(input_string)
        else:
            return input_string


    def calculate(self, quantities=['energy'], keep_files=False):
        """ Run calculation """
        Calculator.write_input(self, self.qmconf)

        self.write_input() # write input

        command = self.program_path +" "+ self.label + '.inp'
        
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
        """clean from previus run """

        if keep_files:
            os.rename(self.calc_dir+'/'+self.prefix + '.inp', self.label + '.inp')

        for f in os.listdir(self.calc_dir):
            f = self.calc_dir+'/'+f

            if os.path.isfile(f):
                os.remove(f)

        os.rmdir(self.calc_dir)
