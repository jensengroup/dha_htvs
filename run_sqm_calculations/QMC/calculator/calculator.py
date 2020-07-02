import os
from subprocess import Popen, PIPE

def get_calculator(name):
    """Return the calculator class."""

    classname = name.title()
    module = __import__('calculator.' + name, {}, None, [classname])
    Calculator = getattr(module, classname)

    return Calculator


implemented_programs = ["gaussian"]


class Calculator:
    """Base class for all calculators"""

    implemented_properties = []
    default_parameters = {}

    def __init__(self, qmconf=None, label=None, **kwargs):

        self.qmconf = None
        self.results = {}
        self.parameters = None

        if 'parameters' in kwargs:
            self.parameters = kwargs.pop('parameters')

        elif self.parameters is None:
            # Use default parameters if they were not read :
            self.parameters = self.default_parameters

        # link qmconf to calculator object, and calculator to qmconf.
        if qmconf is not None:
            qmconf.calc = self   # set calculator in qmconf.
            self.qmconf = qmconf # set qmconf in calculator.


        self.label = None
        self.directory = None
        self.prefix = None

        if not hasattr(self, 'name'):
            self.name = self.__class__.__name__.lower()


    def set_label(self, label):
        """ Set label and convert label to directory and prefix.

        Examples:

        * label = 'abc': (directory='.', prefix='abc')
        * label = 'dir1/abc: (directory='dir1', prefix='abc')

        """

        if self.qmconf is not None:
            self.label = self.qmconf.label
        else:
            self.label = label


        if self.label is None:
            self.directory = None
            self.prefix = None

        else:
             self.directory, self.prefix = os.path.split(self.label)

             if self.directory == '':
                  self.directory = os.curdir


    def add_conf(self, qmconf, use_label=True):
        """ Add qmconf to calculator object"""
        # TODO: is this in use??

        # update self.qmconf
        self.qmconf = qmconf

        # use qmconf labels
        if use_label:
            self.set_label(qmconf.label)

    
    def get_reader(self, fmt):
        """Return reader function from myio.X"""

        out_reader = 'read_' + fmt + '_out'
        reader = __import__('myio.io_' + fmt, {}, None, [out_reader])

        return getattr(reader, out_reader)


    def read_results(self, qmconf, content, quantities):
        """This functions reads results from output. The reader is
        found from the calculator class name. 
        
        returns results as a dictonary.

        Checks i 3 steps. 
            1) did the calculations converge, extract data.
            2) if no convergence, but the results can still be
               found i.e. running out of optimization steps.
               extract last structure, energy etc.
            3) Can't find results. Use old structure and set 
               properties to 99999.9.
        """
        
        properties = ['energy']
        reader = self.get_reader(self.name)
        
        results = {} # to collect results
        if reader(content, quantity='converged'):
            for quant in quantities:
                results[quant] = reader(content, quantity=quant)
            
            results['converged'] = True

        # not converged properly
        else:
            # opt calculation failes, get last step energy and structure.
            try:
                for quant in quantities:
                    results[quant] = reader(content, quantity=quant)
            
            # calculation failed for some other reason. 
            # set properties to 99999.9 and the structure to
            # initial structure.
            except UnboundLocalError:
                for quant in quantities:
                    if quant in properties:
                        results[quant] = 99999.9

                    elif quant == 'structure':
                        results['structure'] = qmconf.structure

                    else:
                        print("output for, {}, is not read properly".format(quant))
               
            results['converged'] = False

        return results


    def write_input(self, qmconf, properties=None):
        """ Write input file(s).

        Call this method first in subclasses so that directories are
        created automatically."""
        
        self.set_label(qmconf) # this is not the pretiest hack, but it works.

        if self.directory != os.curdir and not os.path.isdir(self.directory):
            os.makedirs(self.directory)

        # create directory to contain conf calculations.
        self.calc_dir = '{}/{}'.format(self.directory, self.prefix)
        if not os.path.isdir(self.calc_dir):
            os.makedirs(self.calc_dir)
