import numpy as np
import xyz2mol.xyz2mol as x2m



def read_xtb_out(content, quantity='energy'):
    """Reads gaussian output file

    - quantity = 'structure' - final structure form output.
    - quantity = 'atomic_numbers' - atmoic numbers
    - quantity = 'energy' - final energy from output.
    """

    if quantity == 'structure':
        return read_structure(content)

    elif quantity == 'atomic_numbers':
        return read_atomic_numbers(content)

    elif quantity == 'energy':
        return read_energy(content)

    elif quantity == 'converged':
        return read_converged(content)
    
    elif quantity == 'HOF':
        return read_hof(content)


def read_converged(content):
    """Check if program terminated normally"""
    return True


def read_energy(content):
    """Read total electronic energy """

    for line in content.split('\n'):

        if 'total E' in line:
            energy = float(line.strip().split()[-1])

    return energy


def read_hof(content):
    """Read total electronic energy """

    for line in content.split('\n'):
        if 'HOF' in line:
            hof = float(line.strip().split()[0])

    return hof

def read_structure(content):
    """Read structure from output file """
    
    temp_items = content.split('final structure:')[1:]

    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]

        del lines[:3] # first 3 lines are header lines

        atom_positions = []

        for line in lines:
            
            line = line.strip()
            # if line is empty - mol block ends.
            if set(line) == set(''):
                break
            
            tmp_line = line.split()
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions
            try:
                atom_position = list(map(float, tmp_line[1:]))
            except:
                raise ValueError('Expected a line with one string and three floats.')
            
            atom_positions.append(atom_position)

    return atom_positions


def read_atomic_numbers(content):
    """Read structure from output file """

    temp_items = content.split('final structure:')[1:]

    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]

        del lines[:3] # first 3 lines are header lines

        atom_symbols = []

        for line in lines:
            line = line.strip()

            # is line is empty - mol block ends.
            if set(line) == set(''):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions
            try:
                atom_symbol = str(tmp_line[0])
            except:
                raise ValueError('Expected a line with one string and three floats.')

            atom_symbols.append(atom_symbol)

    # atom symbols to atom numbers
    atomic_numbers = list()
    for atom in atom_symbols:
        atomic_numbers.append(x2m.get_atom(atom))

    return atomic_numbers


if __name__ == '__main__':
    import sys

    with open(sys.argv[1], 'r') as out:
        output = out.read()


    print(read_energy(output))

    #for x in read_atomic_numbers(output):
    #    #print(x)
    #    pass
