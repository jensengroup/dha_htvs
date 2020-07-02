


def read_dftb_out(content, quantity='energy'):
    """Reads DFTB (GAMESS) output

    - quantity = 'structure' - final structure from outputself.
    - quantity = 'atomic_numbers' - atomic numbers.
    - quantity = 'energy' - final electronic energy.
    - quantity = 'converged' - return True if exited normally.
    """

    if quantity == 'structure':
        return read_structure(content)

    elif quantity == 'atomic_numbers':
        return read_atomic_numbers(content)

    elif quantity == 'energy':
        return read_energy(content)

    elif quantity == 'converged':
        return read_converged(content)


def read_energy(content):
    """Read electronic energy """

    for line in content.split("\n"):

        if "NSERCH:" in line:
            electronic_energy = float(line.split()[3])

    return electronic_energy


def read_structure(content):
    """Read structure from indput """

    temp_items = content.split('***** EQUILIBRIUM GEOMETRY LOCATED *****')[1:]

    for item_i in temp_items:
        lines = item_i.split('\n')

        # first 4 lines are headers.
        del lines[:4]

        atom_positions = list()

        for line in lines:
            line = line.strip()

            # if empty line - end of mol block
            if set(line).issubset(set()):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 5:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions:
            try:
                atom_position = list(map(float, tmp_line[2:]))
            except:
                raise ValueError('Expected a line with three integers and three floats.')

            atom_positions.append(atom_position)

    return atom_positions


def read_atomic_numbers(content):
    """Read optimised structure from content. """

    temp_items = content.split('***** EQUILIBRIUM GEOMETRY LOCATED *****')[1:]

    for item_i in temp_items:
        lines = item_i.split('\n')

        # first 4 lines are headers.
        del lines[:4]

        atom_numbers = list()

        for line in lines:
            line = line.strip()

            # if empty line - end of mol block
            if set(line).issubset(set()):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 5:
                raise RuntimeError('Length of line does not match structure!')

            # read atoms and positions:
            try:
                atom_number = int(float(tmp_line[1]))
            except:
                raise ValueError('Expected a line with three integers and three floats.')

            atom_numbers.append(atom_number)

    return atom_numbers


def read_converged(content):
    """Check if program terminated normally"""
    
    converged = False
    for line in content.split("\n"):
        
        if 'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
            converged = True

    return converged


if __name__ == '__main__':
    import sys

    with open(sys.argv[1], 'r') as f:
        output = f.read()


    print(read_dftb_out(output, 'converged'))
