#import xyz2mol.xyz2mol as x2m

def read_orca_out(content, quantity='energy'):
    """Reads ORCA output

    - quantity = 'structure' - final structure form output.
    - quantity = 'atomic_numbers' - atmoic numbers.
    - quantity = 'energy' - final energy from output.
    - quantity = 'ts_scan_structures' - structures from surface scan.
    - quantity = 'ts_scan_energies' - energies from surface scan.
    - quantity = 'ts_guess' - guess ts structure.
    - 
    """

    if quantity == 'structure':
        return read_structure(content)

    elif quantity == 'atomic_numbers':
        return read_atomic_numbers(content)

    elif quantity == 'energy':
        return read_energy(content)

    elif quantity == 'ts_scan_structures':
        return read_ts_scan_structures(content)

    elif quantity == 'ts_scan_energies':
        return read_ts_scan_energies(content)

    elif quantity == 'ts_guess':
        return read_ts_guess(content)
    
    elif quantity == 'ts_guess_energy':
        return read_ts_guess_energy(content)

    elif quantity == 'converged':
        return read_converged(content)


def read_energy(content):
    """Read electronic energy """

    for line in content.split('\n'):

        if 'FINAL SINGLE POINT ENERGY' in line:
            try:
                electronic_energy = float(line.split()[4])
            except:
                raise ValueError('wrong format.')

    return electronic_energy


def read_structure(content):
    """Read optimized structure """

    temp_items = content.split('CARTESIAN COORDINATES (ANGSTROEM)')[1:]

    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]

        # first line is a header
        del lines[:1]

        atom_positions = []
        for line in lines:
            line = line.strip()

            if set(line).issubset(set('-')):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            try:
                atom_position = list(map(float, tmp_line[1:]))
            except:
                raise ValueError('Expected a line with one str and three floats.')

            atom_positions.append(atom_position)

    return atom_positions


def read_atomic_numbers(content):
    """Read atomic numbers """

    temp_items = content.split('CARTESIAN COORDINATES (ANGSTROEM)')[1:]

    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]

        # first line is a header
        del lines[:1]

        atomic_symbols = []
        for line in lines:
            line = line.strip()

            if set(line).issubset(set('-')):
                break

            tmp_line = line.split()
            if not len(tmp_line) == 4:
                raise RuntimeError('Length of line does not match structure!')

            try:
                atom_symbol = tmp_line[0]
            except:
                raise ValueError('Expected a line with one str and three floats.')

            atomic_symbols.append(atom_symbol)

    # atom symbols to atom numbers
    atomic_numbers = list()
    for atom in atomic_symbols:
        atomic_numbers.append(x2m.get_atom(atom))

    return atomic_numbers


def read_ts_scan_structures(content):
    """ Read surface scan structures """

    scan_steps = content.split('ORCA OPTIMIZATION COORDINATE SETUP')[1:]

    step_structures = list()
    for step in scan_steps:
        step_structures.append(read_structure(step))

    return step_structures


def read_ts_scan_energies(content):
    """ """
    temp_split = content.split('RELAXED SURFACE SCAN RESULTS')[1:][0]
    temp_split = temp_split.split('\n')

    del temp_split[:6]

    scan_energies = list()
    for line in temp_split:
        line = line.strip()
        if set(line).issubset(set()):
            break

        scan_energies.append(float(line.split()[1]))

    return scan_energies


def read_ts_guess(content):
    """Read TS guess """

    energies = read_ts_scan_energies(content)
    structures = read_ts_scan_structures(content)

    ts_guess_idx = energies.index(max(energies))

    return structures[ts_guess_idx]

def read_ts_guess_energy(content):
    """ Read ts guess energy """
    return max(read_ts_scan_energies(content))


def read_converged(content):
    """Check if calculation converged"""
    
    converged = False
    for line in content.split('\n'):

        if 'ORCA TERMINATED NORMALLY' in line:
            converged = True
    
    return converged



if __name__ == '__main__':
    import sys

    with open(sys.argv[1], 'r') as f:
        output = f.read()
    
    print(f"{sys.argv[1].split('_')[0]},{read_orca_out(output, quantity='ts_guess_energy')}")
