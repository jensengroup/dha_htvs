import xyz2mol.xyz2mol as x2m

def read_xyz_out(content, quantity='structure'):
    """Read xyz file
    - quantity = 'structure' - structure from rdkit conf
    - quantity = 'atomic_numbers'

    """

    if quantity == 'structure':
        return read_structure(content)

    if quantity == 'atomic_numbers':
        return read_atomic_numbers(content)


def read_structure(content):

    content = content.strip().split('\n')
    
    # first two lines not part of structure.
    del content[:2]

    atomic_positions = []
    for line in content:
        line = line.split()
        atomic_positions.append(list(map(float,line[1:])))

    return atomic_positions


def read_atomic_numbers(content):
    """Read xyz content."""

    content = content.strip().split('\n')

    # first two lines not part of structure.
    del content[:2]

    atomic_symbols = []

    for line in content:
        line = line.split()
        atomic_symbols.append(line[0])

    # atom symbols to atom numbers
    atomic_numbers = list()
    for atom in atomic_symbols:
        atomic_numbers.append(x2m.get_atom(atom))

    return atomic_numbers
