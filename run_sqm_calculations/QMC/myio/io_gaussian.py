

def read_gaussian_out(content, quantity='energy'):
    """Reads gaussian output file
    
    - quantity = 'structure' - final structure form output.
    - quantity = 'atomic_numbers' - atmoic numbers
    - quantity = 'energy' - final energy from output.
    - quantity = 'frequencies' - frequencies if freq is called
    - quantity = 'intensities' - freqency intensities
    - quantity = 'normal_coordinates' - freqency normal coordinates
    """
    
    if quantity == 'structure':
        return read_structures(content)
    
    elif quantity == 'atomic_numbers':
        return read_atomic_numbers(content)

    elif quantity == 'energy':
        return read_energy(content)
    
    elif quantity == 'frequencies':
        return read_frequencies(content)

    elif quantity == 'intensities':
        return read_intensities(content)
    
    elif quantity == 'normal_coordinates':
        return read_normal_coordinates(content)

    elif quantity == 'converged':
        return read_converged(content)


def read_energy(content):
    """Read electronic energy """
    
    for lines in content.split("\n"):
            
        if "E(" in lines:
            energy = float(lines.split()[4])
    
    return energy


def read_structures(content):
    """Read optimised structure from content. """

    # only interested in opt structure, hence -1.
    temp_items = content.split('Standard orientation')[1:] 
    
    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]
        
        #first 5 lines are headers
        del lines[:5]
        
        atom_positions = [] 

        for line in lines:
            line = line.strip()            
            
            #if only - in line it is the end
            if set(line).issubset(set(['-', ' '])):
                break
            
            tmp_line = line.split()
            if not len(tmp_line) == 6:
                raise RuntimeError('Length of line does not match structure!')
            
            # read atoms and positions:
            try:
                atom_position = list(map(float, tmp_line[3:]))
            except:
                raise ValueError('Expected a line with three integers and three floats.')
        
            atom_positions.append(atom_position)
        
    return atom_positions


def read_atomic_numbers(content):
    """Read optimised structure from content. """

    # only interested in opt structure, hence -1.
    temp_items = content.split('Standard orientation')[1:] 
    
    for item_i in temp_items:
        lines = [ line for line in item_i.split('\n') if len(line) > 0]
        
        #first 5 lines are headers
        del lines[:5]
        
        atom_nums = []

        for line in lines:
            line = line.strip()            
            
            #if only - in line it is the end
            if set(line).issubset(set(['-', ' '])):
                break
            
            tmp_line = line.split()
            if not len(tmp_line) == 6:
                raise RuntimeError('Length of line does not match structure!')
            
            atom_n = int(tmp_line[1])
        
            atom_nums.append(atom_n)
        
    return atom_nums



def read_frequencies(content):
    """Read frequencies and IR intensities"""
    
    frequencies = []

    freq_block = content.split('and normal coordinates:')[-1]
    freq_block = [ line for line in freq_block.split('\n') if len(line) > 0]
    
    for line in freq_block:
        line = line.strip() 

        #if only - in line it is the end
        if set(line).issubset(set('-')):
            break
        
        if 'Frequencies' in line:
            frequencies += list(map(float, line.split()[2:]))
       
    return frequencies


def read_intensities(content):
    """Read frequencies and IR intensities"""
    
    intensities = []

    freq_block = content.split('and normal coordinates:')[-1]
    freq_block = [ line for line in freq_block.split('\n') if len(line) > 0]
    
    for line in freq_block:
        line = line.strip() 

        #if only - in line it is the end
        if set(line).issubset(set('-')):
            break
        
        if 'IR Inten' in line:
            intensities += list(map(float, line.split()[3:]))
    
    return intensities

def read_normal_coordinates(content):
    """Read normal coordinates from frequency calculation."""
    
    normal_coordinates = []
    
    # two freq blocks, but only need last one.
    freq_block = content.split('and normal coordinates:')[-1] 
    temp_freq_blocks = freq_block.split('Frequencies --')[1:]

    for block in temp_freq_blocks:
        lines = [ line for line in block.split('\n') if len(line) > 0]
        
        #first 5 lines are headers
        del lines[:5]
        
        # 3 normal coordinates per line. 
        freq1_coords = []
        freq2_coords = []
        freq3_coords = []

        for line in lines:
            line = line.strip().split()
            
            # if not at least 3 normal coords (min len 9) then end.
            if len(line) < 5:
                break

            freq1_coords.append(line[2:5])
            freq2_coords.append(line[5:8])
            freq3_coords.append(line[8:11])
        
        normal_coordinates += [freq1_coords, freq2_coords, freq3_coords]
    
    return normal_coordinates


def read_converged(content):
    """Check if calculation terminated correctly"""
     
    converged = False
    for lines in content.split("\n"):
         
        if 'Normal termination of Gaussian' in lines:
            converged = True
    
    return converged



if __name__ == '__main__':
    import sys

    with open(sys.argv[1], 'r') as f:
        output = f.read()

    print(read_gaussian_out(output, quantity='converged'))


