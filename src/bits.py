char_to_bits = {
    'A': b'00000',
    'R': b'00001',
    'N': b'00010',
    'D': b'00011',
    'C': b'00100',
    'E': b'00101',
    'Q': b'00110',
    'G': b'00111',
    'H': b'01000',
    'I': b'01001',
    'L': b'01010',
    'K': b'01011',
    'M': b'01100',
    'F': b'01101',
    'P': b'01110',
    'S': b'01111',
    'T': b'10000',
    'U': b'10001',
    'W': b'10010',
    'Y': b'10011',
    'V': b'10100',
    'X': b'10101', 
    'B': b'10110',
    'Z': b'10111'
}
bits_to_char = {v: k for k, v in char_to_bits.items()}

def str_to_bits(string: str) -> bin:
    '''
    Convert a string to a bin

    Input:
        string:     (str)   string to convert
    Output:
        bin
    '''
    bit_arr = b'0'
    for c in string:
        if c in char_to_bits:
            bit_arr += char_to_bits[c]
        else:
            raise f'Error: character {c} not known'

    return bit_arr[1:]

def bits_to_str(bits: bin) -> str:
    '''
    Convert a bin to a string

    Inputs:
        bits:   (bin)  the bits to convert
    Ouputs:
        str
    '''
    # if its not divisable by 5, raise error
    if len(bits) % 5 != 0:
        raise f'Error: bit string should be in lengths divisable by 5'
    
    string = ''
    for i in range(0, len(bits), 5):
        string += bits_to_char[bits[i:i+5]]
        
    return string