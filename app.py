import streamlit as st
st.set_page_config(layout="wide")
import io
from contextlib import redirect_stdout


# Backend

dna_pairing = {'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

enzymes = ['BspEI', 'BglII', 'XhoI', 'PaeR7I', 'Eco53kI', 'SacI', 'HindIII', 'EcoRI', 'PstI', 'SalI', 'AccI', 'Acc65I', 'KpnI', 'SacII', 'PspOMI', 'XmaI', 'TspMI', 'ApaI', 'SmaI', 'BamHI']

enzymes = sorted(enzymes)

data = {'Acc65I': {'temp': 37, 'sequence': 'GGTACC', 'overhang': 'GTAC'},
        'AccI': {'temp': 37, 'sequence': 'GTMKAC', 'overhang': 'MK'},
        'ApaI': {'temp': 37, 'sequence': 'GGGCCC', 'overhang': 'GGCC'},
        'BamHI': {'temp': 37, 'sequence': 'GGATCC', 'overhang': 'GATC'},
        'BglII': {'temp': 37, 'sequence': 'AGATCT', 'overhang': 'GATC'},
        'BspEI': {'temp': -1, 'sequence': 'TCCGGA', 'overhang': 'CCGG'},
        'Eco53kI': {'temp': -1, 'sequence': 'GAGCTC', 'overhang': ''},
        'HindIII': {'temp': 37, 'sequence': 'AAGCTT', 'overhang': 'AGCT'},
        'KpnI': {'temp': 37, 'sequence': 'GGTACC', 'overhang': 'GTAC'},
        'PaeR7I': {'temp': -1, 'sequence': 'CTCGAG', 'overhang': 'TCGA'},
        'PspOMI': {'temp': -1, 'sequence': 'GGGCCC', 'overhang': 'GGCC'},
        'PstI': {'temp': 37, 'sequence': 'CTGCAG', 'overhang': 'TGCA'}, 
        'SacI': {'temp': 37, 'sequence': 'GAGCTC', 'overhang': 'AGCT'},
        'SacII': {'temp': 37, 'sequence': 'CCGCGG', 'overhang': 'GC'},
        'SalI': {'temp': 37, 'sequence': 'GTCGAC', 'overhang': 'TCGA'},
        'SmaI': {'temp': 25, 'sequence': 'CCCGGG', 'overhang': ''},
        'TspMI': None,
        'XhoI': {'temp': 37, 'sequence': 'CTCGAG', 'overhang': 'TCGA'},
        'XmaI': {'temp': 37, 'sequence': 'CCCGGG', 'overhang': 'CCGG'},
        'EcoRI': {'temp': 37, 'sequence': 'GAATTC', 'overhang': 'AATT'}}

pairs = {'Acc65I': {'BglII': 100, 'EcoRI': 100, 'PstI': 100, 'SacI': 100, 'SacII': 100, 'SalI': 100, 'SmaI': 100, 'XhoI': 100,
                    'ApaI': 75, 'BamHI': 75, 'HindIII': 75, 'KpnI': 75,
                    'XmaI': 50, 'AccI': 25},
        'AccI': {'ApaI': 50, 'BamHI': 25, 'BglII': 25, 'BspEI': 0,
                 'Eco53kI': 0, 'EcoRI': 25, 'HindIII': 25, 'KpnI': 25, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 25,  'SacI': 50, 'SacII': 50, 'SalI': 25,
                 'SmaI': 25, 'TspMI': 0, 'XhoI': 25, 'XmaI': 50},
        'ApaI': {'BamHI': 75, 'BglII': 75, 'BspEI': 0,
                 'Eco53kI': 0, 'EcoRI': 50, 'HindIII': 50, 'KpnI': 75, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 50,  'SacI': 75, 'SacII': 100, 'SalI': 25,
                 'SmaI': 75, 'TspMI': 0, 'XhoI': 50, 'XmaI': 50},
        'BamHI': {'BglII': 75, 'BspEI': 0,
                 'Eco53kI': 0, 'EcoRI': 75, 'HindIII': 100, 'KpnI': 75, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 50,  'SacI': 100, 'SacII': 75, 'SalI': 50,
                 'SmaI': 75, 'TspMI': 0, 'XhoI': 75, 'XmaI': 75},
        'BglII': {'BspEI': 0,
                 'Eco53kI': 0, 'EcoRI': 50, 'HindIII': 75, 'KpnI': 25, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 50,  'SacI': 25, 'SacII': 75, 'SalI': 100,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 100, 'XmaI': 75},
        'BspEI': {'Eco53kI': 0, 'EcoRI': 0, 'HindIII': 0, 'KpnI': 0, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 0,  'SacI': 0, 'SacII': 0, 'SalI': 0,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 0, 'XmaI': 0},
        'Eco53kI': {'EcoRI': 0, 'HindIII': 0, 'KpnI': 0, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 0,  'SacI': 0, 'SacII': 0, 'SalI': 0,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 0, 'XmaI': 0},
        'EcoRI': {'HindIII': 75, 'KpnI': 0, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 100,  'SacI': 75, 'SacII': 100, 'SalI': 50,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 100, 'XmaI': 50},
        'HindIII': {'KpnI': 50, 'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 50,  'SacI': 100, 'SacII': 75, 'SalI': 25,
                 'SmaI': 50, 'TspMI': 0, 'XhoI': 75, 'XmaI': 100},
        'KpnI': {'PaeR7I': 0,
                 'PspOMI': 0, 'PstI': 25,  'SacI': 100, 'SacII': 25, 'SalI': 25,
                 'SmaI': 100, 'TspMI': 0, 'XhoI': 25, 'XmaI': 50},
        'PaeR7I': {'PspOMI': 0, 'PstI': 0,  'SacI': 0, 'SacII': 0, 'SalI': 0,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 0, 'XmaI': 0},
        'PspOMI': {'PstI': 0,  'SacI': 0, 'SacII': 0, 'SalI': 0,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 0, 'XmaI': 0},
        'PstI': {'SacI': 25, 'SacII': 100, 'SalI': 50,
                 'SmaI': 25, 'TspMI': 0, 'XhoI': 100, 'XmaI': 50}, 
        'SacI': {'SacII': 75, 'SalI': 25,
                 'SmaI': 100, 'TspMI': 0, 'XhoI': 25, 'XmaI': 50},
        'SacII': {'SalI': 50,
                 'SmaI': 0, 'TspMI': 0, 'XhoI': 100, 'XmaI': 50},
        'SalI': {'SmaI': 0, 'TspMI': 0, 'XhoI': 100, 'XmaI': 25},

        'SmaI': {'TspMI': 0, 'XhoI': 0, 'XmaI': 50},
        'TspMI': {'XhoI': 0, 'XmaI': 0},
        'XhoI': {'XmaI': 75},
        'XmaI': {}
        }


def add_pair(e1, e2):

    if e1 == e2:
        return 0 # diagonal

    if not data[e1] or not data[e2]:
        return -10 # no info
    
    elif data[e1]['temp'] == -1 or data[e2]['temp'] == -1:
        return -10 # no info
    
    elif data[e1]['overhang'] == '' or data[e2]['overhang'] == '':
        # print('Blunt!')
        return -30 # blunt
    
    elif data[e1]['temp'] != data[e2]['temp']:
        return -20 # different temps
    
    else:
    
        s1 = data[e1]['sequence']
        s2 = data[e2]['sequence']

        i1 = MCS.find(s1)
        i2 = MCS.find(s2)

        if i1 == i2:
            return -40 # overlap
        elif i2 > i1:

            if i2 - i1 <= len(s1):
                return -40 # overlap
        elif i1 > i2:
            if i1 - i2 <= len(s2):
                return -40 # overlap
            
        max_score = 0
        if e1 in pairs[e2]:
            max_score = max(pairs[e2][e1], max_score)
        if e2 in pairs[e1]:
            max_score = max(pairs[e1][e2], max_score)

        return max_score


def complement(sequence: str):

    comp_sequence = []
    for symbol in sequence:

        if symbol.lower() in dna_pairing:
            comp_sequence.append(dna_pairing[symbol.lower()])
            if symbol.isupper():
                comp_sequence[-1] = comp_sequence[-1].upper()
        else:
            comp_sequence.append(symbol)

    return ''.join(comp_sequence)


def clean(sequence: str):

    clean_sequence = []
    for symbol in sequence:

        if symbol.lower() in dna_pairing:
            clean_sequence.append(symbol)

    return ''.join(clean_sequence)


def add_direction(sequence: str, reverse: bool = False):

    if reverse:
        return "3'-" + sequence + "-5'"
    return "5'-" + sequence + "-3'"


def to_triplets(sequence: str):

    return [sequence[i:i+3] for i in range(0,len(sequence),3)]

MCS = 'TCCGGACTCAGATCTCGAGCTCAAGCTTCGAATTCTGCAGTCGACGGTACCGCGGGCCCGGGATCC'
AQPR5 = """ATGAAGAAGGAGGTGTGCTCC
GTGGCCTTCCTCAAGGCCGTGTTCGCAGAGTTCTTGGCCACCCTCATCTTCGTCTTCTTTGGCCTGGGCTCGGCCCTCAA
GTGGCCGTCGGCGCTGCCTACCATCCTGCAGATCGCGCTGGCGTTTGGCCTGGCCATAGGCACGCTGGCCCAGGCCCTGG
GACCCGTGAGCGGCGGCCACATCAACCCCGCCATCACCCTGGCCCTCTTGGTGGGCAACCAGATCTCGCTGCTCCGGGCT
TTCTTCTACGTGGCGGCCCAGCTGGTGGGCGCCATTGCCGGGGCTGGCATCCTCTACGGTGTGGCACCGCTCAATGCCCG
GGGCAATCTGGCCGTCAACGCGCTCAACAACAACACAACGCAGGGCCAGGCCATGGTGGTGGAGCTGATTCTGACCTTCC
AGCTGGCACTCTGCATCTTCGCCTCCACTGACTCCCGCCGCACCAGCCCTGTGGGCTCCCCAGCCCTGTCCATTGGCCTG
TCTGTCACCCTGGGCCACCTTGTCGGAATCTACTTCACTGGCTGCTCCATGAACCCAGCCCGCTCTTTTGGCCCTGCGGT
GGTCATGAATCGGTTCAGCCCCGCTCACTGGGTTTTCTGGGTAGGGCCCATCGTGGGGGCGGTCCTGGCTGCCATCCTTT
ACTTCTACCTGCTCTTCCCCAACTCCCTGAGCCTGAGTGAGCGTGTGGCCATCATCAAAGGCACGTATGAGCCTGACGAG
GACTGGGAGGAGCAGCGGGAAGAGCGGAAGAAGACCATGGAGCTGACCACCCGCTGA"""
AQPR5 = AQPR5.replace(' ', '').replace('\n', '')


def apply_enzymes(enzyme_1, enzyme_2):

    print('\n' + '-' * 100)

    v = add_pair(enzyme_1, enzyme_2)
    
    s1 = data[enzyme_1]['sequence']
    s2 = data[enzyme_2]['sequence']

    i1 = MCS.find(s1)
    i2 = MCS.find(s2)

    c_MCS = complement(MCS)
    t_MCS = ''.join([str(i % 3) for i in range(1, len(MCS)+1)])

    upper_MCS_split = MCS[:min(i1, i2)+1] + ' ' * (abs(i1 - i2)) + MCS[max(i1, i2)+1:]
    lower_MCS_split = c_MCS[:min(i1, i2)+5] + ' ' * (abs(i1 - i2)) + c_MCS[max(i1, i2)+5:]
    upper_triplets_MCS_split = t_MCS[:min(i1, i2)+1] + ' ' * (abs(i1 - i2)) + t_MCS[max(i1, i2)+1:]
    lower_triplets_MCS_split = t_MCS[:min(i1, i2)+5] + ' ' * (abs(i1 - i2)) + t_MCS[max(i1, i2)+5:]

    first, second = sorted([enzyme_1, enzyme_2], key=lambda x: MCS.index(data[x]['sequence']))

    print('Original MCS domain:\n')
    print(add_direction('EGFP-' + '-'.join(to_triplets(MCS))))
    print(add_direction('EGFP-' + '-'.join(to_triplets(c_MCS)), reverse=True))

    print(f'\n{enzyme_1} and {enzyme_2} are applied:\n')

    masked_triplets_upper = add_direction('EGFP-' + '-'.join(to_triplets(upper_MCS_split)))
    masked_triplets_lower = add_direction('EGFP-' + '-'.join(to_triplets(lower_MCS_split)))

    indicator_upper_1 = masked_triplets_upper.find('-   ')
    if masked_triplets_upper[indicator_upper_1-2] == ' ':
        indicator_upper_1 -= 3
    elif masked_triplets_upper[indicator_upper_1-1] == ' ':
        indicator_upper_1 -= 2
    else:
        indicator_upper_1 -= 1

    indicator_upper_2 = masked_triplets_upper[::-1].find('-   ')
    indicator_upper_2 = len(masked_triplets_upper) - indicator_upper_2 - 1

    if masked_triplets_upper[indicator_upper_2+2] == ' ':
        indicator_upper_2 += 2
    elif masked_triplets_upper[indicator_upper_2+1] == ' ':
        indicator_upper_2 += 1
    else:
        indicator_upper_2 += 0

    upper_strand = add_direction('EGFP-' + '-'.join(to_triplets(MCS)))
    upper_strand = upper_strand[:indicator_upper_1 + 1] + ' | ' + upper_strand[indicator_upper_1 + 1:indicator_upper_2+1] + ' | ' + upper_strand[indicator_upper_2+1:]
    lower_strand = add_direction('EGFP-' + '-'.join(to_triplets(c_MCS)), reverse=True)
    lower_strand = lower_strand[:indicator_upper_1 + 6] + ' | ' + lower_strand[indicator_upper_1 + 6:indicator_upper_2+6] + ' | ' + lower_strand[indicator_upper_2+6:]

    indicator_upper_1 += 2
    indicator_upper_2 += 5

    pointer_upper = [' ' if i not in [indicator_upper_1, indicator_upper_2] else '|' for i in range(len(upper_strand))]
    pointer_middle = [' ' for i in range(len(upper_strand))]
    pointer_lower = [' ' if i not in [indicator_upper_1+5, indicator_upper_2+5] else '|' for i in range(len(upper_strand))]

    for i in list(range(indicator_upper_1, indicator_upper_1+6)) + list(range(indicator_upper_2, indicator_upper_2+6)):
        pointer_middle[i] = '-'

    pointer_upper = ''.join(pointer_upper)
    pointer_middle = ''.join(pointer_middle)
    pointer_lower = ''.join(pointer_lower)

    min_i, max_i = sorted([i1, i2])

    names = (' ' * indicator_upper_1) + first + (' ' * (indicator_upper_2 - indicator_upper_1 - len(first))) + second
    print(names)
    print(pointer_upper)
    print(upper_strand)
    print(pointer_middle)
    print(lower_strand)
    print(pointer_lower)

    print('\nThe sticky ends of MCS look like:\n')

    print(add_direction('EGFP-' + '-'.join(to_triplets(upper_MCS_split))).replace('- ', '  ').replace(' -', '  '))
    print(add_direction('EGFP-' + '-'.join(to_triplets(lower_MCS_split)), reverse=True).replace('- ', '  ').replace(' -', '  '))

    print('\n' + '-' * 100)

    return {'first': first,
            'second': second,
            'upper': upper_MCS_split,
            'upper_triplets': upper_triplets_MCS_split,
            'lower': lower_MCS_split,
            'lower_triplets': lower_triplets_MCS_split}


def design_primer_for_upper_strand(first_enzyme: str, n_matched: int = 21):

    s = data[first_enzyme]['sequence']
    matched = complement(AQPR5)[-n_matched:]
    clamp_upper_1 = '[x]'
    site = s[::-1]
    clamp_upper_2 = '[z]'

    return {'match_for_upper': matched,
            'x': clamp_upper_1,
            'z': clamp_upper_2,
            'site_for_upper': site}


def design_primer_for_lower_strand(second_enzyme: str, n_matched: int = 21):

    s = data[second_enzyme]['sequence']
    matched = AQPR5[:n_matched]
    clamp_upper_1 = '[y]'
    site = s
    clamp_upper_2 = '[z]'

    return {'match_for_lower': matched,
            'y': clamp_upper_1,
            'z': clamp_upper_2,
            'site_for_lower': site}


def print_gene(show_n: int = 24):

    aqpr = AQPR5[:]
    c_aqpr = complement(AQPR5)

    upper_left = aqpr[:show_n]
    upper_right = aqpr[-show_n:]
    lower_left = c_aqpr[:show_n]
    lower_right = c_aqpr[-show_n:]

    upper_left = '-'.join(to_triplets(upper_left))
    lower_left = '-'.join(to_triplets(lower_left))
    upper_right = '-'.join(to_triplets(upper_right[::-1]))[::-1]
    lower_right = '-'.join(to_triplets(lower_right[::-1]))[::-1]

    if upper_right[0] != '-':
        upper_right = '-' + upper_right
    if lower_right[0] != '-':
        lower_right = '-' + lower_right
    if upper_left[-1] != '-':
        upper_left += '-'
    if lower_left[-1] != '-':
        lower_left += '-'
    
    upper = upper_left + '[AQPR5]' + upper_right
    lower = lower_left + '[AQPR5]' + lower_right
    upper = add_direction(upper)
    lower = add_direction(lower, reverse=True)

    return upper, lower


def main(enzyme_1, enzyme_2, n_matched: int=21):

    n_matched = int(n_matched)

    assert 3 <= n_matched <= 100, 'number of primer matching bps must be between 3 and 100'

    if enzyme_1 == enzyme_2:

        raise ValueError('You cannot choose the same enzyme twice. Choose another pair')

    v = add_pair(enzyme_1, enzyme_2)

    if v == 0:

        raise ValueError('Enzymes are not compatible. Choose another pair')
    
    if v == -10:

        raise ValueError('Not enough information about the enzyme compatibility. Choose another pair')
    
    if v == -20:

        raise ValueError('Enzymes operate at different temperatures. Choose another pair')
    
    if v == -30:

        raise ValueError('At least one enzyme does not produce sticky ends (blunt). Choose another pair')
    
    if v == -40:

        raise ValueError('Sequences of enzymes overlap. Choose another pair')

    d = apply_enzymes(enzyme_1, enzyme_2)

    primers = design_primer_for_upper_strand(d['first'], n_matched=n_matched) | design_primer_for_lower_strand(d['second'], n_matched=n_matched)
    
    u, l = print_gene(n_matched+3)

    print(f'\nThe AQPR5 sequence ({n_matched+3} bps shown in every direction):\n')
    print(u)
    print(l)

    print('\n' + '-' * 100)

    print('\nThe primer for the upper strand mimics the lower strand in its complementary part,')
    print('\tbut it also has a clamp [x] of unknown length, enzyme-split sequence and a free clamp [z]')
    upper_match = '-'.join(to_triplets(primers['match_for_upper'][::-1]))[::-1]
    p_upper = upper_match + '-' + primers['x'] + '-' + primers['site_for_upper'] + '-' + primers['z']
    p_upper = add_direction(p_upper, reverse=True)
    print('\n\t', p_upper)

    print('\nThe primer for the lower strand mimics the upper strand in its complementary part,')
    print('\tbut it also has a clamp [y] of unknown length, enzyme-split sequence and a free clamp [z]')
    lower_match = '-'.join(to_triplets(primers['match_for_lower']))
    p_lower = primers['z'] + '-' + primers['site_for_lower'] + '-' + primers['y'] + '-' + lower_match
    p_lower = add_direction(p_lower)
    print('\n\t', p_lower)

    print('\n' + '-' * 100)

    print('\nThe hybridization scheme:')
    print('(arrows show the direction of polymerization)\n')

    i_l = l.find('3')
    i_l_p = p_lower.find('y')

    offset_l_p = 0
    offset_l = i_l_p - i_l
    offset_u = offset_l

    i_u = u.find('-3')
    i_p_u = p_upper.find('-[')
    offset_u_p = offset_u + i_u - i_p_u

    print((' ' * (offset_u_p - 6)) + ' <--- ' + p_upper)
    print((' ' * offset_u) + u)
    print((' ' * offset_l) + l)
    print((' ' * offset_l_p) + p_lower + ' ---> ')

    from_five_prime = p_lower.removesuffix("-3'") + u[len(p_lower) - offset_u - 3:-3]
    from_five_prime += '-[X]-' + complement(primers['site_for_upper']) + "-[z]-3'"

    print('\n' + '-' * 100)

    from_three_prime = []
    skip = True
    for symbol in from_five_prime:

        if symbol == '[':
            from_three_prime.append('[')
            skip = True
        elif symbol == ']':
            from_three_prime.append(']')
            skip = False
        else:
            if skip:

                if symbol.isupper():

                    from_three_prime.append(symbol.lower())
                elif symbol.islower():
                    from_three_prime.append(symbol.upper())
                else:
                    from_three_prime.append(symbol)
            else:
                from_three_prime.append(complement(symbol))

    from_three_prime = ''.join(from_three_prime)
    from_three_prime = '3' + from_three_prime[1:]
    from_three_prime = from_three_prime[:-2] + "5'"

    print('\nThe PCR product looks like:')
    print('(x/X, y/Y, z/Z are complementary sequences)\n')

    print(from_five_prime)
    print(from_three_prime)

    print('\n' + '-' * 100)

    print('\nNow we again apply enzymes (to the PCR product)\n')

    indicator_upper_1, indicator_upper_2 = 7, len(from_five_prime) - 13

    upper_strand = from_five_prime[:indicator_upper_1 + 1] + ' | ' + from_five_prime[indicator_upper_1 + 1:indicator_upper_2+1] + ' | ' + from_five_prime[indicator_upper_2+1:]
    lower_strand = from_three_prime[:indicator_upper_1 + 5] + ' | ' + from_three_prime[indicator_upper_1 + 5:indicator_upper_2+5] + ' | ' + from_three_prime[indicator_upper_2+5:]

    indicator_upper_1 += 2
    indicator_upper_2 += 5

    pointer_upper = [' ' if i not in [indicator_upper_1, indicator_upper_2] else '|' for i in range(len(upper_strand))]
    pointer_middle = [' ' for i in range(len(upper_strand))]
    pointer_lower = [' ' if i not in [indicator_upper_1+4, indicator_upper_2+4] else '|' for i in range(len(upper_strand))]

    for i in list(range(indicator_upper_1, indicator_upper_1+5)) + list(range(indicator_upper_2, indicator_upper_2+5)):
        pointer_middle[i] = '-'

    pointer_upper = ''.join(pointer_upper)
    pointer_middle = ''.join(pointer_middle)
    pointer_lower = ''.join(pointer_lower)

    names = (' ' * indicator_upper_1) + d['first'] + (' ' * (indicator_upper_2 - indicator_upper_1 - len(d['first']))) + d['second']
    print(names)
    print(pointer_upper)
    print(upper_strand)
    print(pointer_middle)
    print(lower_strand)
    print(pointer_lower)

    upper_strand = add_direction(upper_strand[11:-15])
    lower_strand = add_direction(lower_strand[15:-11], reverse=True)

    offset_lower = 4

    print('\n...which results in the following sticky-end slice:\n')

    print(upper_strand)
    print((offset_lower * ' ') + lower_strand)

    print('\n' + '-' * 100)

    print('\nThe product and the recipient vector can be joined together:\n')

    upper_left = '-'.join(to_triplets(d['upper'][:d['upper'].find(' ')]))
    lower_left = '-'.join(to_triplets(d['lower'][:d['lower'].find(' ')]))
    upper_right = d['upper'].split(' ')[-1]
    lower_right = d['lower'].split(' ')[-1]
    upper_right_triplets = d['upper_triplets'].split(' ')[-1]
    lower_right_triplets = d['lower_triplets'].split(' ')[-1]

    if upper_right_triplets[0] == '1':
        upper_right = '-'.join(to_triplets(upper_right))
    elif upper_right_triplets[0] == '2':
        upper_right = upper_right[:2]+'-'+'-'.join(to_triplets(upper_right[2:]))
    else:
        upper_right = upper_right[:1]+'-'+'-'.join(to_triplets(upper_right[1:]))

    if lower_right_triplets[0] == '1':
        lower_right = '-'.join(to_triplets(lower_right))
    elif lower_right_triplets[0] == '2':
        lower_right = lower_right[:2]+'-'+'-'.join(to_triplets(lower_right[2:]))
    else:
        lower_right = lower_right[:1]+'-'+'-'.join(to_triplets(lower_right[1:]))

    u = upper_left.ljust(1 + max(len(upper_left), len(lower_left)))
    l = lower_left.ljust(1 + max(len(upper_left), len(lower_left)))

    u += upper_strand[3:-3] + (offset_lower * ' ')
    l += (offset_lower * ' ') + lower_strand[3:-3]

    u += upper_right.rjust(1 + max(len(upper_right), len(lower_right)))
    l += lower_right.rjust(1 + max(len(upper_right), len(lower_right)))
    
    print(add_direction(u))
    print(add_direction(l, reverse=True))

    print('\nAfter the ligase is applied:')

    u = u.replace(' ', '')
    l = l.replace(' ', '')

    upper_left_triplets = '-'.join(to_triplets(u.split('-[')[0].replace('-', '')))
    lower_left_triplets = '-'.join(to_triplets(l.split('-[')[0].replace('-', '')))

    upper_mid = ']-'.join(('-[' + '-['.join(u.split('-[')[1:])).split(']-')[:-1]) + ']-'
    lower_mid = ']-'.join(('-['+ '-['.join(l.split('-[')[1:])).split(']-')[:-1]) + ']-'

    upper_right_triplets = '-'.join(to_triplets(u.split(']-')[-1].split('-')[0]))
    lower_right_triplets = '-'.join(to_triplets(l.split(']-')[-1].split('-')[0]))

    upper_rightmost = '-'.join(u.split(']-')[-1].split('-')[1:])
    lower_rightmost = '-'.join(l.split(']-')[-1].split('-')[1:])

    u = upper_left_triplets + upper_mid + upper_right_triplets + '-' + upper_rightmost
    l = lower_left_triplets + lower_mid + lower_right_triplets + '-' + lower_rightmost
    print(add_direction(u))
    print(add_direction(l, reverse=True))

    y_clamp = len(u.split('-[')[0].split('-')[-1])
    y_clamp = f'3N + {3-y_clamp}'

    x_clamp = len(u.split(']-')[-1].split('-')[0])
    x_clamp = f'3N + {3-x_clamp}'

    print('\n' + '-' * 100)

    print(f'\nBy observing the frame, we can conclude that:\n\t- Y ~ {y_clamp}\n\t- X ~ {x_clamp}\n')

    print('Result\n')
    before = from_five_prime.upper().split('-[Y]')[0] + '-[Y]'
    after = '[X]-' + from_five_prime.upper().split('[X]-')[1]

    p_upper = p_upper[::-1]
    p_upper = p_upper.replace("'5", "5'").replace("'3", "'3'")
    p_upper = p_upper.replace(']x[', f'[{x_clamp}]').replace(']X[', f'[{x_clamp}]').replace(']y[', f'[{y_clamp}]').replace(']Y[', f'[{y_clamp}]').replace(']z[', f'[Any 6-8 bps]').replace(']Z[', f'[Any 6-8 bps]')
    p_lower = p_lower.replace('[x]', f'[{x_clamp}]').replace('[X]', f'[{x_clamp}]').replace('[y]', f'[{y_clamp}]').replace('[Y]', f'[{y_clamp}]').replace('[z]', f'[Any 6-8 bps]').replace('[Z]', f'[Any 6-8 bps]')
    before = before.replace('[x]', f'[{x_clamp}]').replace('[X]', f'[{x_clamp}]').replace('[y]', f'[{y_clamp}]').replace('[Y]', f'[{y_clamp}]').replace('[z]', f'[Any 6-8 bps]').replace('[Z]', f'[Any 6-8 bps]')
    after = after.replace('[x]', f'[{x_clamp}]').replace('[X]', f'[{x_clamp}]').replace('[y]', f'[{y_clamp}]').replace('[Y]', f'[{y_clamp}]').replace('[z]', f'[Any 6-8 bps]').replace('[Z]', f'[Any 6-8 bps]')

    example_upper = p_upper.replace('[Any 6-8 bps]', 'CGCGCG').replace('-[3N + 0]', '').replace('[3N + 1]', 'A').replace('[3N + 2]', 'AA')
    example_lower = p_lower.replace('[Any 6-8 bps]', 'CGCGCG').replace('-[3N + 0]', '').replace('[3N + 1]', 'A').replace('[3N + 2]', 'AA')
    before = before.replace('[Any 6-8 bps]', 'CGCGCG').replace('-[3N + 0]', '').replace('[3N + 1]', 'A').replace('[3N + 2]', 'AA')
    after = after.replace('[Any 6-8 bps]', 'CGCGCG').replace('-[3N + 0]', '').replace('[3N + 1]', 'A').replace('[3N + 2]', 'AA')

    print(f'Upper strand primer:')
    print('\t', p_upper)
    print(f'\tfor example: {example_upper} ({sum(s.lower() in dna_pairing for s in example_upper)} bps)')
    print(f'Lower strand primer')
    print('\t', p_lower)
    print(f'\tfor example: {example_lower} ({sum(s.lower() in dna_pairing for s in example_lower)} bps)')

    before = ''.join(s for s in before if s.lower() in dna_pairing)
    after = ''.join(s for s in after if s.lower() in dna_pairing)

    product = before + '\n' + '\n'.join([AQPR5[i:i+50] for i in range(0, len(AQPR5), 50)]) + '\n' + after

    print('\n' + '-' * 100)
    
    print('\nFASTA of the PCR product:')
    print('> PCR')
    print(product)

    print('\n' + '-' * 100)


# Frontend

st.title("Bioanalytics A5 Cloning")


length = st.slider("Length of the complementary part of a primer:", min_value=3, max_value=100, value=10, step=1)

options1 = sorted(['BspEI', 'BglII', 'XhoI', 'PaeR7I', 'Eco53kI', 'SacI', 'HindIII', 'EcoRI', 'PstI', 'SalI', 'AccI', 'Acc65I', 'KpnI', 'SacII', 'PspOMI', 'XmaI', 'TspMI', 'ApaI', 'SmaI', 'BamHI'])  # <-- replace with your own
options2 = sorted(['BspEI', 'BglII', 'XhoI', 'PaeR7I', 'Eco53kI', 'SacI', 'HindIII', 'EcoRI', 'PstI', 'SalI', 'AccI', 'Acc65I', 'KpnI', 'SacII', 'PspOMI', 'XmaI', 'TspMI', 'ApaI', 'SmaI', 'BamHI'])  # <-- replace with your own

choice1 = st.selectbox("Enzyme 1:", options1)
choice2 = st.selectbox("Enzyme 2:", options2)

if st.button("Submit"):

    buffer = io.StringIO()
    with redirect_stdout(buffer):

        try:
            main(choice1, choice2, length)
        except BaseException as e:
            print('ERROR:', e)

    output_text = buffer.getvalue()

    if not output_text.strip():
        st.info("The function did not print anything.")
    else:
        st.subheader("Output")
        st.code(output_text, language="text")
