#!/usr/bin/env python3
from pathlib import Path
import defopt

    
def main(*, in_primers: Path, out_primers: Path) -> None:
    '''Reformats unmapped primers for use with fgbio's IdentifyPrimers.

    The input primers file must have no header and the following columns:
    1. The name of the primer.  The name of the primer should be period ('.') delimited, with
    five values.  The first value should be the pair ID, and be unique across all primer pairs.
    The fifth value should be either "SP1" or "SP2" to which primer it belongs (forward/reverse).
    All other values will be ignored.
    2. The primer sequence (in sequencing order).
    3. Ignored (typically a sequence with the primer and target, with context).
    
    The output primers file will have the following header and columns:
    - pair_id, primer_id, sequence, ref_name, start, end, and strand

    Args:
        in_primers: the input of primers
        out_primers: the output reformatted primers for use with fgbio's IdentifyPrimers
    '''

    with in_primers.open('r') as fh_in, out_primers.open('w') as fh_out:
        in_iter = iter(line.rstrip('\r\n').split('\t') for line in fh_in)
        out_header = ['pair_id', 'primer_id', 'sequence', 'ref_name', 'start', 'end', 'strand']
        fh_out.write('\t'.join(out_header) + '\n')
        for tokens in in_iter:
            name_values = tokens[0].split('.')
            strand = '+' if name_values[-1] == 'SP1' else '-' 
            pair_id = name_values[0]
            primer_id = tokens[0] + ('_FOR' if strand == '+' else '_REV')
            sequence = tokens[1].upper()
            ref_name = ''
            start = '0'
            end = '0'
            data = [pair_id, primer_id, sequence, ref_name, start, end, strand]
            fh_out.write('\t'.join(data) + '\n')


if __name__ == '__main__':
    defopt.run(main)

