#!/usr/bin/env python3
from pathlib import Path
import defopt

    
def main(*, in_primers: Path, out_primers: Path) -> None:
    '''Reformats mapped primers for use with fgbio's IdentifyPrimers.

    The input primers file must have the following header and columns in this order:
    - pair_id, primer_id, sequence, chr, start, and end

    The output primers file will have the following header and columns:
    - pair_id, primer_id, sequence, ref_name, start, end, and strand

    The input "chr" column will be used for the output "ref_name" column.  The input "primer_id"
    column will be used to infer the output "strand" column, and must be either "_FOR" ("+" in
    the output) or "_REV" ("-" in the output).

    The start will be incremented to transform the primer coordinates from 0-base exclusive to
    1-based inclusive

    Args:
        in_primers: the input of primers
        out_primers: the output reformatted primers for use with fgbio's IdentifyPrimers
    '''

    with in_primers.open('r') as fh_in, out_primers.open('w') as fh_out:
        in_iter = iter(line.rstrip('\r\n').split('\t') for line in fh_in)
        in_header = next(in_iter)
        out_header = ['pair_id', 'primer_id', 'sequence', 'ref_name', 'start', 'end', 'strand']
        fh_out.write('\t'.join(out_header) + '\n')
        for tokens in in_iter:
            data = dict(zip(in_header, tokens))
            # 0-based exclusive to 1-based inclusive
            data['start'] = str(int(data['start']) + 1)
            data['ref_name'] = data['chr']
            primer_id = data['primer_id']
            if primer_id.endswith('_FOR'):
                data['strand'] = '+'
            elif primer_id.endswith('_REV'):
                data['strand'] = '-'
            else:
                raise ValueError(f'FOR/REV not found in {primer_id}')
            fh_out.write('\t'.join([data[key] for key in out_header]) + '\n')


if __name__ == '__main__':
    defopt.run(main)

