#!/usr/bin/env python3
from pathlib import Path
import pysam
import defopt
import re
import attr
from typing import Dict
from typing import Optional
from collections import OrderedDict


@attr.s(frozen=True)
class PrimerMatch(object):
    '''Stores the primer match found in the f5/r5 tags produced by fgbio's IdentifyPrimer tool'''

    pair_id: str = attr.ib()
    primer_id: str = attr.ib()
    ref_name: str = attr.ib()
    start: int = attr.ib()
    end: int = attr.ib()
    positive_strand: bool = attr.ib()
    read_num: int = attr.ib()
    match_type: str = attr.ib()
    match_type_info: str = attr.ib()

    @staticmethod
    def get(tag_value: str) -> Optional['PrimerMatch']:
        if tag_value == "none":
            return None
        # ex. core_8.1216,core_8.1216_FOR,chr1:2327768-2327785,+,1,0,0,Location,0
        values = tag_value.split(',')
        ref_name, start, end = re.split('[:-]', values[2]) 
        return PrimerMatch(
            pair_id = values[0],
            primer_id = values[1],
            ref_name = ref_name,
            start = int(start),
            end = int(end),
            positive_strand = values[3] == '+',
            read_num = int(values[4]),
            match_type = values[5],
            match_type_info = values[6:]
        )


def _insert_length(f5: PrimerMatch, r5: PrimerMatch) -> int:
    '''Returns the insert length for a given primer pair.
    The primers must be mapped to the same chromosome. 
    '''
    assert f5.ref_name == r5.ref_name
    min_start = min(f5.start, r5.start)
    max_end = max(f5.end, r5.end)
    return max_end - min_start + 1


def _classify_primer_pair(f5: PrimerMatch, r5: PrimerMatch,
                          min_insert_length: int, max_insert_length: int) -> str:
    if f5 is None and r5 is None:  # no primer matches
        return "NoMatch"
    elif f5 is None or r5 is None:  # only one primer match
        return "Single"
    elif f5.pair_id != r5.pair_id:  # from different primer pairs!
        # not from the same pair; cross-dimer if the template/product size is too small, 
        # otherwise non-canonical
        if f5.ref_name != r5.ref_name:
            return "CrossDimer"
        else:
            insert_length = _insert_length(f5, r5) 
            if insert_length < min_insert_length or max_insert_length < insert_length:
                return "CrossDimer"
            else:
                return "NonCanonical"
    elif f5.primer_id == r5.primer_id:  # the same primer pair!
        return "SelfDimer"
    elif f5.positive_strand == r5.positive_strand:  # same primer pair, but the same strand (odd)
        return "NonCanonical"
    else:  # from the same primer pair, but different primers on opposite strands
        return "Canonical"


def _write_metrics(path: Path, label: str, counter: Dict[str, int]) -> None:
    total = 0
    for count in counter.values():
        total += count
    with path.open('w') as fh:
        fh.write(f'{label}\tcount\tfraction\n')
        for key, count in counter.items():
            fraction = count / float(total)
            fh.write(f'{key}\t{count}\t{fraction:.4f}\n')


def main(*, in_bam: Path, out_bam: Path, metrics: Path,
         min_insert_length: int = 50, max_insert_length: int = 250,
         primer_pair_tag: str = "pp",
         dimer_type_tag: str = "dt") -> None:
    '''Classifies how primers from paired end reads relate.

    The input should be the BAM produced fgbio's IdentifyPrimers tool.

    The records in the output BAM will be annotated with the classification given to the
    primer pair match.  For a given primer pair, the classifications are:
    * Canonical: two primers that are from the same "canonical" pair, with one match to the 
    forward and one to the reverse
    * SelfDimer: the same primer matches both ends of a pair.
    * CrossDimer: two primers that map to different chromosomes, or do not fall within the 
    expected product size.
    * NonCanonical: two primers that are neither CrossDimer nor from the same "canonical" pair.
    * Single: a single primer match to one of the two ends of a pair.
    * NoMatch: no primer matches

    Args:
        in_bam: the BAM produced fgbio's IdentifyPrimers tool
        out_bam: the output BAM
        metrics: the output metrics path prefix
    '''

    assert len(primer_pair_tag) == 2, "--primer-pair-tag must be of length 2"

    # Open the input and output files
    source = pysam.AlignmentFile(str(in_bam), "rb", check_sq=False)  # check_sq for unmapped BAMs
    writer = pysam.AlignmentFile(str(out_bam), "wb", template=source)
    
    # Set up the counters
    match_type_counter = OrderedDict()
    for key in ["Canonical", "SelfDimer", "CrossDimer", "NonCanonical", "Single", "NoMatch"]:
        match_type_counter[key] = 0
    dimer_type_counter = OrderedDict()
    for key in ["Dimer", "NotADimer", "MaybeDimer", "NoInformation"]:
        dimer_type_counter[key] = 0
    match_type_to_dimer_type = {
        "Canonical": "NotADimer",
        "SelfDimer": "Dimer",
        "CrossDimer": "Dimer",
        "NonCanonical": "MaybeDimer",
        "Single": "MaybeDimer",
        "NoMatch": "NoInformation"
    }

    # The main loop
    for record in source:
        assert record.has_tag('f5'), f'No f5 tag found on record {record}'
        
        # Extract the primer matches
        f5: Optional[PrimerMatch] = PrimerMatch.get(record.get_tag('f5'))
        r5: Optional[PrimerMatch] = \
            PrimerMatch.get(record.get_tag('r5')) if record.has_tag('r5') else None

        # Classify them
        classification = _classify_primer_pair(f5, r5, min_insert_length, max_insert_length)
        dimer_type = match_type_to_dimer_type[classification]

        # Add the tag
        record.set_tag(primer_pair_tag, classification)
        record.set_tag(dimer_type_tag, dimer_type)
        
        # Add to the metrics (do not double-count)
        if record.is_read1:
            assert classification in match_type_counter
            match_type_counter[classification] += 1
            dimer_type_counter[dimer_type] += 1
        
        # write the record
        writer.write(record)

    # write the metrics
    metrics_match_type = Path(str(metrics) + ".matches.txt")
    metrics_dimer_type = Path(str(metrics) + ".dimers.txt")
    _write_metrics(metrics_match_type, "match_type", match_type_counter)
    _write_metrics(metrics_dimer_type, "dimer_type", dimer_type_counter)

    # clean up pick up put away, clean up everyday
    writer.close()
    source.close()


if __name__ == '__main__':
    defopt.run(main)

