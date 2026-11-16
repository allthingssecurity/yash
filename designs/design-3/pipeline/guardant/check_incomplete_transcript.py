"""Resolve an unavailable cDNA string with a genomic CDS/protein-prefix check."""
from pathlib import Path
import hashlib
import itertools
import json
import sys

import pysam
from pyensembl import Genome

root = Path('/data/results/TN25-112103/matched-normal-20261107')
config = json.loads((root / 'analysis-genome-config.json').read_text())
genome = Genome(**config)
plan = json.loads((root / 'normal-coding-plan.json').read_text())
fasta = pysam.FastaFile(plan['reference'])
transcript = genome.transcript_by_id('ENST00000260943')
chrom = 'chr' + transcript.contig.removeprefix('chr')
ranges = sorted(transcript.coding_sequence_position_ranges)
positions = [p for start, end in ranges for p in range(start, end + 1)]
if transcript.strand == '-':
    positions.reverse()
cds = ''.join(fasta.fetch(chrom, start - 1, end).upper() for start, end in ranges)
if transcript.strand == '-':
    cds = cds.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
codons = dict(zip((''.join(x) for x in itertools.product('TCAG', repeat=3)),
                  'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
protein = transcript.protein_sequence
if protein.startswith('X'):
    # A partial first CDS codon is represented as X in this transcript's protein.
    # Match the complete downstream codons in all three phases, without inventing
    # the missing nucleotide(s) or treating that partial codon as complete.
    matches = []
    for phase in range(3):
        translated = ''.join(codons.get(cds[i:i + 3], 'X') for i in range(phase, len(cds) - 2, 3))
        if translated[:112] == protein[1:113]:
            matches.append(phase)
    assert len(matches) == 1, (matches, protein[:113])
    first_complete_codon_offset = matches[0]
    offset = first_complete_codon_offset + 3 * (113 - 2)
else:
    translated = ''.join(codons.get(cds[i:i + 3], 'X') for i in range(0, len(cds) - 2, 3))
    assert protein[:113] == translated[:113], (protein[:113], translated[:113])
    first_complete_codon_offset = 0
    offset = 3 * (113 - 1)
assert protein[112] == 'D'
codon = cds[offset:offset + 3]
assert codons[codon] == 'D'
alleles = []
for i in range(3):
    for alt in 'ACGT':
        if alt != codon[i] and codons[codon[:i] + alt + codon[i + 1:]] == 'V':
            ref = codon[i]
            if transcript.strand == '-':
                ref = ref.translate(str.maketrans('ACGT', 'TGCA'))
                alt = alt.translate(str.maketrans('ACGT', 'TGCA'))
            alleles.append([chrom, positions[offset + i], ref, alt])
assert alleles == [['chr2', 211947513, 'T', 'A']], alleles
result = dict(status='VERIFIED', transcript_id=transcript.id, gene=transcript.gene_name,
              reported_change='D113V', cDNA_coding_sequence_available=transcript.coding_sequence is not None,
              CDS_prefix_matches_protein_through_reported_position=True,
              partial_first_amino_acid=protein.startswith('X'),
              first_complete_codon_offset=first_complete_codon_offset,
              reference_codon=codon, codon_positions=positions[offset:offset + 3], compatible_alleles=alleles,
              interpretation='Unavailable cDNA string resolved using GTF CDS positions, recorded reference genome and the matching protein prefix; same allele as the two cDNA-mapped transcripts',
              worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(sys.argv[1]).write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2), flush=True)
