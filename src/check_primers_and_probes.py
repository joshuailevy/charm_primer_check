import pandas as pd
from Bio import SeqIO
import regex as re
from Bio.Seq import Seq
import gzip
import shutil
from itertools import product
from pathlib import Path


def edit_distance(seq_a, seq_b):
    """Return the Levenshtein edit distance between two sequences."""
    seq_a = str(seq_a).upper()
    seq_b = str(seq_b).upper()

    previous_row = list(range(len(seq_b) + 1))
    for i, base_a in enumerate(seq_a, start=1):
        current_row = [i]
        for j, base_b in enumerate(seq_b, start=1):
            insertion = current_row[j - 1] + 1
            deletion = previous_row[j] + 1
            substitution = previous_row[j - 1] + (base_a != base_b)
            current_row.append(min(insertion, deletion, substitution))
        previous_row = current_row

    return previous_row[-1]


def expand_primer_options(primer, translator):
    primer = str(primer).upper()
    base_options = [translator.get(base, [base]) for base in primer]
    return ["".join(option) for option in product(*base_options)]


def gunzip_if_needed(path):
    #this will gunzip the multi fasta file if needed
    path = Path(path)
    if path.suffix != ".gz":
        return str(path)

    unzipped_path = path.with_suffix("")
    if (
        not unzipped_path.exists()
        or unzipped_path.stat().st_mtime < path.stat().st_mtime
    ):
        with gzip.open(path, "rb") as zipped_file:
            with open(unzipped_path, "wb") as unzipped_file:
                shutil.copyfileobj(zipped_file, unzipped_file)

    return str(unzipped_path)


genomes = gunzip_if_needed('../sequences/Parainfluenza_virus_1/Parainfluenza_virus_1.fasta.gz')
metadata_path = '../sequences/Parainfluenza_virus_1/Parainfluenza_virus_1_metadata.csv'
# old_fp = "CGAGCGAGTCGATTTATTACCA"
forward_primer = "CGGGCGAGYMGATTTATTACCA"
reverse_primer = "CAATCCGGTTAACATAATTTGT"
probe = "AATATGGCATTAAAAGARGCAGGW"

ambiguous_bases = {'R', 'Y', 'S', 'W',
                    'K', 'M', 'B', 'D',
                    'H', 'V', 'N'}
## expanding ambiguous bases to multiple primers in the ambiguous set
translator = {'Y':['C','T'],'W':['A','T'],'M':['A','C'],'R':['A','G']}

primer_df = pd.DataFrame(
    [
        [primer_left, primer_right]
        for primer_left in expand_primer_options(forward_primer, translator)
        for primer_right in expand_primer_options(reverse_primer, translator)
    ],
    columns=["primer_seq_x", "primer_seq_y"],
)

metadata_df = pd.read_csv(
    metadata_path,
    usecols=["accession", "collection_date", "geo_loc_name"],
).set_index("accession")

maxmismatch=3
all_results = []

for _, row in primer_df.iterrows():
    primer_left = row["primer_seq_x"]
    primer_right = row["primer_seq_y"]
    primer_right = str(Seq(primer_right).reverse_complement())
    
    pattern_left = f"({primer_left}){{s<={maxmismatch}}}"
    pattern_right = f"({primer_right}){{s<={maxmismatch}}}"

    # # Parse sequences in the multifasta
    for genome_record in SeqIO.parse(genomes, "fasta"):
        genome_id = genome_record.id
        genome_seq = str(genome_record.seq)
        genome_metadata = metadata_df.loc[genome_id]
        collection_date = genome_metadata["collection_date"]
        geo_loc_name = genome_metadata["geo_loc_name"]

        ## string matching for now - should switch to something more robust. 
        left_fwd = [m.start() for m in re.finditer(pattern_left,
                                                   genome_seq,
                                                   flags=re.IGNORECASE,
                                                   overlapped=True)]
        right_rev = [m.start() for m in re.finditer(pattern_right,
                                                    genome_seq,
                                                    flags=re.IGNORECASE,
                                                    overlapped=True)]
        if len(left_fwd)>1 or len(right_rev)>1:
            print('multiple matching sites!')
            asdfasdf
        left_fwd_actual = [genome_seq[pos:pos+len(primer_left)]
                           for pos in left_fwd]
        right_rev_actual = [genome_seq[pos:pos+len(primer_right)]
                            for pos in right_rev]
        distance_left = [edit_distance(actual, primer_left)
                         for actual in left_fwd_actual]
        distance_right = [edit_distance(actual, primer_right)
                          for actual in right_rev_actual]

        left_has_ambiguity = [any(base in ambiguous_bases
                                  for base in lfa.upper())
                              for lfa in left_fwd_actual]
        right_has_ambiguity = [any(base in ambiguous_bases
                                   for base in rfa.upper())
                               for rfa in right_rev_actual]
        for jL,lfa in enumerate(left_fwd_actual):
            for jR,rfa in enumerate(right_rev_actual):
                all_results.append([genome_id,collection_date,geo_loc_name,primer_left,primer_right,lfa,rfa,distance_left[jL],left_has_ambiguity[jL],distance_right[jR],right_has_ambiguity[jR]])

result_df = pd.DataFrame(all_results, columns=['sequence','collection_date','geo_loc_name','primer_left','primer_right_rev_comp','left_match','right_rev_comp_match','edit_dist_left','left_amb_bases','edit_dist_right','right_amb_bases'])

left_dist_min = result_df.groupby('sequence')['edit_dist_left'].min()
right_dist_min = result_df.groupby('sequence')['edit_dist_right'].min()

print('Minimum edit distances (all sequences), for left/right primers')
print(left_dist_min.value_counts().sort_index())
print(right_dist_min.value_counts().sort_index())

result_df = result_df[result_df['collection_date'].str.contains('2024')| 
                        result_df['collection_date'].str.contains('2025')| 
                        result_df['collection_date'].str.contains('2026')]
left_dist_min = result_df.groupby('sequence')['edit_dist_left'].min()
right_dist_min = result_df.groupby('sequence')['edit_dist_right'].min()

print('Minimum edit distances (sequences from last two years), for left/right primers')
print(left_dist_min.value_counts().sort_index())
print(right_dist_min.value_counts().sort_index())


probe_results = []
### now do the same for the probe
probe_df = pd.DataFrame([probe_v for probe_v in expand_primer_options(probe, translator)], columns=["probe_seq"])
# # Parse sequences in the multifasta
for _, row in probe_df.iterrows():
    probe_seq = row["probe_seq"]
    # primer_right = row["primer_seq_y"]
    # primer_right = str(Seq(primer_right).reverse_complement())
    pattern_probe = f"({probe_seq}){{s<={maxmismatch}}}"
    # pattern_right = f"({primer_right}){{s<={maxmismatch}}}"
    for genome_record in SeqIO.parse(genomes, "fasta"):
        genome_id = genome_record.id
        genome_seq = str(genome_record.seq)
        genome_metadata = metadata_df.loc[genome_id]
        collection_date = genome_metadata["collection_date"]
        geo_loc_name = genome_metadata["geo_loc_name"]

        ## string matching for now - should switch to something more robust. 
        probe_fwd = [m.start() for m in re.finditer(pattern_probe,
                                                   genome_seq,
                                                   flags=re.IGNORECASE,
                                                   overlapped=True)]
        if len(probe_fwd)>1 :
            print('multiple matching sites!')
            asdfasdf
        probe_fwd_actual = [genome_seq[pos:pos+len(probe)]
                           for pos in probe_fwd]
        distance_probe = [edit_distance(actual, probe_seq)
                         for actual in probe_fwd_actual]

        site_has_ambiguity = [any(base in ambiguous_bases
                                  for base in lfa.upper())
                              for lfa in probe_fwd_actual]
        for jL,pfa in enumerate(probe_fwd_actual):
            probe_results.append([genome_id,collection_date,geo_loc_name,probe_seq,pfa,distance_probe[jL],site_has_ambiguity[jL]])
probe_result_df = pd.DataFrame(probe_results, columns=['sequence','collection_date','geo_loc_name','probe_seq','probe_match','edit_distance_probe','probe_amb_bases'])
probe_dist_min = probe_result_df.groupby('sequence')['edit_distance_probe'].min()

print('Minimum edit distances (all seqs), for probe')
print(probe_dist_min.value_counts().sort_index())

probe_result_df = probe_result_df[probe_result_df['collection_date'].str.contains('2024')| 
                        probe_result_df['collection_date'].str.contains('2025')| 
                        probe_result_df['collection_date'].str.contains('2026')]
probe_dist_min = probe_result_df.groupby('sequence')['edit_distance_probe'].min()

print('Minimum edit distances (sequences from last two years), for probe')
print(probe_dist_min.value_counts().sort_index())
