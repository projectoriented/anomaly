import pandas as pd

# this file is manually created from rna samples sheets doc
meta = pd.read_table(snakemake.input['rna_meta_tsv'], dtype=str).fillna("")
bam_paths = snakemake.input['fexist']


def get_path(rna_id):
    for path in bam_paths:
        if rna_id in path:
            return path


d = {
    'RNA_ID': [],
    'RNA_BAM_FILE': [],
    'DNA_VCF_FILE': [],
    'DNA_ID': [],
    'DROP_GROUP': [],
    'PAIRED_END': [],
    'COUNT_MODE': ['IntersectionStrict'],
    'COUNT_OVERLAPS': ['True'],
    'STRAND': [],
    'HPO_TERMS': [],
    'GENE_COUNTS_FILE': ['No'],
    'ANNOTATION': [],
    'TISSUE': [],
    'GENDER': [],
    'PHENOTYPE': [],
}

for row in meta.itertuples(index=False):
    rna_path = get_path(row.rna_id)
    d['RNA_BAM_FILE'].append(rna_path)

    d['RNA_ID'].append(row.rna_id)
    d['PAIRED_END'].append(row.paired_end)
    d['STRAND'].append(row.strand)
    d['PHENOTYPE'].append(row.phenotype)
    d['GENDER'].append(row.gender)
    d['TISSUE'].append(row.tissue)

    # specify drop groups by tissue
    d['DROP_GROUP'].append(row.tissue)

# padding the table to make it into a pandas object.
repeat_num = 0
for key, value in d.items():
    test = len(value)
    if test > 1:
        repeat_num = test
    elif test == 1:
        d[key] = value * repeat_num
    else:
        d[key] = [''] * repeat_num

# load dictionary into pandas object
df = pd.DataFrame(data=d)

# write out into tabulated form
df.to_csv(snakemake.output[0], sep='\t', index=False)
