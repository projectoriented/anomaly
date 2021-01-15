import pandas as pd

samples = pd.read_table(snakemake.input[0], dtype=str)
d = {
    'RNA_ID': [],
    'RNA_BAM_FILE': [],
    'DNA_VCF_FILE': [],
    'DNA_ID': [],
    'DROP_GROUP': [],
    'PAIRED_END': ['True'],
    'COUNT_MODE': ['IntersectionStrict'],
    'COUNT_OVERLAPS': ['True'],
    'STRAND': ['Reverse'],
    'HPO_TERMS': [],
    'GENE_COUNTS_FILE': ['No'],
    'ANNOTATION': [],
    'TISSUE': [],
    'SEX': [],
    'PHENOTYPE': [],
}


def remove_string_rna(base):
    substring = '-RNA'
    if substring in base:
        return base.replace(substring, '')
    else:
        return base


for row in samples.itertuples(index=False):
    basename_search = f"{row.lane}_{row.sample}_{row.sample_number}"
    abs_path = (f"/home/mei.wu/rna-seq/output/{row.sample}/picard_markdupe/{basename_search}_trim_star_marked.bam")
    d['RNA_ID'].append(basename_search)
    d['RNA_BAM_FILE'].append(abs_path)

    modified = remove_string_rna(row.sample)[-2:]
    d['TISSUE'].append('Fibroblast') if row.sample == '164-I-2A' or row.sample == '90-I-1A' else d['TISSUE'].append('Blood')
    d['SEX'].append('FEMALE') if '2' in modified else d['SEX'].append('MALE')
    d['PHENOTYPE'].append('UNAFFECTED') if 'U' in modified else d['PHENOTYPE'].append('AFFECTED')

d['DROP_GROUP'] = d['TISSUE']

repeat_num = 0
for key, value in d.items():
    test = len(value)
    if test > 1:
        repeat_num = test
    elif test == 1:
        d[key] = value * repeat_num
    else:
        d[key] = [''] * repeat_num


df = pd.DataFrame(data=d)
df.to_csv(snakemake.output[0], sep='\t', index=False)
