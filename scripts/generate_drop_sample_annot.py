import pandas as pd

samples = pd.read_table(snakemake.input[0], dtype=str)
meta = pd.read_table(snakemake.input[1], dtype=str).fillna("").set_index("sample_id") # this file is manually created from rna samples sheets doc

d = {
    'RNA_ID': [],
    'RNA_BAM_FILE': [],
    'DNA_VCF_FILE': [],
    'DNA_ID': [],
    'DROP_GROUP': [],
    'PAIRED_END': ['True'],
    'STRAND': ['Reverse'],
    'COUNT_MODE': ['IntersectionStrict'],
    'COUNT_OVERLAPS': ['True'],
    'HPO_TERMS': [],
    'GENE_COUNTS_FILE': ['No'],
    'ANNOTATION': [],
    # 'TISSUE': [],
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
    abs_path = (f"/home/mei.wu/rna-seq/drop_files/{row.sample}/{row.sample}_trim_star_marked_no_alt.bam")
    d['RNA_ID'].append(row.sample)
    d['RNA_BAM_FILE'].append(abs_path)

    modified = remove_string_rna(row.sample)

    dna_sample_id = meta.loc[remove_string_rna(row.sample)]['dna_sample_id']
    d['DNA_ID'].append(dna_sample_id)

    dna_case = meta.loc[remove_string_rna(row.sample)]['dna_case']

    if not row.sample.isalnum() and dna_case:
        dir_name = row.sample.split("-")[0]
        vcf_path = f"/home/proj/development/rare-disease/rna_data/{dir_name}/dna_vcf/{dna_case}_gatkcomb.vcf.gz"
    elif not dna_case:
        vcf_path = ""
    else:  # e.g. ACC3285A1 is alnum.
        dir_name = "15138"
        vcf_path = f"/home/proj/development/rare-disease/rna_data/{dir_name}/dna_vcf/{dna_case}_gatkcomb.vcf.gz"

    d['DNA_VCF_FILE'].append(vcf_path)

    d['DROP_GROUP'].append('Fibroblast') if row.sample == '164-I-2A' or row.sample == '90-I-1A' or '15001' in row.sample else d['DROP_GROUP'].append(
        'Blood')
    d['SEX'].append('FEMALE') if '2' in modified[-2:] else d['SEX'].append('MALE')
    d['PHENOTYPE'].append('UNAFFECTED') if 'U' in modified[-2:] else d['PHENOTYPE'].append('AFFECTED')

# d['DROP_GROUP'] = d['TISSUE']

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