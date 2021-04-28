#!/usr/bin/env python3

import yaml
import pandas as pd


sample_annot = pd.read_table(snakemake.input[0], dtype=str).fillna("")
cnt_dupes = sample_annot.pivot_table(index=['DROP_GROUP'], aggfunc='size').ge(30) # >= 30


drop_config = {
    'projectTitle': f'"{snakemake.config["proj_title"]}"',
    'root': f"{snakemake.config['drop_proj']}/output",
    'htmlOutputPath': f"{snakemake.config['drop_proj']}/output/html",
    'indexWithFolderName': True,
    'hpoFile': None,
    'sampleAnnotation': snakemake.input[0],
    'geneAnnotation': {'v34': snakemake.config['gtf']},
    'genomeAssembly': 'hg38',
    'exportCounts': {'geneAnnotations': ['v34'], 'excludeGroups': cnt_dupes[cnt_dupes].index.tolist()},
    'aberrantExpression':
        {
            'groups': cnt_dupes[cnt_dupes].index.tolist(),
            'fpkmCutoff': 1,
            'implementation': 'autoencoder',
            'padjCutoff': 0.05,
            'zScoreCutoff': 0,
            'maxTestedDimensionProportion': 3
        },
    'aberrantSplicing':
        {
            'groups': cnt_dupes[cnt_dupes].index.tolist(),
            'recount': False,
            'longRead': False,
            'keepNonStandardChrs': True,
            'filter': True,
            'minExpressionInOneSample': 20,
            'minDeltaPsi': 0.05,
            'implementation': 'PCA-BB-Decoder',
            'padjCutoff': 0.05,
            'zScoreCutoff': 0,
            'deltaPsiCutoff': 0.3,
            'maxTestedDimensionProportion': 6
        },
    'mae':
        {
            'groups': cnt_dupes[cnt_dupes].index.tolist(),
            'genome': snakemake.config['genome'],
            'gatkIgnoreHeaderCheck': True,
            'padjCutoff': 0.05,
            'allelicRatioCutoff': 0.8,
            'addAF': False,
            'maxAF': 0.001,
            'maxVarFreqCohort': 0.05,
            'qcVcf': snakemake.config['qcVCF'],
            'qcGroups': 'mae'
        },
    'tools':
        {
        'gatkCmd': 'gatk',
        'bcftoolsCmd': 'bcftools',
        'samtoolsCmd': 'samtools'
        }
}

with open(snakemake.output[0], 'w') as file:
    data = yaml.dump(drop_config, file, sort_keys=False)