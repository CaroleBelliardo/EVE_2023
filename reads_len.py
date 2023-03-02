#!/usr/bin/env python
# coding: utf-8

"""
  Usage:
    readLen.py -s <FILE> -b <REPO> -l <REPO>

  Options:
    -s, --species <FILE>   species name
    -b, --bed <REPO>   gene list
    -l, --lenRepo <REPO>    tree repositoy
"""

import csv
import os
import sys
import pandas as pd
from docopt import docopt

pd.options.display.float_format = '{:,.0f}'.format


def main1():
    df = pd.DataFrame()
    args = docopt(__doc__)

    with open(args['--species'], 'r') as file:
        liste_species = [line.strip() for line in file]

    for species in liste_species:
        specie = species.split('_')[0]
        condition = species.split('_')[1]

        bedF = args['--bed'] + '/' + species + '.bed'
        bed = pd.read_table(bedF, sep='\t', names=['Contig', 'start', 'stop', 'strand', 'id', 'eval', 'v', 'family'])

        for contig in bed.Contig:

            eveName = contig + '__' + str(
                bed.start[bed.Contig == contig].values[0]) + '__' + str(
                bed.stop[bed.Contig == contig].values[0])

            lenC = args['--lenRepo'] + '/' + specie + '/' + eveName + '_' + condition + '.bam.fastq_len'

            if os.path.exists(lenC) and os.path.getsize(lenC) > 0:
                with open(lenC, 'r') as file:
                    subDf = pd.DataFrame()
                    lenFile = [int(line.strip()) for line in file]
                    nbTotal = len(lenFile)
                    nb_si = len([i for i in lenFile if (i >= 21) and (i <= 23)])
                    nb_pi = len([i for i in lenFile if (i >= 25) and (i <= 32)])
                    subDf['species'] = [specie]
                    subDf['condition'] = [condition]
                    subDf['eveName'] = [eveName]
                    subDf['nbTotal'] = [nbTotal]
                    subDf['nb_si'] = [nb_si]
                    subDf['nb_pi'] = [nb_pi]

                    df = pd.concat([df, subDf], ignore_index=True)

    df.to_csv('len_table1.txt', sep='\t')


def main2():
    len_tab = pd.read_table('len_table1.txt', sep='\t')

    len_tab['type'] = 'none'
    len_tab.type[len_tab.nb_si < len_tab.nb_pi] = 'piRNA'
    len_tab.type[len_tab.nb_pi < len_tab.nb_si] = 'siRNA'

    data = pd.DataFrame({'count': len_tab.groupby(["species", "type"]).size()}).reset_index()

    df_synth = pd.DataFrame()
    df_synth['species'] = []
    df_synth['type'] = []

    for condition in ['FG', 'FS', 'MS', 'MG']:
        subT = len_tab.loc[len_tab['condition'] == condition]
        df_res = pd.DataFrame({condition: subT.groupby(["species", "type"]).size()}).reset_index()
        df_synth = df_synth.merge(df_res, how='outer', on=['species', 'type'])
    df_synth.to_csv('len_reads_FG_FS_MS_MG.txt', sep='\t')

    data = pd.DataFrame({'count': len_tab.groupby(["species", "type"]).size()}).reset_index()
    df_synth = pd.DataFrame()
    df_synth['species'] = []
    df_synth['type'] = []

    for condition in ['FG', 'FS']:
        subT = len_tab.loc[len_tab['condition'] == condition]
        df_res = pd.DataFrame({condition: subT.groupby(["species", "type"]).size()}).reset_index()
        df_synth = df_synth.merge(df_res, how='outer', on=['species', 'type'])

    df_synth.to_csv('len_reads_FG_FS.txt', sep='\t')

    df_synth = pd.DataFrame()
    df_synth['species'] = []
    df_synth['condition'] = []

    for type in ['piRNA', 'siRNA']:
        subT = len_tab.loc[len_tab['type'] == type]
        df_res = pd.DataFrame({type: subT.groupby(["species", "condition"]).size()}).reset_index()
        df_synth = df_synth.merge(df_res, how='outer', on=['species', 'condition'])

    #t = pd.DataFramtotal_SiPie({'': len_tab.groupby(["species", 'condition']).size()}).reset_index()

    #df_synth = df_synth.merge(t, how='outer', on=['species', 'condition'])

    df_synth['total_SiPi']= df_synth['piRNA'] + df_synth['siRNA']
    df_synth['piRNA_Percent'] = df_synth['piRNA'] / df_synth['total_SiPi'] * 100
    df_synth['siRNA_Percent'] = df_synth['siRNA'] / df_synth['total_SiPi'] * 100
    df_synth.to_csv('len_reads_FG_FS.txt2', sep='\t')

    df_synth = pd.DataFrame()
    df_synth['species'] = []

    for type in ['piRNA', 'siRNA']:
        subT = len_tab.loc[len_tab['type'] == type]
        df_res = pd.DataFrame({type: subT.groupby(["species"]).size()}).reset_index()
        df_synth = df_synth.merge(df_res, how='outer', on=['species'])
    print(df_synth)
    t = pd.DataFrame({'total_SiPi': len_tab.groupby(["species"]).size()}).reset_index()
    df_synth = df_synth.merge(t, how='outer', on=['species'])

    df_synth['piRNA_Percent'] = df_synth['piRNA'] /  df_synth['total_SiPi'] *100
    df_synth['siRNA_Percent'] = df_synth['siRNA'] /  df_synth['total_SiPi'] *100
    print(df_synth)
    df_synth.to_csv('len_reads_FG_FS_Mean.txt2', sep='\t')

    import seaborn as sns
    dp=pd.DataFrame()
    dp['Species']=df_synth['species']
    dp['type'] = 'siRNA'
    dp['Percent'] = df_synth.siRNA_Percent

    dp2=pd.DataFrame()
    dp2['Species']=df_synth['species']
    dp2['type'] = 'piRNA'
    dp2['Percent'] = df_synth.piRNA_Percent

    dp=dp.merge(dp2,how='outer')

    dp.to_csv('Plot_len_reads_FG_FS.txt', sep='\t')


if __name__ == '__main__':
    main1()
    main2()
