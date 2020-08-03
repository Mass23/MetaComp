import glob
import pandas as pd
import os
from Bio import SeqIO
import argparse
import pandas as pd

def ReverseSequence(sequence,direction):
    if direction == '-':
        sequence = sequence[::-1]
    return(sequence)

def FindKEGG(row):
    if 'KEGG' in row[8]:
        KEGG_annotation = row[8].split(';KEGG=')[1].split(';')[0]
    else:
        KEGG_annotation = 'Unassigned'
    return(KEGG_annotation)

def AddKEGGDict(KEGG_dict, KEGG_annotation):
    if KEGG_annotation in KEGG_dict.keys():
        KEGG_dict[KEGG_annotation] +=1
    else:
        KEGG_dict[KEGG_annotation] = 1
        if KEGG_annotation not in glob.glob('*'):
            os.mkdir(KEGG_annotation)
    return(KEGG_dict)

def AppendKEGGFasta(KEGG_annotation, sequence, gene_name):
    KEGG_file = open(KEGG_annotation + '/' + KEGG_annotation + '.fasta',"a")
    KEGG_file.write('>' + gene_name + '\n' + str(sequence) + '\n')
    KEGG_file.close()
    return(1)

def Main():
    print('MC parse:')
    print('  - Creating KOs fasta files')
    print('  - Counting KOs annotations\n')
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--InputFile", help="File containing the files to process", required = True)
    args = parser.parse_args()
    input_file = args.InputFile

    samples_file = pd.read_csv(input_file, sep = ',', header = None)
    output_table = pd.DataFrame()

    for index,row in samples_file.iterrows():
        sample = row[0]
        print('Sample: ' + sample + '                        ')
        sample_annotation = pd.read_csv(row[1], sep = '\t', header = None)
        sample_fasta = SeqIO.to_dict(SeqIO.parse(row[2],'fasta'))
        gene_count = 0
        sample_dict = dict()

        for index, row in sample_annotation.iterrows():
            gene_count += 1
            print('  - gene number: ' + str(gene_count), end = '\r')
            contig_name = row[0]
            gene_type = row[2]
            pos_beg = row[3]
            pos_end = row[4]
            direction = row[6]
            gene_name = sample + '_' + str(gene_count)

            sequence = ReverseSequence(sample_fasta[contig_name].seq[pos_beg:pos_end+1], direction)
            KEGG_annotation = FindKEGG(row)
            KEGG_dict = AddKEGGDict(sample_dict, KEGG_annotation)
            AppendKEGGFasta(KEGG_annotation, sequence, gene_name)

        output_table = pd.concat([output_table, pd.DataFrame({sample: [sample_dict[i] for i in sample_dict.keys()]}, index =[i for i in sample_dict.keys()])], axis=1, sort=False)
    output_table = output_table.fillna(0)
    output_table.to_csv('MC_KEGG_table.tsv', sep = '\t', header = True, index = True)

if __name__ == "__main__":
    Main()
