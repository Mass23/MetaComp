def GetKEGGdb(KEGG_id):
    import random
    print('Importing KEGG reference sequences...')
    import urllib3
    http = urllib3.PoolManager()


    KEGG_url = http.request('GET','http://rest.kegg.jp/get/'+KEGG_id+'/')
    KEGG_list = str(KEGG_url.data).split('GENES')[1].split('\\n///\\n')[0].split('\\n            ')

    # Filter only sequences from META or organisms in the list:
    KEGG_list = random.sample(KEGG_list,20)

    seq_dict = dict()
    for gene_line in KEGG_list:
        gene_line_split = gene_line.split(' ')
        gene_line_split = [i for i in gene_line_split if i != '']
        species_name = gene_line_split[0].replace(':','')
        gene_ids = gene_line_split[1:]
        gene_id = random.sample(gene_ids, 1)[0]

        gene_url = http.request('GET','http://rest.kegg.jp/get/' + species_name.lower() + ':' + gene_id.split('(')[0])
        gene_page = str(gene_url.data)
        print('http://rest.kegg.jp/get/' + species_name.lower() + ':' + gene_id.split('(')[0] + '         ', end = '\r')
        gene_sequence = gene_page.split('NTSEQ')[1]
        gene_sequence = ''.join(gene_sequence.replace("\\n///\\n'",'').split('\\n            ')[1:]).upper()
        seq_dict[str(species_name + ':' + gene_id)] = gene_sequence
    print('Kegg reference sequences imported!\n')
    return(seq_dict)

def RefSeqSpecies(species_name):
    refseq_list = []
    if species_name in refseq_list:
        return(True)
    else:
        return(False)

def CreatKEGGdbFasta(kegg_accession, sequence_dict):
    print('Creating KEGG db reference fasta file...')
    with open(kegg_accession + '/' + kegg_accession + '_db.fasta','w') as out_fasta:
        for sequence in sequence_dict.keys():
            out_fasta.write('>' + sequence + '\n' + sequence_dict[sequence] + '\n')
    return(kegg_accession + '/' + kegg_accession + '_db.fasta')

def FastaFilter1(kegg_accession, prop_filter, median_ref_seq):
    print('First sequence filtering...')
    from Bio import SeqIO
    input_file = kegg_accession + '/' + kegg_accession + '.fasta'
    parsed_file = SeqIO.parse(input_file, 'fasta')
    output_file = [record for record in parsed_file if len(record.seq) > (median_ref_seq * prop_filter)]
    n_kept_seq = len(output_file)
    if n_kept_seq > 0:
        SeqIO.write(output_file, input_file.replace('.fasta','_clean1.fasta'), 'fasta')
        return(input_file.replace('.fasta','_clean1.fasta'))
    else:
        return(False)

def AlignMacse(input_file, KEGG_db_fasta):
    print('Alignment:')
    import subprocess
    aln_args = ['java -jar macse_v2.03.jar -prog alignSequences -seq', input_file, '-seq_lr', KEGG_db_fasta]
    subprocess.call(' '.join(aln_args),shell=True)
    output_file = input_file.replace('.fasta','_NT.fasta')
    return(output_file)

def IncompleteCodon(codon):
    if '-' in codon:
        return(True)
    else:
        return(False)

def StopShiftCodon(codon):
    if (codon in ['TAA','TGA','TGA']) or ('!' in codon):
        return(True)
    else:
        return(False)

def TrimSequence(sequence):
    #Trims all nucleotides outside of codons and stop the sequence at stop codons
    trimmed_seq = ''
    stop_shift_codon = False
    for pos in range(0,len(sequence),3):
        codon = sequence[pos:pos+3]

        # Is the gene already finished?
        if stop_shift_codon == True:
            trimmed_seq = trimmed_seq + '---'
        # Is the codon incomplete? Does it contains a stop or frameshift?
        else:
            incomplete = IncompleteCodon(codon)
            if incomplete == True:
                trimmed_seq = trimmed_seq + '---'
            else:
                stopshift = StopShiftCodon(codon)
                if stopshift == True:
                    stop_shift_codon = True
                    trimmed_seq = trimmed_seq + '---'
    return(trimmed_seq)

def FastaFilter2(kegg_accession, prop_filter, median_ref_seq):
    print('Second sequence filtering...')
    from Bio import SeqIO
    input_file = kegg_accession + '/' + kegg_accession + '_clean1_NT.fasta'
    parsed_file = SeqIO.parse(input_file, 'fasta')

    # filtering
    for record in parsed_file:
        record.seq = TrimSequence(record.seq)

    output_file = [record for record in parsed_file if len(record.seq != '-') > (median_ref_seq * prop_filter)]

    n_kept_seq = len(output_file)
    if n_kept_seq > 0:
        SeqIO.write(output_file, input_file.replace('.fasta','_clean2.fasta'), 'fasta')
        return(input_file.replace('.fasta','_clean2.fasta'))
    else:
        return(False)

def GetGC(sequence):
    GC_content = (sequence.count('C') + sequence.count('G')) / (sequence.count('A') + sequence.count('C') + sequence.count('G') + sequence.count('T'))
    return(GC_content)

def Main():
    print('MC align:')
    print('  - Aligning sequences to the KEGG GENES db sequences')
    print('  - Filtering for evolutionary analyses\n')

    import argparse
    import statistics

    # arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--KeggAccession", help="KEGG id for the reference sequences.", required = True)
    parser.add_argument("-l", "--LengthFilter", help="Minimal sequence length for the first filtering: relative proportion of the median length in the KEGG db of this accession.", default = 0.5)

    args = parser.parse_args()
    kegg_accession = args.KeggAccession
    length_filter = float(args.LengthFilter)

    # Step 1: create KEGG sequences db using the KEGG API and create a fasta file
    KEGG_db_sequences = GetKEGGdb(kegg_accession)
    median_db_seq = statistics.median([len(KEGG_db_sequences[seq]) for seq in KEGG_db_sequences.keys()])
    KEGG_db_fasta = CreatKEGGdbFasta(kegg_accession, KEGG_db_sequences)

    # Step 2: open the fasta file and remove the sequences too short (clean1)
    fasta_file_clean1 = FastaFilter1(kegg_accession, length_filter, median_db_seq)
    if fasta_file_clean1 == False:
        print('No sequences passed the length threshold...')
        return(0)

    # Step 3: Align the sequences and the database
    alignment_file = AlignMacse(fasta_file_clean1, KEGG_db_fasta)

    # Step 4: second cleaning of the sequences (clean2)
    fasta_file_clean2 = FastaFilter2(kegg_accession, length_filter, median_db_seq)
    if fasta_file_clean2 == False:
        print('No sequences passed the length threshold...')
        return(0)

if __name__ == "__main__":
    Main()
