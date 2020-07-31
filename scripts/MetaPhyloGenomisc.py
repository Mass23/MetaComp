# usage: python3 MetaPhyloGenomisc.py -f K00163/K00163.fasta -k K00163 -g META_metadata.csv -i 0.3

def GetKEGGdb(KEGG_id):
    print('Importing KEGG reference sequences...')
    import urllib3
    http = urllib3.PoolManager()
    KEGG_url = http.request('GET','http://rest.kegg.jp/get/'+KEGG_id+'/')
    KEGG_list = str(KEGG_url.data).split('GENES')[1].split('\\n///\\n')[0].split('\\n')
    seq_dict = dict()
    for gene_line in KEGG_list:
        gene_line_split = gene_line.split(' ')
        gene_line_split = [i for i in gene_line_split if i != '']
        species_name = gene_line_split[0].replace(':','')
        gene_ids = gene_line_split[1:]
        for gene_id in gene_ids:
            print('http://rest.kegg.jp/get/' + species_name.lower() + ':' + gene_id.split('(')[0] + '         ', end = '\r')
            gene_url = http.request('GET','http://rest.kegg.jp/get/' + species_name.lower() + ':' + gene_id.split('(')[0])
            gene_page = str(gene_url.data)
            gene_sequence = gene_page.split('\\nNTSEQ')[1]
            gene_sequence = ''.join(gene_sequence.replace("\\n///\\n'",'').split('\\n            ')[1:]).upper()
            seq_dict[str(species_name + ':' + gene_id)] = gene_sequence
    print('Kegg reference sequences imported!')
    return(seq_dict)
def CreatKEGGdbFasta(sequence_dict):
    print('Creating KEGG db reference fasta file...')
    with open(kegg_accession + '/' + kegg_accession + '_db.fasta','w') as out_fasta:
        for sequence in sequence_dict.keys():
            out_fasta.write('>' + sequence + '\n' + sequence_dict[sequence] + '\n')
    return(kegg_accession + '/' + kegg_accession + '_db.fasta')
def FastaFilter1(input_file, prop_filter, longest_ref_seq):
    from Bio import SeqIO
    parsed_file = SeqIO.parse(input_file, 'fasta')
    output_file = [record for record in parsed_file if len(record.seq) > (longest_ref_seq * prop_filter)]
    SeqIO.write(output_file, input_file.replace('.fasta','_clean1.fasta'), 'fasta')
    return(input_file.replace('.fasta','_clean1.fasta'))
def AlignMacse(input_file, KEGG_db_fasta):
    import subprocess
    aln_args = ['java -jar macse.jar -prog alignSequences -seq', input_file, '-seq_lr', KEGG_db_fasta]
    subprocess.call(' '.join(aln_args),shell=True)
    output_file = input_file.replace('.fasta','_NT.fasta')
    return(output_file)
def TrimAtStop(sequence):
    trimmed_seq = ''
    for pos in range(0,len(sequence),3):
        codon = sequence[pos:pos+3]
        if codon in ['TAA','TGA','TGA']:
            break
        else:
            trimmed_seq = trimmed_seq + codon
    return(trimmed_seq)
def GetGC(sequence):
    GC_content = (sequence.count('C') + sequence.count('G')) / (sequence.count('A') + sequence.count('C') + sequence.count('G') + sequence.count('T'))
    return(GC_content)
def Main():
    import argparse
    import statistics
    # arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--Fasta", help="Fasta file to process, produced by GFF2Fasta.py.")
    parser.add_argument("-k", "--KeggAccession", help="KEGG id for the reference sequences.")
    parser.add_argument("-g", "--Group", help="File with the group IDs.")
    parser.add_argument("-l", "--LengthFilter", help="Minimal sequence length for the first filtering: relative proportion of the median length in the KEGG db of this accession.")
    parser.add_argument("-i", "--IdentityFilter", help="Minimal sequence identity that a queried sequence should share with one KEGG db sequence.")
    args = parser.parse_args()
    fasta_file = args.Fasta
    kegg_accession = args.KeggAccession
    variable_file = args.Group
    length_filter = args.LengthFilter
    identity_filter = args.IdentityFilter
    # Step 1: create KEGG sequences db using the KEGG API and create a fasta file
    KEGG_db_sequences = GetKEGGdb(kegg_accession)
    longest_db_seq = statistics.median([len(KEGG_db_sequences[seq]) for seq in KEGG_db_sequences.keys()])
    KEGG_db_fasta = CreatKEGGdbFasta(KEGG_db_sequences)
    # Step 2: open the fasta file and remove the sequences too short (clean1)
    fasta_file_clean1 = FastaFilter1(fasta_file, length_filter, longest_db_seq)
    # Step 3: Align the sequences and the database
    alignment_file = AlignMacse(fasta_file_clean1, KEGG_db_fasta)
    # Step 4: second cleaning of the sequences (clean2)
    # Step 5: compare the MG sequences to the db ones
Main()
