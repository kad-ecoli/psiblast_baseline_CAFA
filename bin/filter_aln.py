#!/usr/bin/env python
docstring='''
filter_aln.py seq.fasta blastp.msa homoflag blastp.msa.homoflag
    filter out sequences in FASTA MSA file "blastp.msa" that shares
    sequence similarity with sequence in single FASTA file "seq.fasta".
    output filtered MSA to "blastp.msa.homoflag"

    if homoflag is set to real, preserve all sequence
    if homoflag is set to benchmark, remove sequence sharing >=0.3
        sequence identity with "seq.fasta"
    if homoflag is set to a number, remove sequence sharing at least
        specified sequence idenity with "seq.fasta"
'''
import sys

def read_single_fasta(fasta_file):
    '''read single entry FASTA format sequence file "fasta_file"
    and return the sequence'''
    fp=open(fasta_file,'rU')
    txt=fp.read().strip()
    fp.close()
    if txt.startswith('>'): # fasta format:
        sequence=''.join([line.strip() for line in txt.splitlines()[1:]])
    else: # plain text
        sequence=''.join([line.strip() for line in txt.splitlines()])
    return sequence

def read_multiple_fasta(fasta_file):
    '''read multiple entry FASTA format sequence "fasta_file".
    return a list for headers and a list for sequences
    '''
    header_list=[]
    sequence_list=[]
    fp=open(fasta_file,'rU')
    txt=fp.read().strip().lstrip('>')
    fp.close()
    if not txt:
        return header_list,sequence_list
    for block in txt.split('\n>'):
        lines=block.splitlines()
        header_list.append(lines[0])
        sequence_list.append(''.join([l.strip() for l in lines[1:]]))
    return header_list,sequence_list

def filter_aln(sequence,msa_file,idcut=1):
    '''filter FATSA format MSA file "msa_file" using sequence identity
    cutoff "idcut" '''
    msa_txt=''
    header_list,sequence_list=read_multiple_fasta(msa_file)
    if idcut==1:
        return ''.join(['>'+h+'\n'+s+'\n' for (h,s) in \
            zip(header_list,sequence_list)])

    for h,s in zip(header_list,sequence_list):
        identical_residue_num=len([a for a,b in zip(s,sequence) if a==b])
        seqID=1.*identical_residue_num/len([a for a in s if a!='-'])
        if seqID<idcut:
            msa_txt+='>'+h+'\n'+s+'\n'
    return msa_txt

if __name__=="__main__":
    if len(sys.argv)<3:
        sys.stderr.write(docstring)
        exit()

    homoflag=sys.argv[3] if len(sys.argv)>3 else "real"
    if homoflag=="real":
        idcut=1
    elif homoflag=="benchmark":
        idcut=0.3
    else:
        idcut=float(homoflag)

    sequence=read_single_fasta(sys.argv[1])
    msa_txt=filter_aln(sequence,sys.argv[2],idcut)

    if len(sys.argv)<=4:
        sys.stdout.write(msa_txt)
    else:
        fp=open(sys.argv[4],'w')
        fp.write(msa_txt)
        fp.close()
