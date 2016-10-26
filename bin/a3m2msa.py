#!/usr/bin/env python
docstring='''
a3m2msa.py hhblits.a3m hhblits.msa
    convert a3m format hhblits output file "hhblits.a3m" to FASTA format file
    "hhblits.msa". Assuming that the first sequence is the query sequence
'''
import sys

def hhblits2msa(hhblits_a3m=""):
    '''convert HHsuite a3m format output text into FASTA format text'''
    hhblits_msa=''
    sequence=''
    qlen=0
    for block in hhblits_a3m.split('>')[1:]:
        if not block.strip():
            continue
        Hit_id=block.split()[0]
        if '|' in Hit_id:
            Hit_id=Hit_id.split('|')[1]
        Hsp_hseq=''.join(block.splitlines()[1:])
        Hsp_hseq=''.join([res for res in Hsp_hseq if not ('a'<=res and res<='z')])
        if not sequence:
            sequence=Hsp_hseq
            qlen=len(sequence)
        else:
            aln_len=sum([not res in ".-" for res in Hsp_hseq])
            Hsp_identity=sum([q==t for q,t in zip(sequence,Hsp_hseq)])
            hhblits_msa+=">%s\t%d/%d\n%s\n"%(Hit_id,Hsp_identity,aln_len,Hsp_hseq)
    return hhblits_msa

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    fp=open(sys.argv[1],'rU')
    hhblits_a3m=fp.read()
    fp.close()
    hhblits_msa=hhblits2msa(hhblits_a3m)

    if len(sys.argv)==2:
        sys.stdout.write(hhblits_msa)
    else:
        fp=open(sys.argv[2],'w')
        fp.write(hhblits_msa)
        fp.close()
