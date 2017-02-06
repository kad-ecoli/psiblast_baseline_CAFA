#!/usr/bin/env python
docstring='''
blast2msa.py query.fasta blastp.xml blastp.msa
    convert XML blastp outfile file "blastp.xml" to FASTA file "blastp.msa"
    aligned to query sequnce supplied by FASTA format query sequence 
    "query.fasta"
'''
import sys
import re
from filter_aln import read_single_fasta

# hit sequence name
Hit_id_pattern=re.compile("<Hit_id>([\w\W]+?)</Hit_id>")
# e-value
Hsp_evalue_pattern=re.compile("<Hsp_evalue>([-.e\d]+?)</Hsp_evalue>")
# first aligned residue in query
Hsp_query_from_pattern=re.compile("<Hsp_query\-from>(\d+)</Hsp_query\-from>")
# last aligned residue in query
Hsp_query_to_pattern=re.compile("<Hsp_query\-to>(\d+)</Hsp_query\-to>")
# number of identical residues
Hsp_identity_pattern=re.compile("<Hsp_identity>(\d+)</Hsp_identity>")
# aligned query sequence
Hsp_qseq_pattern=re.compile("<Hsp_qseq>([-\w]+?)</Hsp_qseq>")
# aligned hit sequence
Hsp_hseq_pattern=re.compile("<Hsp_hseq>([-\w]+?)</Hsp_hseq>")

def blast2msa(sequence,blastp_xml=""):
    '''convert ncbi blast+ XML format output text into FASTA format text
    where all blastp hits are aligned to query sequence "sequence"
    '''
    blastp_msa=''
    qlen=len(sequence)
    for block in blastp_xml.split("<Hit>")[1:]:
        Hit_id=Hit_id_pattern.findall(block)[0]
        Hsp_evalue=Hsp_evalue_pattern.findall(block)[0]
        Hsp_identity=Hsp_identity_pattern.findall(block)[0]
        Hsp_query_from=int(Hsp_query_from_pattern.findall(block)[0])
        Hsp_query_to=int(Hsp_query_to_pattern.findall(block)[0])
        Hsp_qseq=Hsp_qseq_pattern.findall(block)[0]
        Hsp_hseq=Hsp_hseq_pattern.findall(block)[0]

        aln_len=Hsp_query_to-Hsp_query_from+1
        Hsp_qseq,Hsp_hseq=zip(*[(q,h) for (q,h) in \
            zip(Hsp_qseq,Hsp_hseq) if q!='-'])
        Hsp_hseq='-'*(Hsp_query_from-1)+''.join(Hsp_hseq)+ \
                 '-'*(qlen-Hsp_query_to)
        
        header=Hit_id+'\t'+Hsp_evalue+'\t'+Hsp_identity+'/'+str(aln_len)
        blastp_msa+='>'+header+'\n'+Hsp_hseq+'\n'
    return blastp_msa

if __name__=="__main__":
    if len(sys.argv)<3:
        sys.stderr.write(docstring)
        exit()
    
    sequence=read_single_fasta(sys.argv[1])
    if sys.argv[2]=='-':
        blastp_xml=sys.stdin.read()
    else:
        fp=open(sys.argv[2],'rU')
        blastp_xml=fp.read()
        fp.close()
    blastp_msa=blast2msa(sequence,blastp_xml)

    if len(sys.argv)==3:
        sys.stdout.write(blastp_msa)
    else:
        fp=open(sys.argv[3],'w')
        fp.write(blastp_msa)
        fp.close()
