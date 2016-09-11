#!/usr/bin/env python
docstring='''
msa2go.py -MFdb=PDB_GOterms.MF -BPdb=PDB_GOterms.BP -CCdb=PDB_GOterms.CC -prefix=blastp_  seq.fasta blastp.msa 
    GO transfer using MSA file "blastp.msa" and COFACTOR format GO mapping
    file PDB_GOterms.*

options:
    -{MFdb,BPdb,CCdb} path to GO mapping file for MF, BP, CC
        if not set, do not predict that catergory
    -prefix file prefix of the output GO prediction file
        if -prefix=blastp_ then output GO-MF file will be blastp_MF

output:
    prefix_{globalID,localID,evalue,GOfreq}_{MF,BP,CC}

globalID,localID,GOfreq,evalue are different scoring scheme:
    
Target       AKLMNWTY
Template1    -KICNWT- GO:0000001 evalue=0.01
Template2    --ICNWT- GO:0000001 evalue=0.1
Template3    AKLCNYTY GO:0000002 evalue=0.001

Confidence Score of GO:0000001
[1] globalID (identical residues number normalized by target length): 
    Cscore(GO:0000001)=4/8
[2] localID (identical residue number normalized by aligned region length):
    Cscore(GO:0000001)=4/6
[3] GOfreq (number of templates annotated with GO divided by template number):
    Cscore(GO:0000001)=2/3
[4] evalue (Sigmoid transform of log evalue):
    Cscore(GO:0000001)=1-Sigmoid(x)=1-1/(1+exp(-x)),x=log(evalue*alpha)/log(alpha)
[5] gwGOfreq (globalID weighted GOfreq):
    Cscore(GO:0000001)=(4/8+3/8)/3
[6] lwGOfreq (localID weighted GOfreq):
    Cscore(GO:0000001)=(4/6+3/5)/3
[7] lwGOfreq (evalue weighted GOfreq):
    Cscore(GO:0000001)=(0.661+0.583)/3
'''
import sys
from filter_aln import read_single_fasta,read_multiple_fasta
from math import exp,log
alpha=10**3

def msa2seqID(GO_list,header_list,sequence):
    '''GO prediction by global/local sequence identity'''
    globalID_pred='' # CSV format GO prediction by global seqID
    localID_pred=''  # CSV format GO prediction by local seqID
    seqlen=len(sequence)

    GO_globalID_pred_dict=dict()
    GO_localID_pred_dict=dict()
    for line,header in zip(GO_list,header_list):
        if not line:
            continue
        seqID=header.strip().split()[-1]
        identical_residue_num=float(seqID.split('/')[0])
        aligned_residue_num=float(seqID.split('/')[1])
        globalID=identical_residue_num/seqlen
        localID=identical_residue_num/aligned_residue_num
        for GO in line:
            if GO in GO_globalID_pred_dict:
                if globalID>GO_globalID_pred_dict[GO]:
                    GO_globalID_pred_dict[GO]=globalID
                if localID>GO_localID_pred_dict[GO]:
                    GO_localID_pred_dict[GO]=localID
            else:
                GO_globalID_pred_dict[GO]=globalID
                GO_localID_pred_dict[GO]=localID
    
    for GO in GO_globalID_pred_dict:
        globalID_pred+="%s\t%.2f\n"%(GO,GO_globalID_pred_dict[GO])
        localID_pred+="%s\t%.2f\n"%(GO,GO_localID_pred_dict[GO])
    return globalID_pred,localID_pred

def scale_evalue(evalue=0):
    '''Sigmoid transform evalue so that it is bound between [0,1]'''
    cscore=1
    if evalue:  # the larger alpha is, the smaller cscore will become
        x=log(evalue*alpha)/log(alpha)
        cscore=1-1/(1+exp(-x))
    return cscore

def msa2evalue(GO_list,header_list,sequence):
    '''GO prediction by evalue. evalue are converted to cscore by
    cscore=1-1/(1+exp(-log(evalue)/log(1000)))
    '''
    evalue_pred='' # CSV format GO prediction by e-value

    evalue_pred_dict=dict()
    for line,header in zip(GO_list,header_list):
        if not line:
            continue
        evalue=float(header.split()[1])
        cscore=scale_evalue(evalue)
        for GO in line:
            if GO in evalue_pred_dict:
                if cscore>evalue_pred_dict[GO]:
                    evalue_pred_dict[GO]=cscore
            else:
                evalue_pred_dict[GO]=cscore
    
    for GO in evalue_pred_dict:
        evalue_pred+="%s\t%.2f\n"%(GO,evalue_pred_dict[GO])
    return evalue_pred

def msa2GOfreq(GO_list):
    '''GO prediction using the frequency of a GO term being annotated in all 
    templates'''
    GOfreq_pred='' # CSV format GO prediction
    
    annotated_template_num=0. # number of templates annotated by at least one GO
    GO_list_cat=[]
    for line in GO_list:
        GO_list_cat+=line
        if line:
            annotated_template_num+=1.
    GO_set=set(GO_list_cat)

    for GO in GO_set:
        GO_count=len([line for line in GO_list if GO in line])
        GO_freq=GO_count/annotated_template_num
        GOfreq_pred+="%s\t%.2f\n"%(GO,GO_freq)
    return GOfreq_pred

def msa2wGOfreq(GO_list,header_list,sequence):
    '''GO prediction by globalID/localID/evalue weighted GOfreq'''
    gwGOfreq_pred='' # CSV format GO prediction by global seqID weighted GOfreq
    lwGOfreq_pred='' # CSV format GO prediction by local seqID weighted GOfreq
    ewGOfreq_pred='' # CSV format GO prediction by evalue weighted GOfreq
    seqlen=len(sequence)

    GO_gwGOfreq_pred_dict=dict()
    GO_lwGOfreq_pred_dict=dict()
    GO_ewGOfreq_pred_dict=dict()
    for line,header in zip(GO_list,header_list):
        if not line:
            continue
        seqID=header.strip().split()[-1]
        identical_residue_num=float(seqID.split('/')[0])
        aligned_residue_num=float(seqID.split('/')[1])
        globalID=identical_residue_num/seqlen
        localID=identical_residue_num/aligned_residue_num
        evalue=float(header.split()[1])
        cscore=scale_evalue(evalue)
        for GO in line:
            if GO in GO_gwGOfreq_pred_dict:
                GO_gwGOfreq_pred_dict[GO]+=globalID
                GO_lwGOfreq_pred_dict[GO]+=localID
                GO_ewGOfreq_pred_dict[GO]+=cscore
            else:
                GO_gwGOfreq_pred_dict[GO]=globalID
                GO_lwGOfreq_pred_dict[GO]=localID
                GO_ewGOfreq_pred_dict[GO]=cscore

    annotated_template_num=1.*len([line for line in GO_list if line])

    for GO in GO_gwGOfreq_pred_dict:
        gwGOfreq_pred+="%s\t%.2f\n"%(GO,GO_gwGOfreq_pred_dict[GO]/annotated_template_num)
        lwGOfreq_pred+="%s\t%.2f\n"%(GO,GO_lwGOfreq_pred_dict[GO]/annotated_template_num)
        ewGOfreq_pred+="%s\t%.2f\n"%(GO,GO_ewGOfreq_pred_dict[GO]/annotated_template_num)
    return gwGOfreq_pred,lwGOfreq_pred,ewGOfreq_pred

def parse_GOdb(GOdb=''):
    '''parse COFACTOR format GO mapping file "GOdb"
    return a dict, whose key is protein name and value is a list of GO'''
    GOdict=dict()
    fp=open(GOdb,'rU')
    txt=fp.read()
    fp.close()
    for line in txt.splitlines():
        if not line.strip():
            continue
        name,GO_str=line.split()
        GOdict[name]=GO_str.split(',')
    return GOdict

def label_GO_to_template(GOdb='',header_list=[]):
    '''assign GO terms to sequence listed in "header_list" using COFACTOR
    format GO mapping file "GOdb"'''
    GOdict=parse_GOdb(GOdb)
    GO_list=[]
    for header in header_list:
        header=header.split()[0]
        if header in GOdict:
            GO_list.append(GOdict[header])
        else:
            GO_list.append([])
    return GO_list

if __name__=="__main__":
    MFdb=''
    BPdb=''
    CCdb=''
    prefix=''

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-MFdb="):
            MFdb=arg[len("-MFdb="):]
        elif arg.startswith("-BPdb="):
            BPdb=arg[len("-BPdb="):]
        elif arg.startswith("-CCdb="):
            CCdb=arg[len("-CCdb="):]
        elif arg.startswith("-prefix="):
            prefix=arg[len("-prefix="):]
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)
    GOdb_dict={"MF":MFdb,"BP":BPdb,"CC":CCdb}

    if len(argv)<2:
        sys.stderr.write(docstring)
        exit()

    sequence=read_single_fasta(argv[0])
    header_list,sequence_list=read_multiple_fasta(argv[1])

    for Aspect in GOdb_dict:
        GO_list=label_GO_to_template(GOdb_dict[Aspect],header_list)
        globalID_pred,localID_pred=msa2seqID(GO_list,header_list,sequence)
        GOfreq_pred=msa2GOfreq(GO_list)
        evalue_pred=msa2evalue(GO_list,header_list,sequence)
        gwGOfreq_pred,lwGOfreq_pred,ewGOfreq_pred=msa2wGOfreq(
            GO_list,header_list,sequence)

        fp=open(prefix+"globalID_"+Aspect,'w')
        fp.write(globalID_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()

        fp=open(prefix+"localID_"+Aspect,'w')
        fp.write(localID_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()

        fp=open(prefix+"GOfreq_"+Aspect,'w')
        fp.write(GOfreq_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()

        fp=open(prefix+"evalue_"+Aspect,'w')
        fp.write(evalue_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()

        fp=open(prefix+"gwGOfreq_"+Aspect,'w')
        fp.write(gwGOfreq_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()

        fp=open(prefix+"lwGOfreq_"+Aspect,'w')
        fp.write(lwGOfreq_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()

        fp=open(prefix+"ewGOfreq_"+Aspect,'w')
        fp.write(ewGOfreq_pred.replace("\t","\t%s\t"%Aspect[-1]))
        fp.close()
