#!/usr/bin/env python
docstring='''
combineGOfreq.py datadir
    Combine GOfreq score from blast and psiblast output.

input:
    psiblast_localID_{MF,BP,CC}
    psiblast_GOfreq_{MF,BP,CC}
    blastp_GOfreq_{MF,BP,CC}

output:
    combine_GOfreq_{MF,BP,CC}

MF,BP,CC for molecular function, biological process, cellular component.
Result of blastp-GOfreq and psiblast-GOfreq are combined by:
    GOfreq=(seqID^t)*GOfreq_blastp+(1-seqID^t)*GOfreq_psiblast
    where t=1 for MF and BP and t=0.5 for CC
    seqID is the maximum local seqID for psiblast hit
'''
import sys,os

t_dict={'MF':0.5,'BP':0.5,'CC':1.0}

def parseGOfreq(GOfreq_file="psiblast_GOfreq_MF"):
    '''parse GOfreq prediction data file. return a dict whose key is 
    GOterm and value is cscore
    '''
    GOfreq_dict=dict()
    fp=open(GOfreq_file,'rU')
    txt=fp.read().strip()
    fp.close()
    for line in txt.splitlines():
        line=line.split()
        GOfreq_dict['\t'.join(line[:-1])]=float(line[2])
    return GOfreq_dict

def combineGOfreq(t=1,seqID_file="psiblast_localID_MF",
    psiblast_file="psiblast_GOfreq_MF",blast_file="blastp_GOfreq_BP"):
    '''combine PSIBLAST-GOfreq and BLAST-GOfreq score.'''
    GOfreq_pred=''

    seqID_dict=parseGOfreq(seqID_file)
    psiblast_dict=parseGOfreq(psiblast_file)
    blastp_dict=parseGOfreq(blastp_file)

    seqID=max([0]+seqID_dict.values())

    for GOterm in sorted(psiblast_dict):
        if not GOterm in blastp_dict:
            blastp_dict[GOterm]=0
        cscore=(seqID**t)*  blastp_dict[GOterm]+ \
             (1-seqID**t)*psiblast_dict[GOterm]
        GOfreq_pred+="%s\t%.2f\n"%(GOterm,cscore)
    return GOfreq_pred

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()
    
    for datadir in sys.argv[1:]:
        for Aspect in ["MF","BP","CC"]:
            seqID_file=os.path.join(datadir,"psiblast_localID_"+Aspect)
            psiblast_file=os.path.join(datadir,"psiblast_GOfreq_"+Aspect)
            blastp_file  =os.path.join(datadir,"blastp_GOfreq_"+Aspect)

            GOfreq_pred=combineGOfreq(t_dict[Aspect],
                seqID_file,psiblast_file,blastp_file)

            fp=open(os.path.join(datadir,"combine_GOfreq_"+Aspect),'w')
            fp.write(GOfreq_pred)
            fp.close()
