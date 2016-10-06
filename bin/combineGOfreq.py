#!/usr/bin/env python
docstring='''
combineGOfreq.py datadir
    Combine GOfreq score from blast and psiblast output.

input:
    psiblast_localID_{MF,BP,CC}
    psiblast_GOfreq_{MF,BP,CC}
    blastp_GOfreq_{MF,BP,CC}

output:
    combine_GOfreq_{MF,BP,CC}.{2,1,0.5}

MF,BP,CC for molecular function, biological process, cellular component.
Result of blastp-GOfreq and psiblast-GOfreq are combined by:
    GOfreq=(seqID^t)*GOfreq_blastp+(1-seqID^t)*GOfreq_psiblast
    where t=2,1,0.5
    seqID is the maximum local seqID for psiblast hit
'''
import sys,os

t_list=[1.0]
seqID_cutoff=0 # seqID below which only PSI-BLAST is considered

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

def combineGOfreq(seqID_file="psiblast_localID_MF",
    psiblast_file="psiblast_GOfreq_MF",blast_file="blastp_GOfreq_BP"):
    '''combine PSIBLAST-GOfreq and BLAST-GOfreq score.
    return a dict whose k is parameter t and value is combined GOfreq score
    '''
    GOfreq_pred_dict=dict()

    seqID_dict=parseGOfreq(seqID_file)
    psiblast_dict=parseGOfreq(psiblast_file)
    blastp_dict=parseGOfreq(blastp_file)

    seqID=max([0]+seqID_dict.values())

    for t in t_list:
        if seqID<=seqID_cutoff:
            fp=open(psiblast_file)
            GOfreq_pred_dict[t]=fp.read()
            fp.close()
            continue

        GOfreq_pred_dict[t]=''
        for GOterm in sorted(psiblast_dict):
            if not GOterm in blastp_dict:
                blastp_dict[GOterm]=0
            cscore=(seqID**t)*  blastp_dict[GOterm]+ \
                 (1-seqID**t)*psiblast_dict[GOterm]
            GOfreq_pred_dict[t]+="%s\t%.2f\n"%(GOterm,cscore)
    return GOfreq_pred_dict

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()
    
    for datadir in sys.argv[1:]:
        for Aspect in ["MF","BP","CC"]:
            seqID_file=os.path.join(datadir,"psiblast_localID_"+Aspect)
            #seqID_file=os.path.join(datadir,"blastp_localID_"+Aspect)
            psiblast_file=os.path.join(datadir,"psiblast_GOfreq_"+Aspect)
            blastp_file  =os.path.join(datadir,"blastp_GOfreq_"+Aspect)

            GOfreq_pred_dict=combineGOfreq(
                seqID_file,psiblast_file,blastp_file)

            for t in GOfreq_pred_dict:
                fp=open(os.path.join(datadir,
                    "combine_GOfreq_%s.%.1f"%(Aspect,t)),'w')
                fp.write(GOfreq_pred_dict[t])
                fp.close()
