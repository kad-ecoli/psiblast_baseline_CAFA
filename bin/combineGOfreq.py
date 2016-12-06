#!/usr/bin/env python
docstring='''
combineGOfreq.py datadir
    Combine GOfreq score from blast and psiblast output.

input:
    psiblast_localID_{MF,BP,CC}
    psiblast_GOfreq_{MF,BP,CC}
    blastp_GOfreq_{MF,BP,CC}

    hhblits_GOfreq_{MF,BP,CC}    (optional)
    hhblits_globalID_{MF,BP,CC}  (optional)

output:
    combine_GOfreq_{MF,BP,CC}

    hhbGOfreq_GOfreq_{MF,BP,CC}   (optional)
    hhbglobalID_GOfreq_{MF,BP,CC} (optional)

MF,BP,CC for molecular function, biological process, cellular component.
Result of blastp-GOfreq and psiblast-GOfreq are combined by:
    GOfreq=(seqID^t)*GOfreq_blastp+(1-seqID^t)*GOfreq_psiblast
    where t=1 for MF and BP and t=0.5 for CC
    seqID is the maximum local seqID for psiblast hit

    hhbGOfreq_GOfreq=1-(1-hhblits_GOfreq)(1-GOfreq)
    hhbseqID_GOfreq=1-(1-hhblits_globalID)(1-GOfreq)
'''
import sys,os

t_dict={'MF':1.0,'BP':1.0,'CC':1.0}

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

    seqID_dict=parseGOfreq(seqID_file)
    psiblast_dict=parseGOfreq(psiblast_file)
    blastp_dict=parseGOfreq(blastp_file)

    seqID=max([0]+seqID_dict.values())
    
    GOfreq_pred_list=[]
    for GOterm in sorted(psiblast_dict):
        if not GOterm in blastp_dict:
            blastp_dict[GOterm]=0
        cscore=(seqID**t)*  blastp_dict[GOterm]+ \
             (1-seqID**t)*psiblast_dict[GOterm]
        if float('%.2f'%cscore):
            GOfreq_pred_list.append((cscore,GOterm))
    GOfreq_pred=''.join(['%s\t%.2f\n'%(GOterm,cscore) for cscore,GOterm \
        in sorted(GOfreq_pred_list,reverse=True)])
    return GOfreq_pred,GOfreq_pred_list

def hhbseqID_GOfreq(GOfreq_pred_list=[],
    hhblits_seqID_file="hhblits_globalID_MF"):
    '''combine the combineGOfreq score in GOfreq_pred_list and HHBLITS-globalID 
    score'''
    seqID_dict=parseGOfreq(hhblits_seqID_file)
    
    GOfreq_dict=dict()
    for cscore,GOterm in GOfreq_pred_list:
        GOfreq_dict[GOterm]=cscore

    hhbseqID_GOfreq_pred_list=[]
    for GOterm in list(set(seqID_dict.keys()+GOfreq_dict.keys())):
        if not GOterm in seqID_dict:
            cscore=GOfreq_dict[GOterm]
        elif not GOterm in GOfreq_dict:
            cscore=seqID_dict[GOterm]
        else:
            cscore=1-(1-GOfreq_dict[GOterm])*(1-seqID_dict[GOterm])
        if float('%.2f'%cscore):
            hhbseqID_GOfreq_pred_list.append((cscore,GOterm))
    hhbseqID_GOfreq_pred=''.join(['%s\t%.2f\n'%(GOterm,cscore) for cscore,GOterm \
        in sorted(hhbseqID_GOfreq_pred_list,reverse=True)])
    return hhbseqID_GOfreq_pred

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()
    
    for datadir in sys.argv[1:]:
        for Aspect in ["MF","BP","CC"]:
            #seqID_file=os.path.join(datadir,"psiblast_localID_"+Aspect)
            #psiblast_file=os.path.join(datadir,"psiblast_GOfreq_"+Aspect)
            #blastp_file  =os.path.join(datadir,"blastp_GOfreq_"+Aspect)

            #GOfreq_pred,GOfreq_pred_list=combineGOfreq(t_dict[Aspect],
                #seqID_file,psiblast_file,blastp_file)
            #fp=open(os.path.join(datadir,"combine_GOfreq_"+Aspect),'w')
            #fp.write(GOfreq_pred)
            #fp.close()

            seqID_file=os.path.join(datadir,"psiblast_globalID_"+Aspect)
            psiblast_file=os.path.join(datadir,"psiblast_gwGOfreq_"+Aspect)
            blastp_file  =os.path.join(datadir,"blastp_gwGOfreq_"+Aspect)

            GOfreq_pred,GOfreq_pred_list=combineGOfreq(t_dict[Aspect],
                seqID_file,psiblast_file,blastp_file)

            fp=open(os.path.join(datadir,"combine_gwGOfreq_"+Aspect),'w')
            fp.write(GOfreq_pred)
            fp.close()

            continue # skip hhblits combination

            hhblits_seqID_file=os.path.join(datadir,"hhblits_globalID_"+Aspect)
            if os.path.isfile(hhblits_seqID_file):
                hhbseqID_GOfreq_pred=hhbseqID_GOfreq(
                    GOfreq_pred_list, hhblits_seqID_file)

                fp=open(os.path.join(datadir,"hhbseqID_GOfreq_"+Aspect),'w')
                fp.write(hhbseqID_GOfreq_pred)
                fp.close()

            hhblits_GOfreq_file=os.path.join(datadir,"hhblits_GOfreq_"+Aspect)
            if os.path.isfile(hhblits_seqID_file):
                hhbGOfreq_GOfreq_pred=hhbseqID_GOfreq(
                    GOfreq_pred_list, hhblits_GOfreq_file)

                fp=open(os.path.join(datadir,"hhbGOfreq_GOfreq_"+Aspect),'w')
                fp.write(hhbGOfreq_GOfreq_pred)
                fp.close()
