#!/usr/bin/env python
docstring='''cscore2csv.py combine_GOfreq_MF combine_GOfreq_MF.csv
    clean GOfreq format or COFACTOR format GO prediction file by:
    [1] mapping alternaive GOterm ID to current GOterm
    [2] removing root of 3 Aspect and "protein binding"
    [3] Convert Cscore into two decimal digits
    [4] Remove all prediction with zero confidence
    [5] removing redundant predictions
    [6] rank GOterm by Cscore in descending order

option:
    -infmt={GOfreq,COFACTOR,GoFDR} input format
        GOfreq: GOterm Aspect Cscore
        COFACTOR: GOterm Cscore Name
        GoFDR: Target GOterm Cscore
    -outfmt={default,CAFA} output format
        default: GOterm Aspect Cscore Name
        CAFA: GOterm Cscore
    -target='' target name
        If not empty, it will show at the beginning of every line 
'''
import sys
from module import obo2csv
from module.fetch import wget
obo_url="http://geneontology.org/ontology/go-basic.obo"

def cscore2csv(GOfreq_txt='',obo_dict=dict(),infmt="GOfreq",outfmt="default",
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575",
    target=''):
    '''parse GO predictions and output to CSV format'''
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    
    report_dict=dict()
    for line in GOfreq_txt.strip().splitlines():
        if infmt=="GOfreq":
            GOterm,Aspect,Cscore=line.split()
        elif infmt=="COFACTOR":
            GOterm,Cscore=line.split()[:2]
        elif infmt=="GoFDR":
            GOterm,Cscore=line.split()[-2:]
        else:
            sys.stderr.write("ERROR! Unknown format %s\n"%infmt)

        ## [1] mapping alternaive GOterm ID to current GOterm ##
        GOterm=obo_dict.alt_id(GOterm)
        if not GOterm:
            continue # obsolete GO term

        ## [2] remove root of 3 aspects and "protein binding"
        if GOterm in excludeGO:
            continue

        ## [3] Convert Cscore into two decimal digits
        Cscore=float("%.2f"%float(Cscore))

        ## [4] Remove all prediction with zero confidence
        if Cscore<=0:
            continue

        for Aspect in obo_dict:
            if GOterm in obo_dict[Aspect]["Term"]:
                if not Aspect in report_dict:
                    report_dict[Aspect]=dict()
                report_dict[Aspect][GOterm]=Cscore
    
    report_txt=''
    for Aspect in report_dict:
        ## [5] removing redundant predictions
        for GOterm in report_dict[Aspect].keys():
            if not GOterm in report_dict[Aspect]:
                continue # Term already removed
            parent_node_set=set(obo_dict.is_a(GOterm,direct=False,name=False,
                number=False).strip().split())-set([GOterm])
            for parent_GOterm in list(parent_node_set):
                if parent_GOterm in report_dict[Aspect] and report_dict[
                    Aspect][parent_GOterm]<=report_dict[Aspect][GOterm]:
                    del report_dict[Aspect][parent_GOterm]

        ## [6] rank GOterm by Cscore in descending order
        report_list=sorted([(report_dict[Aspect][GOterm],GOterm)
            for GOterm in report_dict[Aspect]],reverse=True)
    
        for Cscore,GOterm in report_list:
            if outfmt=="default":
                Name=obo_dict[Aspect]["Term"][GOterm].name
                line='%s\t%s\t%.2f\t%s\n'%(GOterm,Aspect,Cscore,Name)
            elif outfmt=="CAFA":
                line="%s\t%.2f\n"%(GOterm,Cscore)
            else:
                sys.stderr.write("ERROR! Unknown format %s\n"%outfmt)
            report_txt+=(target+'\t')*(target!='')+line
    return report_txt

if __name__=="__main__":
    infmt="GOfreq"
    outfmt="default"
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575"
    target=''

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-infmt="):
            infmt=arg[len("-infmt="):]
        elif arg.startswith("-outfmt="):
            outfmt=arg[len("-outfmt="):]
        elif arg.startswith("-excludeGO="):
            excludeGO=arg[len("-excludeGO="):]
        elif arg.startswith("-target="):
            target=arg[len("-target="):]
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown format %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    #### parse GO hierachy ####
    fp=open(wget(obo_url,show_url=True),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)

    #### parse GO prediction result ####
    fp=open(argv[0],'rU')
    GOfreq_txt=fp.read()
    fp.close()

    report_txt=cscore2csv(GOfreq_txt,obo_dict,infmt,outfmt,excludeGO,target)

    #### write output ####
    fp=sys.stdout
    if len(argv)==2:
        fp=open(sys.argv[1],'w')
    fp.write(report_txt)
    if len(argv)==2:
        fp.close()
