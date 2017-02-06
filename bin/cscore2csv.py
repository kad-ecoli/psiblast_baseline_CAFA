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
    -T={eps,ps,svg,png} plot directed acyclic directed graph (DAG) by
        graphviz in the specified format. default is do not plot DAG.
        The output image file will be named after the input file.
    -execpath=dot
        path to graphviz's "dot" executable. default is at the same folder
        of this script.
'''
import sys,os
from module import obo2csv
from module.fetch import wget
import subprocess
import textwrap

obo_url="http://geneontology.org/ontology/go-basic.obo"
color_list=["red","orange","yellow","green","cyan","magenta","white"]

def cscore2csv(GOfreq_txt='',obo_dict=dict(),infmt="GOfreq",outfmt="default",
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575",
    target='',min_cscore=0):
    '''parse GO predictions and output to CSV format

    options:
        infmt - (input format)
            GOfreq: (default) GOterm Aspect Cscore
            COFACTOR: GOterm Cscore Name
            GoFDR: Target GOterm Cscore
        outfmt - output format
            default: GOterm Aspect Cscore Name
            CAFA: GOterm Cscore
        excludeGO - list of GO terms to be excluded.
            default is to exclude root terms and protein binding
        target - target name
        min_cscore - minimum cscore below which terms are not considered
    '''
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
        if Cscore<=min_cscore:
            continue

        for Aspect in obo_dict:
            if GOterm in obo_dict[Aspect]["Term"]:
                if not Aspect in report_dict:
                    report_dict[Aspect]=dict()
                report_dict[Aspect][GOterm]=Cscore
    
    report_txt=''
    GOterm_dict=dict() # a dict whose key is 
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
        if not report_dict[Aspect]:
            del report_dict[Aspect]
            continue

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
    return report_txt,report_dict

def detect_graphviz(execpath=''):
    '''locate the path of graphviz's "dot" executable'''
    if execpath:
        execpath=os.path.abspath(execpath)
    elif os.path.isfile(os.path.join(os.path.dirname(__file__),"dot")):
        execpath=os.path.join(os.path.dirname(__file__),"dot")
    if not os.path.isfile(execpath):
        execpath="dot"
    return execpath

def make_GOterm_node(GOterm,obo_dict,cscore=0):
    '''make node in graphviz graph for a GO term'''
    color=color_list[min([int((1-cscore)*10),len(color_list)-1])]
    GVtxt='"%s"[label="%s" shape=rectangle fillcolor=%s style=filled];'%(
        GOterm,GOterm+'\n'+textwrap.fill(obo_dict.Term(GOterm).name,25),color)
    return GVtxt

def draw_GO_DAG(report_dict,dag_img,obo_dict,T,execpath,
    cscore_min=0.3,color_term_list=[]):
    '''plot directed acyclic graph listed for GO terms listed in "report_dict".
    write the image file to "dag_img". use "T" as output format. Use graphviz
    executable "execpath" for image compilation.

    option:
        cscore_min - minimum cscore to be plotted, default is 0.3
        color_term_list - list of explicit GO terms to be colored
    '''
    GOterm_list=[]
    is_a_list=[]
    GVtxt='digraph G{ graph[splines=true,rankdir="BT"];'
    for Aspect in report_dict:
        for GOterm in report_dict[Aspect]:
            if not GOterm in GOterm_list:
                cscore=report_dict[Aspect][GOterm]
                if cscore<cscore_min:
                    continue
                if not color_term_list or GOterm in color_term_list:
                    GVtxt+=make_GOterm_node(GOterm,obo_dict,cscore)
                else:
                    GVtxt+=make_GOterm_node(GOterm,obo_dict)
                GOterm_list.append(GOterm)

            # set of parent nodes whose is_a relation is plotted,
            parent_node_set=set(obo_dict.is_a(GOterm,direct=False,
                name=False,number=False).strip().split())-set([GOterm])
            # list of parent nodes whose is_a relation is plotted
            plotted_parent_list=[GOterm]

            while(parent_node_set):
                for GOterm in list(plotted_parent_list):
                    direct_parent_node_set=set(obo_dict.is_a(GOterm,
                        direct=True,name=False,number=False).strip(
                        ).split())-set([GOterm])
                    parent_node_set-=direct_parent_node_set
                    plotted_parent_list+=list(direct_parent_node_set)
                    for direct_parent in list(direct_parent_node_set):
                        if not (GOterm,direct_parent) in is_a_list:
                            GVtxt+='"%s"->"%s";'%(GOterm,direct_parent)
                            is_a_list.append((GOterm,direct_parent))
                            if not direct_parent in GOterm_list:
                                if not direct_parent in color_term_list or \
                                   not direct_parent in report_dict[Aspect]:
                                    GVtxt+=make_GOterm_node(direct_parent,obo_dict)
                                else:
                                    GVtxt+=make_GOterm_node(direct_parent,obo_dict,
                                        report_dict[Aspect][direct_parent])
                                GOterm_list.append(GOterm)
    GVtxt+="}"

    p=subprocess.Popen("%s -T%s"%(execpath,T),shell=True,
        stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate(input=GVtxt)
    fp=open(dag_img,'w')
    fp.write(stdout)
    fp.close()
    return GVtxt

if __name__=="__main__":
    infmt="GOfreq"
    outfmt="default"
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575"
    T=""        # do not plot DAG
    execpath="" # path of graphviz executable
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
        elif arg.startswith("-T="):
            T=arg[len("-T="):]
        elif arg.startswith('-'):
            sys.stderr.write("ERROR! Unknown format %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<1:
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

    report_txt,report_dict=cscore2csv(
        GOfreq_txt,obo_dict,infmt,outfmt,excludeGO,target)

    #### write output ####
    fp=sys.stdout
    if len(argv)>1:
        fp=open(argv[1],'w')
    fp.write(report_txt)
    if len(argv)>1:
        fp.close()

    #### plot DAG using graphviz ####
    if T:
        execpath=detect_graphviz(execpath)
        dag_img=argv[0]+'.'+T
        if len(argv)>2:
            dag_img=argv[2]
        draw_GO_DAG(report_dict,dag_img,obo_dict,T,execpath)
