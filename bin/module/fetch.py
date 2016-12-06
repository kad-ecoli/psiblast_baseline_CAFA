#!/usr/bin/env python
import os,sys
import urllib

def wget(url='',outfile='',no_err=False,show_url=False):
    '''retrieve file from internet if not exists at current directory.
    return file output filename if download successfully.
    return empty string otherwise.

    outfile - output file name. By default it is the basename
    no_err - whether supress downloading error
    show_url - whether print url at stdout
    '''
    if not outfile:
        outfile=os.path.basename(url)
    if not os.path.isfile(outfile):
        if show_url:
            sys.stderr.write("fetching %s\n"%url)
        try:
            urllib.urlretrieve(url, outfile)
        except Exception,err:
            if not err:
                sys.stderr.write(str(err)+'\n')
            return ''
    elif show_url:
        sys.stderr.write("%s already exists\n"%outfile)
    return outfile
