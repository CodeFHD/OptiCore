import re

from .surfaces import *

def _split_opened_zmx(f):
    lines = f.readlines()
    headerlines = []
    surflines = {} # dict, for each surface index, a list of the lines
    footerlines = []
    lc = 0 # linecounter

    # loop over header
    for line in lines:
        line = line.strip()
        if line.startswith('SURF'):
            break
        headerlines.append(line)
        lc += 1

    # loop over surfaces
    surfidx = -1
    for line in lines[lc:]:
        if re.match(r'\w', line) and not line.startswith('SURF'):
            # break when something else than a surf comes up
            break
        line = line.strip()
        lc += 1
        if line.startswith('SURF'):
            if surfidx > -1:
                surflines[surfidx] = curlines 
            surfidx = int(line.split()[1])
            curlines = []
        else:
            curlines.append(line)
    # finish last surface
    surflines[surfidx] = curlines 

    # loop over footer
    for line in lines[lc:]:
        line = line.strip()
        footerlines.append(line)

    return headerlines, surflines, footerlines

def split_zmx_file(filename):
    # check encoding of file
    f = open(filename, 'rb')
    data = f.read(2)
    if data == b'\xff\xfe':
        encoding = 'utf-16-le'
    else:
        encoding = 'utf-8'
    f.close()

    # split he(ader), su(rfaces), fo(oter)
    with open (filename, 'r', encoding=encoding) as f:    
        he, su, fo = _split_opened_zmx(f)

    return he, su, fo