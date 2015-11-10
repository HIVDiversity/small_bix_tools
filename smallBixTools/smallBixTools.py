def dct_to_fasta(d, fn):
    import re, sys, traceback, os
    fileName, fileExtension = os.path.splitext(fn)
    try:
        assert fileExtension.lower() in [".fasta", ".fa", ".fas", ".fna", ".ffn", ".faa", ".frn"]
    except AssertionError:
        _, _, tb = sys.exc_info()
        traceback.print_tb(tb) # Fixed format
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]
        print(('An error occurred on line {} in statement {}'.format(line, text)))
        exit(1)
    try:
        with open(fn, "w") as fw:
            for k, v in list(d.items()):
                fw.write(">"+k+"\n"+v+"\n")
        return True
    except Exception as e:
        print(e)
        return False

def fasta_to_dct(fn):
    from Bio import SeqIO
    import re, sys, traceback, os
    fileName, fileExtension = os.path.splitext(fn)
    try:
        assert fileExtension.lower() in [".fasta", ".fa", ".fas", ".fna", ".ffn", ".faa", ".frn"]
    except AssertionError:
        _, _, tb = sys.exc_info()
        traceback.print_tb(tb) # Fixed format
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]
        print(('An error occurred on line {} in statement {}'.format(line, text)))
        exit(1)
    dct = {}
    for sequence in SeqIO.parse(open(fn), "fasta"):
        dct[sequence.description.replace(" ", "_")] = str(sequence.seq)
    return dct

def hamdist(str1, str2):        #use after aligning the seqs
#counts the number of differences btwn equal length str1 and str2
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs +=1
    return float(diffs)/float(len(str1))

from itertools import groupby
from operator import itemgetter
def find_ranges(data):
    ranges = []
    for k, g in groupby(enumerate(data), lambda i_x:i_x[0]-i_x[1]):
        ranges.append(list(map(itemgetter(1), g)))
    for rng in ranges:
        if len(rng) == 1:
            ranges.remove(rng)
    return ranges

def get_regions_from_panel(pnl, regions, wd, outfn):
    p_dct = fasta_to_dct(pnl)
    fw = open(wd + "/"+outfn, "w")
    for k, v in list(p_dct.items()):
        p_seq = v
        p_joined = ""
        for rgn in regions:
            p_joined += p_seq[rgn[0]:rgn[1]]
        fw.write(">"+k+"\n"+p_joined+"\n")
    fw.close()

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    return node_path[-2]

def main():
    print("in main of smallBixTools.py")
