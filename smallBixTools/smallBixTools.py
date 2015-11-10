from itertools import groupby
from operator import itemgetter
import sys, traceback, os
from Bio import SeqIO


def dct_to_fasta(d, fn):
    """
    :param d: dictionary in the form: {sequence_id: sequence_string, id_2: sequence_2, etc.}
    :param fn: The file name to write the fasta formatted file to.
    :return: Returns True if successfully wrote to file.
    """
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
    """
    :param fn: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
             {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
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


def hamdist(str1, str2):
    """
    Use this after aligning sequences.
    This counts the number of differences between equal length str1 and str2
    The order of the input sequences does not matter.
    :param str1: The first sequence.
    :param str2: The second sequence.
    :return: Returns a float value of the number of differences divided by the length of the first input argument.
    """
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs +=1
    return float(diffs)/float(len(str1))


def find_ranges(data):
    """
    Find contiguous ranges in a list of numerical values.
    eg: data = [1,2,3,4,8,9,10]
        find_ranges(data) will return:
        [[1, 2, 3, 4], [8, 9, 10]]
    :param data: a list of numerical values. (eg: int, float, long, complex)
    :return: a list of lists, each is a contiguous list of values.
    """
    ranges = []
    for k, g in groupby(enumerate(data), lambda i_x:i_x[0]-i_x[1]):
        ranges.append(list(map(itemgetter(1), g)))
    for rng in ranges:
        if len(rng) == 1:
            ranges.remove(rng)
    return ranges


def get_regions_from_panel(in_fn, regions, wd, outfn):
    """
    Slices regions out of a fasta formatted file, joins them together, and writes the resulting fasta file to the given location.
    an example call might be: get_regions_from_panel("test.fasta", [[0, 10], [20, 30]], "/tmp", "outfile.fasta")
    which would, for each sequence in the input file: "test.fasta", take the region from 0 to 10 joined with the
    region from 20 to 30, and write the result to the file: "/tmp/outfile.fasta".
    :param in_fn: the source / input fasta formatted file.
    :param regions: a list of lists. each sub-list has a start and a stop value. these demote the "regions" to
    use / slice. eg: [[0, 10], [20, 30]].
    :param wd: the directory where the output file will be written to.
    :param outfn: the output file name.
    :return: no return.
    """
    p_dct = fasta_to_dct(in_fn)
    fw = open(os.path.join(wd, outfn), "w")
    for k, v in list(p_dct.items()):
        p_seq = v
        p_joined = ""
        for rgn in regions:
            p_joined += p_seq[rgn[0]:rgn[1]]
        fw.write(">"+k+"\n"+p_joined+"\n")
    fw.close()


def get_parent(tree, child_clade):
    """
    Not used. removing in next commit.
    :param tree:
    :param child_clade:
    :return:
    """
    node_path = tree.get_path(child_clade)
    return node_path[-2]


def main():
    print("Call to main in smallBixTools.py. Nothing to do in the main.")
