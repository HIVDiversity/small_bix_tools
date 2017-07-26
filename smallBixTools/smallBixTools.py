from itertools import groupby
from operator import itemgetter
import sys, traceback, os
from Bio import SeqIO
import operator


def cluster_counter_from_fasta(orig_fasta, cluster_fn):
    print("The source fasta file which was used for clustering: {}".format(orig_fasta))
    print("The .uc cluster file produced from software like vsearch: {}".format(cluster_fn))


def convert_count_to_frequency_on_fasta(source_fasta_fn, target_fasta_fn):
    """
    when running vsearch as such:
    vsearch --cluster_fast {} --id 0.97 --sizeout --centroids {}
    We get a centroids.fasta file with seqid header lines like:
    >ATTCCGGTATCT_9;size=1432;
    >CATCATCGTAAG_14;size=1;
    etc.
    This method converts those count values into frequencies.
    Notes: The delimiter between sections in the sequence id must be ";".
    There must be a section in the sequence id which has exactly: "size=x" where x is an integer.
    This must be surrounded by ";"'s
    :param source_fasta_fn: the input fasta file. Full path required.
    :param target_fasta_fn: the output fasta file. Full path required. If this is the same as the input, the input will
                            be over-written.
    :return: No return value.
    """
    dct = fasta_to_dct(source_fasta_fn)
    total_count = 0
    for k in dct.keys():  # step over every sequence id, and tally up the total number of sequences in the orig. file.
        splt = k.split(";")
        for item in splt:
            if "size" in item:
                total_count += int(item.split("=")[1])
    dct2 = {}
    for k in dct.keys():  # step over each seqid, converting its count to a frequency, saving into a new dictionary.
        splt = k.split(";")
        new_k = ""
        for item in splt:
            if item == '':
                continue
            if "size" in item:
                new_k += "freq=" + str(int(item.split("=")[1]) / total_count) + ";"
            else:
                new_k += item + ";"
        new_k.replace(";;", ";")
        dct2[new_k] = dct[k]
    dct_to_fasta(dct2, target_fasta_fn)  # Write the new dictionary to the same filename which was fed into the
    print("Converted absolute counts to frequencies on file: {}".format(target_fasta_fn))


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)


def py2_fasta_iter(fasta_name):
    """
    from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def size_selector(in_fn, out_fn, min, max):
    dct = fasta_to_dct(in_fn)
    dct2 = {}
    for k, v in dct.items():
        if len(v) < min:
            continue
        if len(v)>max:
            dct2[k] = v[:max]
        else:
            dct2[k] = v
    dct_to_fasta(dct2, out_fn)


def split_file_into_timepoints(infile):
    # if we have one fasta file with multiple time points, and we want to split this into multiple files, each with
    # only one time point.
    # The sequence ids in the fasta file have to follow the naming schema:
    # >CAP177_1100_002wpi_C2C3_xxxx (where xxxx is variable)
    print("TODO.splitting file into multiple time points.TODO")


def own_cons_maker(infile):
    dct = fasta_to_dct(infile)
    seqs = [v for k,v in dct.items()]
    cons = ""
    for i in range(len(seqs[0])):
        pos = {}
        for k in range(len(seqs)):
            letter = seqs[k][i]
            if letter in pos.keys():
                pos[letter] += 1
            else:
                pos[letter] = 1
        cons += max(pos.iteritems(), key=operator.itemgetter(1))[0]
    return cons


def build_cons_seq(infile):
    # https://www.biostars.org/p/14026/
    from Bio import AlignIO
    from Bio.Align import AlignInfo

    alignment = AlignIO.read(open(infile), "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    return consensus


def auto_duplicate_removal(in_fn, out_fn):
    print("Trying to automatically remove duplicate sequences from the in file: %s" %in_fn)
    try:
        dct = {}
        for sequence in SeqIO.parse(open(in_fn), "fasta"):
            new_key = sequence.description.replace(" ", "_")
            new_seq = str(sequence.seq)

            if new_key in dct.keys():
                print("Duplicate sequence ids found. seq_id: %s" %new_key)
                print("Checking if the sequences are the same...")
                if new_seq == dct[new_key]:
                    print("The sequences are the same, the keys are the same. We update the record for this duplicate.")
                elif new_seq != dct[new_key]:
                    print("The sequences are NOT the same, but they have the same sequence id. Cannot auto resolve.")
                    print("Exiting...")
                    sys.exit()
            dct[new_key] = str(sequence.seq)
        dct_to_fasta(dct, out_fn)
        return True
    except Exception as e:
        print(e)
        print("Failed to auto remove duplicates")
        raise

def hyphen_to_underscore_fasta(fn, out_fn):
    print("Cleaning a fasta file to have underscores in the sequence names, not hyphens.")
    dct = fasta_to_dct(fn)
    cleaned_dct = {}
    for key, value in dct.items():
        cleaned_key = key.replace("-", "_")
        cleaned_dct[cleaned_key] = value
    dct_to_fasta(cleaned_dct, out_fn)
    print("Finished cleaning. Wrote to file: %s" %out_fn)


def customdist(s1, s2):
    assert len(s1) == len(s2)
    dist = 0
    for c1, c2 in zip(s1, s2):
        if c1 != c2:
            dist += 1
    diff = 0
    for i in range(len(s1)-1):
        if s1[i] != s2[i]:
            if (s1[i] == "-" and s1[i+1] == "-" and s2[i+1] != "-") or (s2[i] == "-" and s2[i+1] == "-" and s1[i+1] != "-"):
                diff += 1
    return (dist-diff)


def find_duplicate_ids(fn):
    # and squash if safe
    dct = {}
    for seq in SeqIO.parse(open(fn), "fasta"):
        seq_id = str(seq.description)
        if seq_id in dct.keys():
            dct[seq_id] += 1
        else:
            dct[seq_id] = 1
    for k, v in dct.items():
        if v != 1:
            print("%s was found %s times." %(k, v))
#    dct2 = {}
#    safe = []
#    for k, v in dct.items():
#        if v != 1:
#            seqs = []
#            for seq in SeqIO.parse(open(fn), "fasta"):
#                if str(seq.description) == k:
#                    seqs.append(str(seq.seq).replace("-", ""))
#            if len(list(set(seqs))) == 1:
#                print k
#                print "safe to squash..."
#                safe.append(k)
#    final_dct = {}
#    for seq in SeqIO.parse(open(fn), "fasta"):
#        seq_id = str(seq.description)
#        if seq_id in final_dct.keys():
#            if seq_id in safe:
#                final_dct[seq_id] = str(seq.seq)
#            else:
#                print(seq_id)
#                raw_input("rawr")
#        else:
#            final_dct[seq_id] = str(seq.seq)
#    return final_dct


def dct_to_fasta(d, fn):
    """
    :param d: dictionary in the form: {sequence_id: sequence_string, id_2: sequence_2, etc.}
    :param fn: The file name to write the fasta formatted file to.
    :return: Returns True if successfully wrote to file.
    """
    fileName, fileExtension = os.path.splitext(fn)
#    try:
#        assert fileExtension.lower() in [".fasta", ".fa", ".fas", ".fna", ".ffn", ".faa", ".frn"]
#    except AssertionError:
#        _, _, tb = sys.exc_info()
#        traceback.print_tb(tb) # Fixed format
#        tb_info = traceback.extract_tb(tb)
#        filename, line, func, text = tb_info[-1]
#        print(('An error occurred on line {} in statement {}'.format(line, text)))
#        exit(1)
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
#    try:
#        assert fileExtension.lower() in [".fasta", ".fa", ".fas", ".fna", ".ffn", ".faa", ".frn"]
#    except AssertionError:
#        _, _, tb = sys.exc_info()
#        traceback.print_tb(tb) # Fixed format
#        tb_info = traceback.extract_tb(tb)
#        filename, line, func, text = tb_info[-1]
#        print(('An error occurred on line {} in statement {}'.format(line, text)))
#        exit(1)
    dct = {}
    # for sequence in SeqIO.parse(open(fn), "fasta"):
    #     new_key = sequence.description.replace(" ", "_")
    #     if new_key in dct.keys():
    #         print("Duplicate sequence ids found. Exiting")
    #         raise KeyError("Duplicate sequence ids found")
    #     dct[new_key] = str(sequence.seq)
    if sys.version_info[0] < 3:
        my_gen = py2_fasta_iter(fn)
    elif sys.version_info[0] >= 3:
        my_gen = py3_fasta_iter(fn)
    for k, v in my_gen:
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = v
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
