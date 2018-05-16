import subprocess
from itertools import groupby
from operator import itemgetter
import sys, traceback, os
from Bio import SeqIO
import operator
import hashlib
from Bio import pairwise2


__version__ = "0.0.20"
__author__ = "David Matten"
__credits__ = ["David Matten", "Colin Anthony", "Carolyn Williamson"]
__license__ = "GPL"
__maintainer__ = "David Matten"
__email__ = "david.matten@uct.ac.za"
__status__ = "Development"


def test_cmd_present(cmd):
    '''
    This tests if cmd is executable. This can be used to check if some
    software is available at the given path and executable on a machine.
    For example, if you want to use mafft, and the user says its located at:
    "/opt/not_really_here/mafft", we an test and return False.
    :param cmd: the command to test.
    :return: Bool. True if command is available at path and executable. False if not.
    '''
    from shutil import which
    return which(cmd) is not None


def find_start(s):
    '''
    Given a string, find the position of the first character which is not a gap "-" hyphen.
    :param s: String to search for the first non gap character in
    :return: position of first non-gap character
    '''
    import re
    s = s.upper()
    start = re.search(r'[^-]', s).start()
    return start


def find_end(s):
    '''
    Given a string, s, find the position of the last character which is not a gap "-" hyphen.
    :param s: String to search for the last non gap character in
    :return: position of last non-gap character
    '''
    import re
    # reverse the string.
    s = s[::-1]
    s = s.upper()
    # the last character which is not a gap is the first one in the reversed string.
    # get its position by taking length minus that number.
    end = len(s) - re.search(r'[^-]', s).start()
    return end


def combine_align_trim(align_to_me, fasta_file_to_be_added, out_file, mafft_call='mafft'):
    '''
    This method adds the second argument file (fasta_file_to_be_added) onto the first. It uses mafft add by default
    to align the second fasta file onto the first.
    :param align_to_me: This is the first source file to be aligned onto. The file which will be trimmed to. Typically
    This is a file with shorter reads.
    :param fasta_file_to_be_added:
    :param out_file: The path to the fasta file where results should be written. This file will be an aligned combined
                        version of the two other input source files.
    :param mafft_call: A non-required argument. How should we call mafft on this machine? Typical mafft installs allow
    for the default mafft call: "mafft"
    :return: No return.
    '''
    if not test_cmd_present(mafft_call):
        print("We require mafft to run this method. Please make sure it is installed, or specify how to run it on this"
              "machine. Get it from: https://mafft.cbrc.jp/alignment/software/")
    if not os.path.isfile(align_to_me):
        print("Specified file is not a file.")
        raise ValueError("Specified file is not a file. {}".format(align_to_me))
    if not os.path.isfile(fasta_file_to_be_added):
        print("Specified file is not a file.")
        raise ValueError("Specified file is not a file. {}".format(fasta_file_to_be_added))
    if os.path.isfile(out_file):
        print("The output file will be overwriten: {}".format(out_file))

    cmd = "mafft --thread -1 --add {} {} > {}".format(align_to_me, fasta_file_to_be_added, out_file)
    try:
        subprocess.call(cmd, shell=True)
    except Exception as e:
        print(e)
        print("Error trying to call mafft to add the OGV sequences onto the existing NGS alignment., or removing or "
              "copying files.")
        raise

    align_to_me_dct = fasta_to_dct(align_to_me)

    starting_positions = [find_start(align_to_me_dct[k]) for k in align_to_me_dct.keys()]
    start_pos = int(min(starting_positions))
    ending_positions = [find_end(align_to_me_dct[k]) for k in align_to_me_dct.keys()]
    end_pos = int(max(ending_positions)) - 1
    print("The starting position found is: {}".format(start_pos))
    print("The ending position found is: {}".format(end_pos))
    print("Trimming to this region.")

    out_dct = fasta_to_dct(out_file)
    for k in out_dct.keys():
        out_dct[k] = out_dct[k][start_pos:end_pos]

    dct_to_fasta(out_dct, out_file)


def find_best_global_between_fastas(target_fn, query_fn, csv_out_fn):
    """
    For every sequence in the .fasta formatted query_fn, find the best matching sequence from the target_fn based on
    a global alignment from the BioPython package.
    :param target_fn: The fasta file to check against. Sequences from the query file will be searched against this file.
    :param query_fn: The file to check from. Each one of these sequences will have its best match in the target file
    searched for.
    :param csv_out_fn: Each sequence from the query file will have a single row in this csv file. Column 1 is the query
     seqid. Column 2 is the best found target seqid. Column 3 is an identity count.
    :return: No return. Writes output to csv_out_fn
    """
#    seq_align_call = "/home/dave/Software/seq-align/bin/needleman_wunsch"
    target_dct = fasta_to_dct(target_fn)
    query_dct = fasta_to_dct(query_fn)

    for query_seqid, query_seq in query_dct.items():
        for target_seqid, target_seq in target_dct.items():
            alignments = pairwise2.align.globalxx(query_seq, target_seq)
            print(len(alignments))
            print(alignments[0])

    # TODO complete. Writing out to csv.


def sanitize_fasta_seqids(infile, outfile, valid_chars):
    """
    Read a fasta formatted file. Remove all characters which are not in the valid_chars string.
    :param infile: input file name to check.
    :param outfile: output file name to write to.
    :param valid_chars: string of valid characters.
                        Typically one might use "{}{}".format(string.ascii_letters, string.digits)
    :return: No return. Writes to the outfile specified.
    """
    dct = fasta_to_dct(infile)
    new_dct = {}
    for k, v in dct.items():
        new_k = ""
        for char in k:
            if char in valid_chars:
                new_k += char
        if len(new_k) > 0:
            new_dct[new_k] = v
        else:
            print("If we remove all the illegal characters, there is nothing left to use as a sequence id.")
            raise Exception
    dct_to_fasta(new_dct, outfile)


def compare_seqs_of_fasta_files(fn1, fn2):
    """
    Compares the sequences of two fasta files. Does not compare seqids. Prints a list of the seqids for the sequences
    which are found in one fasta file, but not the other, for both comparison directions. Also returns this a tuple
    of lists.
    :param fn1: fasta file 1
    :param fn2: fasta file 1
    :return: returns a tuple of lists. ([in fn1 and not in fn2], [in fn2 and not in fn1], [in both])
    """
    dct1 = fasta_to_dct(fn1)
    dct2 = fasta_to_dct(fn2)
    all_fn1_seqs = list(dct1.values())
    all_fn2_seqs = list(dct2.values())
    print(all_fn1_seqs)
    in1_not_in2 = []
    in2_not_in1 = []
    in_both = []
    for seqid1, seq1 in dct1.items():
        if seq1 in all_fn2_seqs:
            in_both.append(seqid1)
        else:
            in1_not_in2.append(seqid1)
    for seqid2, seq2 in dct2.items():
        if seq2 not in all_fn1_seqs:
            in2_not_in1.append(seqid2)

    print("All sequences from file 1, which were not found in file 2, their sequence ids are: ")
    print(in1_not_in2)
    print("All sequences from file 2, which were not found in file 1, their sequence ids are: ")
    print(in2_not_in1)

    print("All sequences found in both file 1 and file 2, their sequence ids are: ")
    print(in_both)

    return in1_not_in2, in2_not_in1, in_both


def make_hash_of_seqids(src_fn, out_fn):
    """
    When calling mafft - sequence ids over 253 in length are truncated. This can result in non-unique ids if the first
    253 characters of the seqid are the same, with a difference following that.
    To get around this - we can has the sequence ids, and write a new .fasta file for mafft to work on, then
    translate the sequence ids back afterwards.

    This function does the hashing and writing to file.

    This is a sibling function to: unmake_hash_of_seqIDS

    Will raise an exception on error

    :param src_fn: the src file to operate on
    :param out_fn: an output file is produced, with the modified sequence ids.
    :return: returns a lookup dictionary for finding the new sequence id name. This dictionary can be used with the
    sibling function: "unmake_hash_of_seqIDS".
    """
    try:
        dct = fasta_to_dct(src_fn)
        new_dct = {}
        translation_dct = {}
        for k, v in dct.items():
            new_k = hashlib.md5(k.encode()).hexdigest()
            new_dct[new_k] = v
            translation_dct[k] = new_k
        dct_to_fasta(new_dct, out_fn)
    except Exception as e:
        print(e)
        raise

    return translation_dct


def unmake_hash_of_seqids(lookup_dict, src_fn, out_fn):
    """
    When calling mafft - sequence ids over 253 in length are truncated. This can result in non-unique ids if the first
    253 characters of the seqid are the same, with a difference following that.
    To get around this - we can has the sequence ids, and write a new .fasta file for mafft to work on, then
    translate the sequence ids back afterwards.

    This function does the translation back afterwards.

    This is a sibling function to: make_hash_of_seqIDS.

    Will raise an exception on error

    :param lookup_dict:
    :param src_fn:
    :param out_fn:
    :return: no return value.
    """
    try:
        dct = fasta_to_dct(src_fn) # the dictionary with hashed seqids.
        back_translated_dct = {} # the 'original' sequence ids are back.
        for old_k, hash_k in lookup_dict.items():
            back_translated_dct[old_k] = dct[hash_k]
        dct_to_fasta(back_translated_dct, out_fn)
    except Exception as e:
        print(e)
        raise


def compare_fasta_files(file1, file2):
    """
    Compares two fasta files, to see if they contain the same data. The sequences must be named the same. We check if
    sequence A from file 1 is the same as sequence A from file 2.
    The order in the files does not matter.
    Gaps are considered.
    :param file1: first fasta file
    :param file2: second fasta file
    :return: True if the files contain the same data. False if the files contain different data.
    """
    dct1 = fasta_to_dct(file1)
    dct2 = fasta_to_dct(file2)
    for seqid, seq1 in dct1.items():
        if seqid not in dct2.keys():
            return False
        seq2 = dct2[seqid]
        if seq1 != seq2:
            return False
    return True


def countNinPrimer(primer_seq):
    """
    Motifbinner2 requires values to be specified for primer id length and primer length. Its tiresome to have to
    calculate this for many strings. So, I wrote this to help myself.
    An example of a primer sequence might be: NNNNNNNAAGGGCCAAAGGAACCCTTTAGAGACTATG
    And we would like to know how many N's there are, how many other characters there are, and what the combined
    total lenght is.
    :param primer_seq: the primer sequence to have calculations performed on.
    :return: nothing. prints to stdout are done.
    """
    n_count = primer_seq.count("N")
    other_count = primer_seq.count("A") + primer_seq.count("C") + primer_seq.count("G") + primer_seq.count("T")
    print("N: {}".format(n_count))
    print("!N: {}".format(other_count))
    print("total: {}".format(n_count + other_count))
    min_score = 0.85 * (n_count + other_count)
    print("suggested min score: {}".format(min_score))


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
    :param fasta_name: The fasta formatted file to parse.
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def py2_fasta_iter(fasta_name):
    """
    from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    :param fasta_name: The fasta formatted file to parse.
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq


def try_call(cmd, logging=None):
    """
    Try a subprocess call on a string cmd. Raises an error on exceptions. Logging is allowed.
    :param cmd: The string to try calling.
    :param logging: Can optionally take a logging arg. This is the default logging for Python.
    :return: Returns the return code from calling the string command.
    """
    rtrn_code = None
    try:
        captured_stdout = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        if logging != None:
            logging.info(captured_stdout)
    except Exception as e:
        if logging != None:
            logging.warning(e)
        print(e)
        raise
    return rtrn_code


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
    seqs = [v for k, v in dct.items()]
    cons = ""
    for i in range(len(seqs[0])):
        pos = {}
        for k in range(len(seqs)):
            letter = seqs[k][i]
            if letter in pos.keys():
                pos[letter] += 1
            else:
                pos[letter] = 1
        if sys.version_info[0] < 3:
            cons += max(pos.iteritems(), key=operator.itemgetter(1))[0]
        elif sys.version_info[0] >= 3:
            cons += max(pos.items(), key=operator.itemgetter(1))[0]
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
    """
    Attempts to automatically remove duplicate sequences from the specifed file. Writes results to output file
    specified. Uses BioPython SeqIO to parse the in file specified. Replaces spaces in the sequence id with underscores.
    Itterates over all sequences found - for each one, checking if its key already exists in an accumulating, if it
    does: check if the sequence which each specifies is the same. If they have the same key, and the same sequence -
    then keep the second instance encountered. Once the file has been parsed - write to the output file specified all
     sequences found which
    Will raise an exception if an error occurs during execution.
    :param in_fn: The file to check. Full path to file is required.
    :param out_fn: Output
    :return: No return value.
    """
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
    """
    A distance measure between two iterables. Typically meant for DNA sequence strings. eg: ATCG and A-CG would be a
    distance of 1. ATTCG to A--CG distance 1 also.
    Gap scores: gap opening: -1. gap extension: penalty of zero. This means any length gap counts the same as a 1 length gap.
    A second disconnected gap counts as an additional -1.
    :param s1: first iterable to compare.
    :param s2: second iterable to compare.
    :return: returns the distance. int.
    """
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
    Checks which version of Python is being used, calls the appropriate iterator.
    Spaces in the sequence ids are replaced with underscores.
    Duplicate sequence ids are not allowed. An error will be raised.
    :param fn: The fasta formatted file to read from.
    :return: A dictionary of the contents of the fasta file specified. The dictionary in the format:
             {sequence_id: sequence_string, sequence_id2: sequence_2, etc.}
    """
    dct = {}
    if sys.version_info[0] < 3:
        my_gen = py2_fasta_iter(fn)
    elif sys.version_info[0] >= 3:
        my_gen = py3_fasta_iter(fn)
    for k, v in my_gen:
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = str(v).replace("~", "_")
    return dct


def hamdist(str1, str2):
    """
    Use this after aligning sequences.
    This counts the number of differences between equal length str1 and str2
    The order of the input sequences does not matter.
    :param str1: The first sequence.
    :param str2: The second sequence.
    :return: Returns an int count of the number of differences between the two input strings, considered per position.
    """
    return sum(el1 != el2 for el1, el2 in zip(str1, str2))


def normalized_hamdist(str1, str2):
    """
    Use this after aligning sequences.
    This counts the number of differences between equal length str1 and str2
    The order of the input sequences does not matter.
    :param str1: The first sequence.
    :param str2: The second sequence.
    :return: Returns a float value of the number of differences divided by the length of the first input argument.
    """
    return sum(el1 != el2 for el1, el2 in zip(str1, str2))/float(len(str1))


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
