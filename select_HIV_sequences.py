import argparse
from argparse import RawTextHelpFormatter
import sys
import os
from smallBixTools import smallBixTools as st


def main(in_fn, region, ft, out_good_fn, out_bad_fn):
    # example command /uct/dev/code/small_bix_tools$ python3.4 select_HIV_sequences.py -in /uct/dev/source/smallBixTools/swipe_hiv_region_selection/2016_01_27_MurrayLogan_Murray-Logan_SFF_nomatch.fasta -region gag -out_matching /uct/dev/source/smallBixTools/swipe_hiv_region_selection/out_matching -out_nonmatching /uct/dev/source/smallBixTools/swipe_hiv_region_selection/out_nonmatching

    # software location
    # /home/dave/Software/swipe-2.0.5/swipe

    # database: made from makeblastdb. we have made nucleotide databases
    # -d /uct/ref/own_made_blastdb/hxb2_gag/hxb2_gag_blastdb

    # is the databse a protein / nucleotide database. In our case we are using nucleotide - so we use 0 here.
    # -p 0

    # input file. We are using a fasta file for input here
    # -i /uct/dev/source/smallBixTools/swipe_hiv_region_selection/2016_01_27_MurrayLogan_Murray-Logan_SFF_nomatch.fasta

    # number of threads to use. running on a machine which reports to have 8 cores, I run with 7 threads,
    # and it does not max any of them. In fact, each only gets to about 40%. Running 4 instances of swipe,
    # simultaniously on the same machine with 8 cores, each instance being told it can use 7 threads, all the cores
    # get to about 70% each. -- Makes me want to try specifying many more threads (32) on an 8 core machine.
    # -a 7

    # output format [0,7-9=plain,xml,tsv,tsv+]. valid optinos are: 0, 7, 8, 9.
    # -m 7

    # strand. we want to search forward, and reverse. It takes longer, but we can be more sure.
    # -S 3

    # the output file
    # -o /uct/dev/source/smallBixTools/swipe_hiv_region_selection/swipe_outfile

    print("main")
    print("input file: %s" %in_fn)
    print("gene region: %s" %region)
    print("ft: %s" %ft)
    print("output matching file: %s" %out_good_fn)
    print("output non-matching file: %s" %out_bad_fn)
    print("\n\n\n")

    # first part:
    # TODO

    # swipe outfile: /uct/dev/source/smallBixTools/swipe_hiv_region_selection/swipe_outfile_m8
    #swipe_fn = "/uct/dev/source/smallBixTools/swipe_hiv_region_selection/swipe_outfile_m8"
    swipe_fn = "/uct/dev/source/smallBixTools/swipe_hiv_region_selection/swipe_outfile_m8_adjusted"
    swipe_fn = "/uct/Breakthrough/analysis/HVTN503/2016-07-07_data_transfer/29/demultiplex/swipe_pol_out_m8"
    # second part.
    # splitting a fasta file (based on swipe outfile) into "good" vs. "bad".

    dct = st.fasta_to_dct(in_fn)
    dct2 = {}
    for k, v in dct.items():
        dct2[k] = {'seq': v,
                   'len': len(v),}

    print(len(dct))

    # parsing swipe output
    swipe_data = []
    with open(swipe_fn, "r") as fh:
        for line in fh:
            swipe_data.append(line.strip().split("\t"))
    print(len(swipe_data))
    print(swipe_data[0])

    longest_aln_swipe_dct = {}
    for row in swipe_data:
        seq_id = row[0]
        aln_len = row[3]
        if seq_id in longest_aln_swipe_dct.keys():
            if aln_len > longest_aln_swipe_dct[seq_id][3]:
                longest_aln_swipe_dct[seq_id] = row
        else:
            longest_aln_swipe_dct[seq_id] = row

    print("we have selected only the alingments from the swipe results, where if there were two alignments, "
          "the longest is represented here.")
    print("len of this: %s" %(len(longest_aln_swipe_dct)))
    cleaned_swipe_data = []


    # get all the sequences from the fasta file where the alignment length is longer than 80% of read length
    longer_alignments = []
    shorter_alignments = []
    prcnt = 60
    for k, v in longest_aln_swipe_dct.items():
        seq_id = v[0]
        ident = v[2]
        aln_len = v[3]
        try:
            aln_len = int(aln_len)
        except Exception as e:
            print(e)
            sys.exit()

        seq_len = dct2[seq_id]['len']

        if (aln_len / seq_len) < prcnt/100.0:
            shorter_alignments.append(v)
        else:
            longer_alignments.append(v)

    with open("short_alignments.fasta", "w") as fw:
        for row in shorter_alignments:
            fw.write(">" + row[0] + "\n" + dct[row[0]] + "\n")
    with open("longer_alignments.fasta", "w") as fw:
        for row in longer_alignments:
            fw.write(">" + row[0] + "\n" + dct[row[0]] + "\n")

    print("of all the longest alignments, the number which were below %s percent alignmed are: %s"
          %(prcnt, len(shorter_alignments)))
    print("the number which were above %s percent aligned are: %s" %(prcnt, len(longer_alignments)))

    print("ending")
    sys.exit(0)
# python3.4 select_HIV_sequences.py -in /uct/dev/source/smallBixTools/swipe_hiv_region_selection/2016_01_27_MurrayLogan_Murray-Logan_SFF_nomatch.fasta -region gag -out_matching /uct/dev/source/smallBixTools/swipe_hiv_region_selection/out_matching -out_nonmatching /uct/dev/source/smallBixTools/swipe_hiv_region_selection/out_nonmatching
# vsearch --cluster_fast short_alignments.fasta --id 0.8 --centroids centroids.fas


if __name__ == "__main__":
    descrip = '''
    Input file (fastq / fasta)
    Input region of HIV (gene = gag, pol, env, nef)
    Input file type (if not specified, we will try to get this from the input file. Only .fasta and .fastq allowed)
    Output matching filename
    Output not matching filename

    We take either a fasta or fastq file, and use the software swipe to pairwise align each of its sequences against
    a own made blastdb of HIV, and only sequences of high enough identity, which cover a long enough region will
    be written to the "matching" file, everything else gets written to the "nonmatching" file.


    '''

    parser = argparse.ArgumentParser(description='' + descrip, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-in', '--in_filename', type=str,
                help='The path to the input file to have HIV sequences selected from it.', required=True)
    parser.add_argument('-region', '--region', type=str, help='The gene region which you wish to search against.',
                        choices=['gag', 'pol', 'env', 'nef'], required=True)
    parser.add_argument('-ft', '--filetype', type=str, help='If the infile does not have a file extension, you can '
                                                            'specify the file type here.',
                        required=False, choices=['fasta', 'fastq'])
    parser.add_argument('-out_matching', '--out_matching', type=str,
                help='The path to the output file to be populated with matching HIV sequences from the infile.',
                        required=True)
    parser.add_argument('-out_nonmatching', '--out_nonmatching', type=str,
                help='The path to the output file to be populated with matching HIV sequences from the infile.',
                        required=True)
    options = parser.parse_args()

    in_fn = options.in_filename
    region = options.region
    ft = options.filetype
    out_good_fn = options.out_matching
    out_bad_fn = options.out_nonmatching

    if not ft:
        if os.path.splitext(in_fn)[1] not in [".fasta", ".fastq"]:
            print("File extension found to be: " + os.path.splitext(in_fn)[1])
            print("We don't recognize the file type extension. Please either use: '.fastq', or '.fasta' only.")
            sys.exit()

    main(in_fn, region, ft, out_good_fn, out_bad_fn)
