import argparse
from argparse import RawTextHelpFormatter
from smallBixTools import smallBixTools as st
import sys


def main(in_fn, out_fn, side, prcnt):

    # if side == 'fwd':
    #     print("trimming from the front of the alignment.")
    # elif side == 'rev':
    #     print("trimming from the end of the alignment.")
    # elif side == 'both':
    #     print("trimming from both the front and the end.")

    dct = st.fasta_to_dct(in_fn)
    seqs = list(dct.values())
    seq_count = len(seqs)
    #print("we have %s sequences in the input file. " %seq_count)

    #prcnt = 90
    for i in range(len(seqs[0])-1):#columns
        pos_gap_count = 0
        for j in range(seq_count):#rows
            if seqs[j][i] == "-":
                pos_gap_count += 1
        if ((pos_gap_count/seq_count*100.0) < (prcnt)):
            break
    front_trim = i
    #print("counting from the front, we got to col.: %s" %front_trim)

    for i in range(len(seqs[0])-1, 0, -1): #cols, counting down. ie: from back to front.
        pos_gap_count = 0
        for j in range(seq_count): #rows
            if seqs[j][i] == "-":
                pos_gap_count += 1
        if ((pos_gap_count/seq_count*100.0) < (prcnt)):
            break
    rev_trim = i
    #print("counting from the end, we got to col.: %s" %rev_trim)

    dct2 = {}
    for k, v in dct.items():
        if side == 'fwd':
            dct2[k] = v[front_trim:]
        if side == 'rev':
            dct2[k] = v[:rev_trim+1]
        if side == 'both':
            dct2[k] = v[front_trim:rev_trim+1]

    st.dct_to_fasta(dct2, out_fn)

    #print("ending")
    sys.exit(0)


if __name__ == "__main__":
    descrip = '''
    This tool can be used to trip the ends of alignments down to where the majority of the alignmnet lives.
    GCGCTAGATGATCGCTAGCATCGT
    --------TGATCGC--------
    --------TGCTCGC--------
    --------TGACCGC--------
    >>
    TGATCGC
    TGATCGC
    TGCTCGC
    TGACCGC
    We do this by selecting where the frequency of gaps is > 90%, from either the front, end or both sides of the
    sequence.
    '''

    parser = argparse.ArgumentParser(description='' + descrip, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-in', '--in_filename', type=str,
                help='The path to the input file to be trimmed', required=True)
    parser.add_argument('-out', '--out_filename', type=str,
                help='The path to the output file which has been trimmed', required=True)
    parser.add_argument('-side', '--side', type=str,
                help='which side of the alignment do you want trimmed?', required=True, choices=['fwd','rev','both'])
    parser.add_argument('-prcnt', '--percent', type=int,
                help='The percent of the combined alignment which needs to be gap in a column for that column to be '
                     'trimmed out.', required=True)
    options = parser.parse_args()

    in_fn = options.in_filename
    out_fn = options.out_filename
    side = options.side
    prcnt = options.percent

    main(in_fn, out_fn, side, prcnt)