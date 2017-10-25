# smallBixTools
a few small functions for bioinformatics



# smallBixTools a few small functions for bioinformatics.

See readme for full details.

Repo location:

https://bitbucket.org/DavidMatten/small_bix_tools

available on VM.

List of functions:

get_regions_from_panel:

Slices regions out of a fasta formatted file, joins them together, and writes the resulting fasta file to the given location.
an example call might be: get_regions_from_panel("test.fasta", 0, 10], [20, 30, "/tmp", "outfile.fasta")
which would, for each sequence in the input file: "test.fasta", take the region from 0 to 10 joined with the
region from 20 to 30, and write the result to the file: "/tmp/outfile.fasta".

find_ranges

Find contiguous ranges in a list of numerical values.
eg: data = [1,2,3,4,8,9,10]
find_ranges(data) will return:
1, 2, 3, 4], [8, 9, 10

hamdist

Use this after aligning sequences.
This counts the number of differences between equal length str1 and str2
The order of the input sequences does not matter.

fasta_to_dct

a dictionary of the contents of the file name given. Dictionary in the format:
{sequence_id: sequence_string, id_2: sequence_2, etc.}

dct_to_fasta

:param d: dictionary in the form: {sequence_id: sequence_string, id_2: sequence_2, etc.}
:param fn: The file name to write the fasta formatted file to.
:return: Returns True if successfully wrote to file.

find_duplicate_ids

customdist

hyphen_to_underscore_fasta

auto_duplicate_removal

Attempts to automatically remove duplicate sequences from the specifed file. Writes results to output file
specified. Uses BioPython SeqIO to parse the in file specified. Replaces spaces in the sequence id with underscores.
Itterates over all sequences found - for each one, checking if its key already exists in an accumulating, if it
does: check if the sequence which each specifies is the same. If they have the same key, and the same sequence -
then keep the second instance encountered. Once the file has been parsed - write to the output file specified all
sequences found which
Will raise an exception if an error occurs during execution.

build_cons_seq

# https://www.biostars.org/p/14026/

own_cons_maker

split_file_into_timepoints

size_selector

py2_fasta_iter

from Brent Pedersen: https://www.biostars.org/p/710/#1412
given a fasta file. yield tuples of header, sequence

py3_fasta_iter

modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
given a fasta file. yield tuples of header, sequence

convert_count_to_frequency_on_fasta

when running vsearch as such:
vsearch –cluster_fast {} –id 0.97 –sizeout –centroids {}
We get a centroids.fasta file with seqid header lines like:
>ATTCCGGTATCT_9;size=1432;
>CATCATCGTAAG_14;size=1;
etc.
This method converts those count values into frequencies.
Notes: The delimiter between sections in the sequence id must be ";".
There must be a section in the sequence id which has exactly: "size=x" where x is an integer.
This must be surrounded by ";"'s

countNinPrimer

Motifbinner2 requires values to be specified for primer id length and primer length. Its tiresome to have to
calculate this for many strings. So, I wrote this to help myself.
An example of a primer sequence might be: NNNNNNNAAGGGCCAAAGGAACCCTTTAGAGACTATG
And we would like to know how many N's there are, how many other characters there are, and what the combined
total lenght is.

compare_fasta_files

Compares two fasta files, to see if they contain the same data. The sequences must be named the same. We check if
sequence A from file 1 is the same as sequence A from file 2.
The order in the files does not matter.
Gaps are considered.

unmake_hash_of_seqids

When calling mafft - sequence ids over 253 in length are truncated. This can result in non-unique ids if the first
253 characters of the seqid are the same, with a difference following that.
To get around this - we can has the sequence ids, and write a new .fasta file for mafft to work on, then
translate the sequence ids back afterwards.

This function does the translation back afterwards.

This is a sibling function to: make_hash_of_seqIDS.

Will raise an exception on error

make_hash_of_seqids

When calling mafft - sequence ids over 253 in length are truncated. This can result in non-unique ids if the first
253 characters of the seqid are the same, with a difference following that.
To get around this - we can has the sequence ids, and write a new .fasta file for mafft to work on, then
translate the sequence ids back afterwards.

This function does the hashing and writing to file.

This is a sibling function to: unmake_hash_of_seqIDS

Will raise an exception on error
