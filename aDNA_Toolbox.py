#!/usr/bin/env python3

'''
Author: Gustaw Eriksson
Date: 2019-06-03

aDNA Toolbox is a program including a variation of functions with which allows the user run a limited set of analysis
on BAM/SAM-files. The program is to be used when running a proof-of-concept analysis on ancient DNA when studying the
prevelance of CT-substitution on aDNA sequence end, determining distance shortest distance between 3' and 5' ends of
overlapping sequences, outputting a graphical outline of sequence alignment between read and reference sequence as well
as reconstructing reference sequences of the ancient DNA molecule from the MD-tag and CIGAR-string.

The script requires a BAM-file, the Samtools program installed locally or in an environment and it is runned by using
argparse. For further information on how to run the program using argparse, please write "./aDNA_Toolbox -h" in the
command terminal.

In case of further questions or bugs, please contact me at erikssongustaw@gmail.com
'''
### Import sys to be able to exit script during error
import sys

### Import collection to set up OrderedDict and calculate frequence of distance to closest 3' on reverse strand
import collections

### Importing argparse for command line interface and subprocess to run samtools
import argparse, subprocess

### Importing numpy and matplotlib to create plots
import matplotlib.pyplot as plt
import numpy as np

### ARGPARS BLOCK ###

### Args.parse to input/output files as well as showing info to user
usage ='Run the program by adding flags and arguments for the desired function. \
For information about each function, see information about each flag.'

parser = argparse.ArgumentParser(description=usage)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 1.0'
    )
parser.add_argument(
    '-b', '--b',
    help='User BAM input file',
    metavar='[BAM_FILE]',
    dest='bam_file',
    )
parser.add_argument(
    '-rr', '--rr',
    help='Output reconstructed reference seq. If output file is not supplied, output is printed in terminal window',
    dest='output_reconstruct_ref',
    action='store_true'
    )
parser.add_argument(
    '-c', '--c',
    help='Output consensus sequence of reads. If output file is not supplied, output is printed in terminal window',
    dest='output_consensus_seq',
    action='store_true'
    )
parser.add_argument(
    '-ct', '--ct',
    help='Output plot illustrating C to T substitution frequeny close to read ends. Outputs plot in user window',
    dest='call_CT_substitution',
    action='store_true'
    )
parser.add_argument(
    '-ag', '--ag',
    help='Output plot illustrating A to G substitution frequeny close to read ends. Outputs plot in user window',
    dest='call_AG_substitution',
    action='store_true'
    )
parser.add_argument(
    '-freq', '--freq',
    help='Output plot illustrating nucleotide frequeny close to read ends. Outputs plot in user window',
    dest='call_nt_freq',
    action='store_true'
    )
parser.add_argument(
    '-frw', '--frw',
    help='Allow forward reads to construct consensus sequence',
    dest='forward_reads_allow',
    action='store_true'
    )
parser.add_argument(
    '-acl', '--acl',
    help='Allow all consensus sequence lengths, i.e allow not fully covered consensus sequences',
    dest='all_cons_lengths',
    action='store_true'
    )
parser.add_argument(
    '-rf', '--rf',
    help='Allow reads with insertions and deletions in analysis',
    dest='remove_filter',
    action='store_true'
    )
parser.add_argument(
    '-over', '--over',
    help="Output figure of 3'-5' distance of overlapping sequences. Outputs plot in user window",
    dest='output_figure',
    action='store_true'
    )
parser.add_argument(
    '-r', '--r',
    help="X-axis range of output 3'-5' distance of overlapping sequences histogram. Default is -100-100",
    metavar='[X-AXIS_RANGE]',
    dest='range_overlap',
    )
parser.add_argument(
    '-o', '--o',
    help='Output file',
    metavar='[OUTPUT_FILE]',
    dest='output_file',
    )

args = parser.parse_args()

### FUNCTIONS BLOCK #####

### Function to seperate letters and numbers
def seperate_letters_numbers(cigar_md):
    ## String variables to store digits and letters
    digit = ''
    alpha = ''

    ## List to store seperated digits and letters char
    cigar_md_list = []
    ## Excludes 'MD:Z:'
    if cigar_md[:5] == 'MD:Z:':
        cigar_md = cigar_md[5:]

    ## Looping over each character in the MD-tag string
    for char in cigar_md:

        ## Seperate and store digits in digit as string
        if char.isdigit() == True:
            digit += char

            ## '0' means zero matches and is excluded
            if digit[0] == '0':
                digit = ''

        ## Seperate and store letters '^' in alpha
        elif char.isalpha() == True or char == '^':
            alpha += char
            alpha = str(alpha)

            ## Letters fall after digits, therefore append to MD-list
            ## is done after letter seperation
            if digit != '':

                ## Digit is turned to integer for future analysis
                digit = int(digit)

                ## Appending to MD-list and cleaning variable
                cigar_md_list.append(digit)
                digit = ''

                cigar_md_list.append(alpha)
                alpha = ''

            ## If no digit is found before the letter
            else:
                ## Appending to MD-list and cleaning variable
                cigar_md_list.append(alpha)
                alpha = ''

    ## If the last character of the md is a digit, it is turned to a integer
    ## and appended to the MD-list
    if digit != '':
        digit = int(digit)
        cigar_md_list.append(digit)

    ## Function returns the MD-list
    return cigar_md_list

### Function to iterate over list of seperated MD-tag and read sequence to return a reconstructed
### reference sequence
def reference_seq(sep_md, sep_cigar, read_seq):

    ## Setting up variables. ref_seq is a empty string for the created sequence, while
    ## start- and end_seq_position are position variables when creating the sequence.
    ref_seq = ''
    start_seq_position = 0
    end_seq_position = 0
    item_index = 0
    ref_seq_extension_count = 0
    ## Variables to handle insertions.
    insertion_seq_position_count = 0
    insertion_seq_position = []
    insertion_seq_count = []
    insertion = False
    insertions_present = 0
    insertions_inserted = 0
    ## Deletion is a flag used when deletion tag '^' is detected in list loop
    deletion = False
    deletion_seq_count = []
    deletion_inserted = 0

    #If an insertion is present and found in the CIGAR string:
    if 'I' in sep_cigar:
        ## Make insertion flag True
        insertion = True
        ## Loop over each item in CIGAR list
        for item in sep_cigar:
            ## If digit is found
            if isinstance(item, int) == True:
                ## The digit is assigned to the digit variable
                digit = item
                ## The sequence position is defined by the sum of digit values
                insertion_seq_position_count += item

            ## If a letter is found and it is a insertion 'I'
            elif item == 'I':
                ## Count the number of insertions found
                insertions_present += 1

                ## The qurrent insertion_seq_count is appended to a list if insertion if found
                insertion_seq_position.append(insertion_seq_position_count)
                ## The number of insertion, i.e. digit before 'I' is appended to a list
                insertion_seq_count.append(digit)

                for count in insertion_seq_count:
                    ref_seq_extension_count += count


    #If an deletion is present and found in the CIGAR string:
    if 'D' in sep_cigar:
        ## Make deletion flag True

        ## Loop over each item in CIGAR list
        for item in sep_cigar:
            ## If digit is found
            if isinstance(item, int) == True:
                ## The digit is assigned to the digit variable
                digit = item

            ## If a letter is found and it is a deletion 'D'
            elif item == 'D':
                ## The number of deletions, i.e. digit before 'D' is appended to a list
                deletion_seq_count.append(digit)


    ## Loop over each item in MD-list
    for item in sep_md:
        item_index += 1

        ## The digits in the list are integers, so true if digit in item. The digits
        ## also imply that the reference and read sequences are matching
        if isinstance(item, int) == True:

            ## The end_seq_position will be the current item
            end_seq_position += item

            if item_index == len(sep_md):
                end_seq_position += ref_seq_extension_count

            ## The ref_seq string will be extended with the matched read sequence
            ## from the start to end position which are are applied on read sequence string
            ref_seq += read_seq[start_seq_position:end_seq_position]

            ## The start position is updated by adding the item
            start_seq_position += item

            ## The deletion flag is turned back to False because it goes from letter to digit
            deletion = False

        ## If the '^' deletion tag is detected, the deletion flag is turned to True
        elif item == '^':
            deletion = True

        ## The Indels in the list are string, so true if letter in item. The letters
        ## imply that the nucleotide will be in the reference sequprint(item_index, sep_md[len(sep_md)-1])ence but not in read sequence
        elif item.isalpha() == True:

            ## If there is a deletion, the flag is true. The nucleotide is in the reference but not
            ## in the read sequence
            if deletion == True:
                ref_seq += item

                ## Check CIGAR if number of deletions match it of the MD-tag. If it does not match. The next
                ## letter will also be treated as if being a deletion. If they do match, the deletion flag
                ## is turned to False and the number of deletion positions is extended by 1.
                if deletion_seq_count[deletion_inserted] == 1:
                    deletion = False
                    deletion_inserted += 1
                    deletion_flag_working = 'OK'

            ## The deletion is False so the nucleotide is a insertion. The nucleotide is in the reference
            ## but does not match the nucleotide found in the read sequence. Start and end sequence
            ## are added with 1.
            elif deletion == False:

                if insertion == True and start_seq_position >= insertion_seq_position[insertions_inserted]:

                    #Extend the start and end seq position with number of insertions
                    end_seq_position += insertion_seq_count[insertions_inserted]
                    ref_seq += read_seq[start_seq_position:end_seq_position]
                    start_seq_position += insertion_seq_count[insertions_inserted]

                    insertions_inserted += 1

                    if insertions_inserted == len(insertion_seq_count):
                        insertion = False

                start_seq_position += 1
                end_seq_position += 1

                ref_seq += item

    return ref_seq

## Function to edit the reference sequence with the CIGAR string
def cigar_modification_ref(sep_cigar, ref_seq):

    ## Setting up variable to keep count of CIGAR position
    seq_position = 0
    mod_ref_seq = ref_seq

    ## Loop over each item in the CIGAR list
    for item in sep_cigar:
        ## The digits in the list are integers, so true if digit in item. The digits
        ## are found both before matches, insertions and deletions. Therefore they are
        ## stored in a variable till letter is identified.
        if isinstance(item, int) == True:
            digit = item

        ## The if statement checks if the item is a letter.
        elif item.isalpha() == True:

            ## If the item is either M or D, the sequence position is extended with the digit
            ## found before the M or D letter
            if item == 'M' or item == 'D':
                seq_position += digit

            ## If the item is I i.e. insertion
            elif item == 'I':

                #To be able to edit the mod_ref_seq, it is turned to a list
                mod_ref_seq = list(mod_ref_seq)

                #If the digit before the I
                if digit == 1:
                    ## The current position in the reference seq, which is where the insertion is,
                    ## will be changed from current nucleotide to '*' to show insertion
                    mod_ref_seq[seq_position] = '*'
                    ## The edited list is returned to being a string
                    mod_ref_seq = ''.join(mod_ref_seq)
                    ## The seq position if extended with 1
                    seq_position += digit

                elif digit > 1:

                    ## end_insertion_position is used because if there is >1 insertions, a range of
                    ## nucleotides have to be changed from nucleotide to '*'
                    end_insertion_position = seq_position + digit
                    ## The reference sequence positions that will be changed are set in the variable
                    ## target position which is a list of position from seq to end insertion position
                    target_position = list(range(seq_position, end_insertion_position))

                    ## Looping over target position list
                    for seq_position in target_position:
                        ## For each position in target position list, the same position in the mod ref seq list
                        ## will be changed from nucleotide to '*' to mark the insertion
                        mod_ref_seq[seq_position] = '*'

                    ## The edited list is returned to being a string
                    mod_ref_seq = ''.join(mod_ref_seq)
                    ## The seq position if extended with the number of insertions
                    seq_position += digit

    return mod_ref_seq

def cigar_modification_read(sep_cigar, read_seq):
    ## Setting up variable to keep count of CIGAR position
    seq_position = 0
    mod_read_seq = read_seq

    ## Loop over each item in the CIGAR list
    for item in sep_cigar:
        ## The digits in the list are integers, so true if digit in item. The digits
        ## are found both before matches, insertions and deletions. Therefore they are
        ## stored in a variable till letter is identified.
        if isinstance(item, int) == True:
            digit = item

        ## Check is item is alphabetic
        elif item.isalpha() == True:

            ## Check CIGAR tag
            if item == 'M' or item == 'I':
                seq_position += digit

            elif item == 'D':

                mod_read_seq = list(mod_read_seq)

                if digit == 1:

                    mod_read_seq.insert(seq_position, '*')
                    mod_read_seq = ''.join(mod_read_seq)
                    seq_position += digit

                elif digit > 1:

                    insertion_list = ['*'] * digit

                    mod_read_seq[seq_position:seq_position] = insertion_list

                    mod_read_seq = ''.join(mod_read_seq)
                    seq_position += digit

    return mod_read_seq

## Function to output graphical alignment between reference and read sequence after
## MD-tag and CIGAR-string modification
def graphic_match_seq(mod_ref_seq, mod_read_seq):

    graphic_match_seq = ''

    mod_read_seq = list(mod_read_seq)
    mod_ref_seq = list(mod_ref_seq)

    ## Parsing over sequences and matching nucleotides
    for ref, read in zip(mod_ref_seq, mod_read_seq):

        if ref == read:
            graphic_match_seq += '.'

        elif ref == '*':
            graphic_match_seq += read

        elif read == '*':
            graphic_match_seq += '*'

        elif ref != read:
            graphic_match_seq += read

    return graphic_match_seq

## Function to calculate the shortest distance between 3' and 5' ends of overlapping sequences
def calculate_distance_closest_3(samfile, range_overlap):

    ### Set length of samfile list
    len_samfile = len(samfile)

    ### Creating list to store reads of same chromosome
    same_chr_reads = []

    ### Creating two list to store forward and reverse reads
    frw_reads = []
    rev_reads = []

    ### List to store the distance to closest 3' on the reverse strand
    distance_closest_3 = []
    #########number_of_0 = 0

    ### Setting flags and variables to seperate reads and control first instances
    first_chr = True
    first_read = False
    seperate_by_flag = False
    new_chr_read = ''
    n_same_chr_reads = 0
    n_f_read = 0
    len_same_chr_reads = 0
    len_frw_reads = 0
    stop_script = False
    first_start_pos_frw_read = 0
    first_frw_read = True
    first_end_pos_rev_read = 0
    first_rev_read = True

    ### Iterating over the list
    for read in samfile:

    ### When last relevant chromosome has been, flag is turned to False and
    ### script is terminated. Script only handles autosomal, sex and mt chromosome
        if stop_script == False:

        ### Setting list item to string and spliting item in list after '\\tt'.
            read = str(read)
            read = read.split('\\t')

        ### Assigning the chromosome to a variable called chr
            chr = read[2]

            ## Seperate reads according to chromosome
            if first_read == True or first_chr == True or chr == current_chr:

                ## Append new_chr_read to list, i.e. the first read with the new chr
                if first_read == True:
                    same_chr_reads.append(new_chr_read)
                    new_chr_read = ''
                    ## Turning first_read flag to False
                    first_read = False

                ## Setting current chr which is seperated
                current_chr = chr

                ## Appending the read to list which stores reads of same chr
                ## Filter out empty items
                same_chr_reads.append(read)

                ## Turning first_chr flag to False
                first_chr = False


            ## When the next chromosome is reached in the bam, the former chromosome is seperated by flag
            elif chr != current_chr:

                ## Save the read in a variable which is then appended to the same_chr_reads list when currenct chr changes
                new_chr_read = read

                ## Turning first_read flag to True
                first_read = True

                ## Change current_chr to new chr
                current_chr = chr

                ## Set length of same_chr_reads list
                len_same_chr_reads = len(same_chr_reads)

                ## When chr has changed, start looping over the
                for chr_read in same_chr_reads:

                    ## Filter out empty items
                    if chr_read != '':

                        if len(current_chr) > 2:
                            stop_script = True

                        ## Keep count on number of reads
                        n_same_chr_reads += 1

                        ## Assigning the start position and seq to seperate variables
                        flag = int(chr_read[1])
                        start_pos_seq = int(chr_read[3])
                        seq = chr_read[9]

                        ### Length of the sequence is determined
                        length_seq = len(seq)

                        ### Depending on if the read is forward or reverse read, which is seen by the flag (read[1]),
                        ### the end position of the sequence will be (-) the start position sequence in reverse reads.
                        ### The end position of the sequence is later used to determine the distance to closest 3'
                        ### on reverse strand.

                        ## If it is the last read of the same chr read list, then next step is to parse thorugh the seperate frw and rev lists
                        if n_same_chr_reads == len_same_chr_reads:
                            ## Restore n_same_chr_reads and len_same_chr_reads to 0
                            n_same_chr_reads = 0
                            len_same_chr_reads = 0
                            ## Last read also have to be seperated
                            if flag is 0:

                                ## Append forward reads to seperate list
                                frw_reads.append(chr_read)

                            elif flag is 16:

                                ## Seperate reverse from forward reads by appending to reverse list
                                rev_reads.append(chr_read)

                                ### Loop over every forward read in the forward list and match it to a reverse read which differs
                                ### the least in regard to distance to 3' on the reverse read. Furthermore, save the distance to
                                ### to later output this in an histogram.

                            ## Set length of frw_reads list
                            len_frw_reads = len(frw_reads)

                            for f_read in frw_reads:
                                ## Variables to count distance to 3' on reverse strand
                                current_distance = None
                                shortest_distance = None
                                shortest_abs_distance = None

                                first_match = True

                                ## Count number of f_reads
                                n_f_read += 1
                                ## Store the start position and chromosome in variables
                                f_chr = f_read[2]
                                start_pos_seq = int(f_read[3])

                                ## For every f_read we loop over the reverse list
                                for r_read in rev_reads:

                                    ## Store the chromosome in variables
                                    r_chr = r_read[2]

                                    ## Store the end position in a variable
                                    end_pos_seq = int(r_read[3])
                                    ## The distance is calculated bu substracting end position with start position.
                                    #current_distance = int(start_pos_seq - end_pos_seq)
                                    current_distance = int(end_pos_seq - start_pos_seq)

                                    ## Abs is used to turn eventuall negative numbers to positive.
                                    current_abs_distance = abs(current_distance)

                                    if current_distance <= -range_overlap:
                                        rev_reads.remove(r_read)

                                    elif current_distance >= range_overlap:
                                        break

                                    else:
                                        ## The first match flag lets the first read be the shortest_distance until shorter is found
                                        if first_match == True or current_abs_distance <= shortest_abs_distance:
                                            ## First match flag is turned to false
                                            first_match = False

                                            ## If the current read distance is shorter than current shortest, then it will be assigned
                                            ## as the current shortest.
                                            shortest_distance = current_distance

                                            ## Turn shortest distance to positive numbers so it can be used in the if-statement
                                            shortest_abs_distance = abs(shortest_distance)

                                if n_f_read == len_frw_reads:
                                    print('CHROMOSOME', f_chr)
                                    ## Restore n_same_chr_reads and len_same_chr_reads to 0:
                                    n_f_read = 0
                                    len_frw_reads = 0

                                    ## Restore first_frw_read flag to True
                                    first_frw_read = True

                                    ## Restore the same_chr, frw and rev list to empty states. This is done to decrease memory and increase speed
                                    same_chr_reads.clear()
                                    frw_reads.clear()
                                    rev_reads.clear()

                                ## Append the shortest distance to 3' on reverse strand to the list of shortest distance to 3'
                                if shortest_distance != None:
                                    distance_closest_3.append(shortest_distance)

                        ## First check if flag is 0 (frw read) or 16 (rev read) to filter out unmapped reads (Flag = 4)
                        elif flag is 0:
                            frw_reads.append(chr_read)

                        elif flag is 16:
                            ## Seperate reverse from forward reads by appending to reverse list
                            rev_reads.append(chr_read)

    return(distance_closest_3)

## Function to call IUPAC nucleotide codes and construct consensus sequence between
## two overlapping reads
def consensus_seq_2_reads(read_nt, overlap_nt):

    if read_nt == overlap_nt:
        consensus_nt = read_nt

    elif overlap_nt == ' ':
        consensus_nt = ' '

    elif read_nt == '*' or overlap_nt == '*':

        if read_nt == '*':
            consensus_nt = overlap_nt

        elif overlap_nt == '*':
            consensus_nt = read_nt

        else:
            consensus_nt = '*'

    elif read_nt != overlap_nt:
        non_matching_nt = read_nt + overlap_nt

        if non_matching_nt == 'AG' or non_matching_nt == 'GA':

            consensus_nt = 'R'

        elif non_matching_nt == 'CT' or non_matching_nt == 'TC':

            consensus_nt = 'Y'

        elif non_matching_nt == 'GC' or non_matching_nt == 'CG':

            consensus_nt = 'S'

        elif non_matching_nt == 'AT' or non_matching_nt == 'TA':

            consensus_nt = 'W'

        elif non_matching_nt == 'GT' or non_matching_nt == 'TG':

            consensus_nt = 'K'

        elif non_matching_nt == 'AC' or non_matching_nt == 'CA':

            consensus_nt = 'M'

        elif 'N' in non_matching_nt:

            consensus_nt = 'N'

    return(consensus_nt)

## Function to call IUPAC nucleotide codes and construct consensus sequence between
## more than two overlapping reads by allowing more IUPAC codes. Uses the function
## above which creates a consensus sequence of the two first overlapping sequences.
def consensus_seq_many_reads(consensus_nt, overlap_nt, first_read_seq, nt_position):

    if consensus_nt == overlap_nt:

        new_consensus_nt = consensus_nt

    elif overlap_nt == ' ' and consensus_nt == ' ':

        new_consensus_nt = first_read_seq[nt_position-1]

    elif overlap_nt == ' ':

        new_consensus_nt = consensus_nt

    elif consensus_nt == ' ':

        new_consensus_nt = consensus_seq_2_reads(first_read_seq[nt_position-1], overlap_nt)

    elif consensus_nt != overlap_nt:

        non_matching_nt = consensus_nt + overlap_nt

        if non_matching_nt == 'AG' or non_matching_nt == 'GA':

            new_consensus_nt = 'R'

        elif non_matching_nt == 'CT' or non_matching_nt == 'TC':

            new_consensus_nt = 'Y'

        elif non_matching_nt == 'GC' or non_matching_nt == 'CG':

            new_consensus_nt = 'S'

        elif non_matching_nt == 'AT' or non_matching_nt == 'TA':

            new_consensus_nt = 'W'

        elif non_matching_nt == 'GT' or non_matching_nt == 'TG':

            new_consensus_nt = 'S'

        elif non_matching_nt == 'AC' or non_matching_nt == 'CA':

            new_consensus_nt = 'M'

        elif non_matching_nt == 'ST' or non_matching_nt == 'YG' or non_matching_nt == 'KC':

            new_consensus_nt = 'B'

        elif non_matching_nt == 'RT' or non_matching_nt == 'WG' or non_matching_nt == 'KA':

            new_consensus_nt = 'D'

        elif non_matching_nt == 'MT' or non_matching_nt == 'WC' or non_matching_nt == 'YA':

            new_consensus_nt = 'H'

        elif non_matching_nt == 'MG' or non_matching_nt == 'RC' or non_matching_nt == 'SA':

            new_consensus_nt = 'V'

        elif 'B' in non_matching_nt or 'D' in non_matching_nt or 'H' in non_matching_nt or 'V' in non_matching_nt:

            new_consensus_nt = 'N'

        else:

            new_consensus_nt = consensus_nt

    return(new_consensus_nt)

## Function to output a consensus sequence from overlapping sequences.
def output_consencus_seq(overlapping_reads, len_overlap_reads):

    ## Setting variables and flags
    first_overlap = True
    first_consensus = True
    n_overlap_reads = 0
    edited_seq = None
    nr_overlap_read = 0

    ## Creating lists for modified overlapping reads and consensus list
    overlapping_seq = []
    overlapping_seq_edited = []
    consensus_seq_list = []

    ## Creating OrderedDict to output read ID and sequence as dictionary
    reads_consensus_dict = collections.OrderedDict()

    for overlap_read in overlapping_reads:

        ## Count number of loops of overlapping reads
        n_overlap_reads += 1

        read = overlap_read

    ### Assining variable to each needed element of the item
        cigar = read[5]
        read_seq = read[9]

    ### Run CIGAR through a function which seperates letters and digits and return a list
        sep_cigar = seperate_letters_numbers(cigar)

    ### Modify the read sequence with CIGAR, i.e. add '-' for deletion position in CIGAR
        mod_read_seq = cigar_modification_read(sep_cigar, read_seq)

    ### Append the modified read seq to list of overlapping seq
        overlapping_seq.append(mod_read_seq)

        ## When the loop has finished, start looping over overlapping_seq and read list
        ## and restore overlapping read count variables
        if len_overlap_reads == n_overlap_reads:

            ## Restoring counting variables
            ##len_overlap_reads = 0
            n_overlap_reads = 0

            ## Looping over overlapping_seq and overlapping_reads lists
            for read, seq in zip(overlapping_reads, overlapping_seq):

                ## If it is the first overlap read/seq, then it is the current target read/seq
                if first_overlap == True:

                    ## Turn flag to false
                    first_overlap = False

                    ## Save the read and seq to two seperate variables to be used later
                    first_read_head = read
                    first_read_seq = seq

                    ## Assing start and end position of the target read
                    first_read_start = int(first_read_head[3])
                    first_read_end = first_read_start + len(first_read_seq)

                    ## Appending first read head and seq to reads_consensus_dict
                    reads_consensus_dict['TARGET SEQ'] = first_read_seq

                ## If it is not the first overlap read/seq, i.e. read which overlaps the target read
                elif first_overlap == False:

                    ## Counting number of overlaps to add in reads_consensus_dict
                    n_overlap_reads += 1

                    ## Calculating difference between start/end position of target and overlap read
                    overlap_diff_start = (first_read_start - int(read[3]))
                    overlap_diff_end = (first_read_end - (int(read[3]) + len(seq)))

                    ## Editing the sequence string using overlap diff variables

                    ## Checking overlap diff at start position,
                    ## If overlap diff at start is >0, then edit overlapping seq
                    if overlap_diff_start > 0:

                        ## Cut out seq before overlap positon
                        edited_seq = seq[overlap_diff_start:]

                    ## If overlap diff at start is <0 then edit the target seq.
                    ## The editing happens in later stage when all other seq have been edited,
                    ## therefore the overlap is appended to a list.
                    elif overlap_diff_start < 0:

                        edited_seq = (' ' * abs(overlap_diff_start)) + seq

                    ## Checking overlap diff at end position
                    ## If overlap diff at the end is >0 then edit the target seq.
                    ## The editing happens in later stage when all other seq have been edited,
                    ## therefore the overlap is appended to a list.
                    if overlap_diff_end > 0:

                        ###overlap_diff_end_first_read.append(overlap_diff_end)
                        if edited_seq == None:

                            edited_seq = seq + (' ' * overlap_diff_end)

                        elif edited_seq != None:

                            edited_seq = edited_seq + (' ' * overlap_diff_end)

                    ## If overlap diff at end is <0, then edit overlapping seq
                    elif overlap_diff_end < 0:

                        ## If the seq has not been edited at its start
                        if edited_seq == None:

                            ## Cut out seq after overlap position
                            edited_seq = seq[:(len(seq) + overlap_diff_end)]


                        ## If the seq has been edited at its start
                        elif edited_seq != None:

                            ## Cut out seq after overlap position
                            edited_seq = edited_seq[:len(edited_seq) + overlap_diff_end]


                    if overlap_diff_start == 0 and overlap_diff_end == 0:

                        edited_seq = seq

                    ## Appending edited sequence to second list of edited sequences

                    overlapping_seq_edited.append(edited_seq)

                    ## Appending seq to reads_consensus_dict
                    overlap_nr = str(n_overlap_reads)
                    overlap_id = 'OVERLAP SEQ' + overlap_nr
                    reads_consensus_dict[overlap_id] = edited_seq

                    ## After appending the edited seq to list of edited seq, reset the edited seq variable
                    edited_seq = None

            for overlap_seq in overlapping_seq_edited:

                if first_consensus == True:

                    first_consensus = False

                    for read_nt, overlap_nt in zip(first_read_seq, overlap_seq):

                        consensus_nt = consensus_seq_2_reads(read_nt, overlap_nt)

                        consensus_seq_list.append(consensus_nt)

                elif first_consensus == False:

                    ## Nucleotide position counter set to 0
                    nt_position = 0

                    for consensus_nt, overlap_nt in zip(consensus_seq, overlap_seq):

                        ## Count nucleotide position in seq
                        nt_position += 1

                        new_consensus_nt = consensus_seq_many_reads(consensus_nt, overlap_nt, first_read_seq, nt_position)

                        consensus_seq_list.append(new_consensus_nt)

                consensus_seq = ''.join(consensus_seq_list)
                consensus_seq_list *= 0

            ## Add consensus seq to reads_consensus_dict
            reads_consensus_dict['CONSENSUS SEQ'] = consensus_seq

    return(reads_consensus_dict)

## Determining shortest distance between overlapping reads to output in file
def set_shortest_3_distance(overlapping_reads):

    first_read = True
    first_distance = True

    for read in overlapping_reads:

        if first_read == True:

            five_end = int(read[3])

            first_read = False

        elif first_read == False:

            three_end = int(read[3])

            distance = three_end - five_end

            if first_distance == True:

                first_distance = False

                shortest_distance = distance

            elif first_distance == False:

                abs_shortest_distance = abs(shortest_distance)
                abs_distance = abs(distance)

                if abs_distance < abs_shortest_distance:

                    shortest_distance = distance

                else:
                    continue

    return(shortest_distance)

### SCRIPT BLOCK ###

if args.bam_file == None or '.bam' not in args.bam_file:
    print('PLEASE ADD BAM FILE')
    sys.exit()

### Subprocess to for BAM to SAM conversion
samfile = subprocess.check_output(['samtools', 'view', args.bam_file])

### Splitlines splits the samfile string along /n creating a list in which each item is a read
samfile = samfile.splitlines()

if args.output_figure:

    if args.range_overlap == None:
        range_overlap = 100

    elif args.range_overlap != None:
        ## Set range for the histogram and distance measure
        range_overlap = int(args.range_overlap)

    distance_closest_3 = calculate_distance_closest_3(samfile, range_overlap)

    ## Order list after increasing distance
    distance_closest_3.sort()

    counter=collections.Counter(distance_closest_3)

    ## Plotting block
    plt.bar(counter.keys(), counter.values(), color='red')
    plt.title("Number of sequences at different distances to closest 3' on reverse strand")
    plt.xlabel("Distance to closest 3' on reverse strand")
    plt.ylabel('Number of sequences')
    plt.grid(True)
    plt.show()

if args.output_consensus_seq or args.output_reconstruct_ref or args.call_CT_substitution or args.call_nt_freq:

    stop_script = False

    ### Creating list to store reads of same chromosome
    same_chr_reads = []
    ### Creating two list to store forward and reverse reads
    frw_reads = []
    rev_reads = []

    ## Create list to store overlapping reads
    overlapping_reads = []

    ## Create list to store distance from ends to CT substitution
    CT_substitution_distance_list_5_end = []
    CT_substitution_distance_list_3_end = []

    ## Create list to store distance from ends to CT substitution
    GA_substitution_distance_list_5_end = []
    GA_substitution_distance_list_3_end = []

    ## Create list to store distance from ends to CT substitution
    Other_substitution_distance_list_5_end = []
    Other_substitution_distance_list_3_end = []

    ## Create list to store frequency of nt close to read ends
    A_freq_list_5_end = []
    A_freq_list_3_end = []
    T_freq_list_5_end = []
    T_freq_list_3_end = []
    C_freq_list_5_end = []
    C_freq_list_3_end = []
    G_freq_list_5_end = []
    G_freq_list_3_end = []
    Other_freq_list_5_end = []
    Other_freq_list_3_end = []


    ## Setting up variables and flags
    first_chr = True
    first_read = False
    new_chr_read = ''
    n_same_chr_reads = 0
    n_f_read = 0
    len_same_chr_reads = 0
    num_bins = 0
    X_axis_range = 10
    nt_sub_reads = 0

    ### Iterating over the list
    for read in samfile:

    ### When last relevant chromosome has been, flag is turned to False and
    ### script is terminated. Script only handles autosomal, sex and mt chromosome
        if stop_script == False:

            ## Saving original form of read ID before splitting
            orig_header = read

        ### Setting list item to string and spliting item in list after '\\tt'.
            read = str(read)

            read = read.split('\\t')

        ## MD-tag found in both column 21, 22 and 24. Within the columns, there are str that is not MD-tag.
        ## Additional filtering is required to extract MD-tag
            cigar = read[5]
            pre_md = read[21:]

        ## Parse column 21- and extract the MD-tag column based on 'MD:Z:' string
            for item in pre_md:
                if 'MD:Z:' in item:
                    md = item

        ## Remove reads contaning soft clips, hard clips or padding
            if 'S' in cigar or 'H' in cigar or 'N' in cigar or 'P' in md:
                continue

            if args.output_reconstruct_ref:

            ### Assining variable to each needed element of the item
                full_id = read[0:]
                id = read[0]
                flag = read[1]
                read_seq = read[9]

            ### Run MD through a function which seperates letters and digits and return a list
                sep_md = seperate_letters_numbers(md)

            ### Run CIGAR through a function which seperates letters and digits and return a list
                sep_cigar = seperate_letters_numbers(cigar)

            ### Iterate over MD-list containing seperated letters and digits to return a reconstructed reference sequence
                ref_seq = reference_seq(sep_md, sep_cigar, read_seq)

            ### Modify the reconstructed reference sequence with CIGAR, i.e. add '-' for insertion position in CIGAR
                mod_ref_seq = cigar_modification_ref(sep_cigar, ref_seq) #Maybe add IF statement to check if read even has a 'I' in it to increase speed

            ### Modify the read sequence with CIGAR, i.e. add '-' for deletion position in CIGAR
                mod_read_seq = cigar_modification_read(sep_cigar, read_seq)

            ### Return a string which illustrates the relationship between read and reference sequence. Used in output
                output_graphic_seq = graphic_match_seq(mod_ref_seq, mod_read_seq)

                if args.output_file:

                    with open(args.output_file, 'a') as fout:

                        print('{}\n{}\n{}\n{}\n'.format(id, mod_ref_seq, output_graphic_seq, mod_read_seq), file=fout)

                else:
                    print('{}\n{}\n{}\n{}\n'.format(id, mod_ref_seq, output_graphic_seq, mod_read_seq))

            elif args.output_consensus_seq or args.call_CT_substitution or args.call_AG_substitution or args.call_nt_freq:

                ## Default filters out deletions and insertions in reads
                if '^' in md or 'I' in cigar or 'D' in cigar:

                    ## Allow insertions and deletions if remove filter flag is added
                    if args.remove_filter:
                        pass

                    else:
                        continue

            ### Assigning the chromosome to a variable called chr
                chr = read[2]

                ## Seperate reads according to chromosome
                if first_read == True or first_chr == True or chr == current_chr:

                    ## Append new_chr_read to list, i.e. the first read with the new chr
                    if first_read == True:
                        same_chr_reads.append(new_chr_read)
                        new_chr_read = ''
                        ## Turning first_read flag to False
                        first_read = False

                    ## Setting current chr which is seperated
                    current_chr = chr

                    ## Appending the read to list which stores reads of same chr
                    ## Filter out empty items
                    same_chr_reads.append(read)

                    ## Turning first_chr flag to False
                    first_chr = False


                ## When the next chromosome is reached in the bam, the former chromosome is seperated by flag
                elif chr != current_chr:

                    ## Save the read in a variable which is then appended to the same_chr_reads list when currenct chr changes
                    new_chr_read = read

                    ## Turning first_read flag to True
                    first_read = True

                    ## Change current_chr to new chr
                    current_chr = chr

                    ## Set length of same_chr_reads list
                    len_same_chr_reads = len(same_chr_reads)

                    ## When chr has changed, start looping over the list of same chromosome
                    for chr_read in same_chr_reads:

                        ## Filter out empty items
                        if chr_read != '':

                            ## Stopping the script from running over other lines than reads
                            if len(current_chr) > 2:
                                stop_script = True

                            ## Keep count on number of reads
                            n_same_chr_reads += 1

                            ## Assigning the start position ########and seq to seperate variables
                            flag = int(chr_read[1])
                            seq = chr_read[9]

                            ## Length of the sequence is determined
                            seq = chr_read[9]
                            length_seq = len(seq)

                            ### Depending on if the read is forward or reverse read, which is seen by the flag (read[1]),
                            ### the end position of the sequence will be (-) the start position sequence in reverse reads.
                            ### The end position of the sequence is later used to determine the distance to closest 3'
                            ### on reverse strand.

                            ## If it is the last read of the same chr read list, then next step is to parse thorugh the seperate frw and rev lists
                            if n_same_chr_reads == len_same_chr_reads:

                                ## Restore n_same_chr_reads and len_same_chr_reads to 0
                                n_same_chr_reads = 0
                                len_same_chr_reads = 0
                                ## Last read also have to be seperated
                                if flag is 0:
                                    ## Append forward reads to seperate list
                                    frw_reads.append(chr_read)
                                elif flag is 16:

                                    ## Seperate reverse from forward reads by appending to reverse list
                                    rev_reads.append(chr_read)
                                    ### Loop over every forward read in the forward list and match it to a reverse read which differs
                                    ### the least in regard to distance to 3' on the reverse read. Furthermore, save the distance to
                                    ### to later output this in an histogram.

                                ## Set length of frw_reads list
                                len_frw_reads = len(frw_reads)

                                for f_read in frw_reads:

                                    n_f_read += 1

                                    ## Store forward read key values in variables
                                    f_read_id = f_read[0]
                                    f_chr = f_read[2]
                                    f_start_pos_seq = int(f_read[3])
                                    f_seq = f_read[9]
                                    f_seq_len = len(f_read[9])
                                    f_end_pos_seq = f_start_pos_seq + f_seq_len

                                    ## Appending end position to the end of the read header
                                    f_read.append(f_end_pos_seq)

                                    ## Appending f_read to list of overlapping_reads
                                    overlapping_reads.append(f_read)

                                    ## Each forward read is thereafter looped against the list again, to
                                    ## match it to any overlapping forward reads.

                                    if args.forward_reads_allow:

                                        for f_2_read in frw_reads:
                                            ## Store the start position,end position and chromosome in variables
                                            f_2_chr = f_2_read[2]
                                            f_2_start_pos_seq = int(f_2_read[3])
                                            f_2_seq_len = len(f_2_read[9])
                                            f_2_end_pos_seq = f_2_start_pos_seq + f_2_seq_len

                                            if f_2_read == f_read:
                                                continue

                                            ## If the end position of the second frw read is before the start position
                                            ## of the current (first frw read), there is no overlap. There will neither
                                            ## be any overlap with future frw reads because the start position is increasing
                                            ## when looping over the frw read list. Therefore it is deleted to save memory.
                                            elif f_2_end_pos_seq < f_start_pos_seq:
                                                frw_reads.remove(f_2_read)

                                            ## If the end position of the second frw read is after the start position
                                            ## of the current read (first frw read) and the start position of the second frw
                                            ## read is before the end position of the current read, there is overlap between
                                            ##the reads. The read is therefore appended to the list of overlapping_reads.
                                            elif f_2_end_pos_seq >= f_start_pos_seq and f_2_start_pos_seq <= f_end_pos_seq:
                                                overlapping_reads.append(f_2_read)

                                            ## If the start position of the second frw read falls after the end position of
                                            ## the current read (f_read), the loop is broken because there will be no more
                                            ## overlaps for any of the reads with the current read.
                                            elif f_2_start_pos_seq > f_end_pos_seq:
                                                break

                                    ## Each forward read is looped against a reverse read to determine if there is any
                                    ## overlap. In case of overlap, the read is appended to the overlapping_reads list.
                                    for r_read in rev_reads:

                                        ## Store the start position,end position and chromosome in variables
                                        r_chr = r_read[2]
                                        r_end_pos_seq = int(r_read[3])
                                        r_seq_len = len(r_read[9])
                                        r_start_pos_seq = r_end_pos_seq + r_seq_len

                                        ## If the start position of the reverse read is before the forward read start position,
                                        ## there can be no overlap and the reverse read is deleted from the rev_reads list to
                                        ## save memory.
                                        if r_start_pos_seq < f_start_pos_seq:
                                            rev_reads.remove(r_read)

                                        ## If the end position of the second frw read is after the start position
                                        ## of the current read (first frw read) and the start position of the second frw
                                        ## read is before the end position of the current read, there is overlap between
                                        ##the reads. The read is therefore appended to the list of overlapping_reads.
                                        elif r_start_pos_seq >= f_start_pos_seq and r_end_pos_seq <= f_end_pos_seq:
                                            overlapping_reads.append(r_read)

                                        ## If the start position of the second frw read falls after the end position of
                                        ## the current read (f_read), the loop is broken because there will be no more
                                        ## overlaps for any of the reads with the current read.
                                        elif r_start_pos_seq > f_end_pos_seq:
                                            break

                                    if len(overlapping_reads) > 1:
                                        ## Set length of overlapping reads list
                                        len_overlap_reads = len(overlapping_reads)

                                        consensus_seq = output_consencus_seq(overlapping_reads, len_overlap_reads)

                                        if args.all_cons_lengths:

                                            with open(args.output_file, 'a') as fout:

                                                for id, seq in consensus_seq.items():
                                                    print('{}\t{}'.format(id, seq), file=fout)

                                        elif ' ' not in (consensus_seq['CONSENSUS SEQ']):

                                            if args.output_consensus_seq:

                                                if args.output_file:

                                                    with open(args.output_file, 'a') as fout:

                                                        print(f_read_id, file=fout)
                                                        for id, seq in consensus_seq.items():
                                                            print('{}\t{}'.format(id, seq), file=fout)
                                                        print('\n', file=fout)

                                                else:
                                                    print(f_read_id)
                                                    for id, seq in consensus_seq.items():
                                                        print('{}\t{}'.format(id, seq))
                                                    print('\n')

                                            elif args.call_CT_substitution or args.call_nt_freq:

                                                ## Filter out sequences below a length of 20 as the functions analysing read ends require
                                                ## minimum length of 10.
                                                if len(f_seq) < X_axis_range * 2:
                                                    continue

                                            ### Assining variable to each needed element of the item
                                                f_cigar = f_read[5]
                                                f_pre_md = f_read[21:]

                                                ## Exctract MD-tag
                                                for item in f_pre_md:
                                                    item = str(item)
                                                    if 'MD:Z:' in item:
                                                        f_md = item
                                                        break

                                            ### Run MD through a function which seperates letters and digits and return a list
                                                f_sep_md = seperate_letters_numbers(f_md)

                                            ### Run CIGAR through a function which seperates letters and digits and return a list
                                                f_sep_cigar = seperate_letters_numbers(f_cigar)

                                            ### Iterate over MD-list containing seperated letters and digits to return a reconstructed reference sequence
                                                f_ref_seq = reference_seq(f_sep_md, f_sep_cigar, f_seq)

                                            ### Modify the reconstructed reference sequence with CIGAR, i.e. add '-' for insertion position in CIGAR
                                                f_mod_ref_seq = cigar_modification_ref(f_sep_cigar, f_ref_seq)

                                            ### Determine shortest 5' to 3' on reverse strand distance
                                                shortest_3_distance = set_shortest_3_distance(overlapping_reads)

                                            ## Setting up counting variables
                                                nt_sub_pos = 0
                                                nt_sub_reads += 1
                                                nt_freq_amount = 0

                                                ## Iterating over forward reference seq and consensus seq nt
                                                for ref_nt, cons_nt in zip(f_mod_ref_seq, consensus_seq['TARGET SEQ']):

                                                    nt_sub_pos += 1
                                                    nt_sub_pos_neg = abs(len(f_mod_ref_seq)-nt_sub_pos) + 1

                                                    ## Only use 10 first and last nucleotides of the read
                                                    if X_axis_range >= nt_sub_pos or (len(f_mod_ref_seq)-X_axis_range+1) <= nt_sub_pos:

                                                        ## Check nucleotide frequency at read ends
                                                        if args.call_nt_freq:

                                                            if nt_sub_pos > len(f_mod_ref_seq)/2:

                                                                if cons_nt == 'A':
                                                                    A_freq_list_3_end.append(nt_sub_pos_neg)

                                                                elif cons_nt == 'T':
                                                                    T_freq_list_3_end.append(nt_sub_pos_neg)

                                                                elif cons_nt == 'C':
                                                                    C_freq_list_3_end.append(nt_sub_pos_neg)

                                                                elif cons_nt == 'G':
                                                                    G_freq_list_3_end.append(nt_sub_pos_neg)
                                                                else:
                                                                    Other_freq_list_3_end.append(nt_sub_pos_neg)

                                                            else:

                                                                if cons_nt == 'A':
                                                                    A_freq_list_5_end.append(nt_sub_pos)

                                                                elif cons_nt == 'T':
                                                                    T_freq_list_5_end.append(nt_sub_pos)

                                                                elif cons_nt == 'C':
                                                                    C_freq_list_5_end.append(nt_sub_pos)

                                                                elif cons_nt == 'G':
                                                                    G_freq_list_5_end.append(nt_sub_pos)
                                                                else:
                                                                    Other_freq_list_5_end.append(nt_sub_pos)

                                                        ## If reference do not match consensus seq, check for substitutions of C to T, A to G or other
                                                        if ref_nt != cons_nt:

                                                            if args.call_CT_substitution and ref_nt == 'C' and cons_nt == 'T' or cons_nt == 'Y':

                                                                f_read_sub_pos = f_start_pos_seq + nt_sub_pos

                                                                if nt_sub_pos > len(f_mod_ref_seq)/2:


                                                                    CT_substitution_distance_list_3_end.append(nt_sub_pos_neg)

                                                                else:

                                                                    CT_substitution_distance_list_5_end.append(nt_sub_pos)

                                                            elif args.call_AG_substitution and ref_nt == 'G' and cons_nt == 'A' or cons_nt == 'R':


                                                                f_read_sub_pos = f_start_pos_seq + nt_sub_pos

                                                                if nt_sub_pos > len(f_mod_ref_seq)/2:

                                                                    GA_substitution_distance_list_3_end.append(nt_sub_pos_neg)

                                                                else:

                                                                    GA_substitution_distance_list_5_end.append(nt_sub_pos)

                                                            else:

                                                                f_read_sub_pos = f_start_pos_seq + nt_sub_pos

                                                                if nt_sub_pos > len(f_mod_ref_seq)/2:

                                                                    Other_substitution_distance_list_3_end.append(nt_sub_pos_neg)

                                                                else:

                                                                    Other_substitution_distance_list_5_end.append(nt_sub_pos)

                                                            if args.output_file != None:
                                                                with open(args.output_file, 'a') as fout:

                                                                    print('{}\t{}\t{}\t{}\t{}'.format(f_read_id, f_chr, f_read_sub_pos, nt_sub_pos, shortest_3_distance), file=fout)

                                                        else:
                                                            continue

                                                    else:
                                                        continue

                                    overlapping_reads.clear()

                                    if n_f_read == len_frw_reads:
                                        print('CHROMOSOME', f_chr)
                                        ## Restore n_same_chr_reads and len_same_chr_reads to 0:
                                        n_f_read = 0
                                        len_frw_reads = 0

                                        ## Restore first_frw_read flag to True
                                        first_frw_read = True

                                        ## Restore the same_chr, frw and rev list to empty states. This is done to decrease memory and increase speed
                                        same_chr_reads.clear()
                                        frw_reads.clear()
                                        rev_reads.clear()

                            ## First check if flag is 0 (frw read) or 16 (rev read) to filter out unmapped reads (Flag = 4)
                            elif flag is 0:
                                frw_reads.append(chr_read)

                            elif flag is 16:
                                ## Seperate reverse from forward reads by appending to reverse list
                                rev_reads.append(chr_read)

## Outputting plots block ##
if args.call_CT_substitution and args.call_nt_freq:

    X = list(range(1, X_axis_range+1))

    ## Setting up X and Y variables
    CT_substitution_distance_list_5_end.sort()
    CT_substitution_distance_list_3_end.sort()

    counter_CT_5 = collections.Counter(CT_substitution_distance_list_5_end)
    counter_CT_3 = collections.Counter(CT_substitution_distance_list_3_end)

    counter_CT_5 = {nt:nr/nt_sub_reads for nt, nr in counter_CT_5.items()}
    counter_CT_3 = {nt:nr/nt_sub_reads for nt, nr in counter_CT_3.items()}

    sub_lists_CT_5 = sorted(counter_CT_5.items()) # sorted by key, return a list of tuples
    x_CT_5, y_CT_5 = zip(*sub_lists_CT_5) # unpack a list of pairs into two tuples

    sub_lists_CT_3 = sorted(counter_CT_3.items())
    x_CT_3, y_CT_3 = zip(*sub_lists_CT_3)

    if args.call_AG_substitution:

        GA_substitution_distance_list_5_end.sort()
        GA_substitution_distance_list_3_end.sort()

        counter_GA_5 = collections.Counter(GA_substitution_distance_list_5_end)
        counter_GA_3 = collections.Counter(GA_substitution_distance_list_3_end)

        counter_GA_5 = {nt:nr/nt_sub_reads for nt, nr in counter_GA_5.items()}
        counter_GA_3 = {nt:nr/nt_sub_reads for nt, nr in counter_GA_3.items()}

        sub_lists_GA_5 = sorted(counter_GA_5.items())
        x_GA_5, y_GA_5 = zip(*sub_lists_GA_5)

        sub_lists_GA_3 = sorted(counter_GA_3.items())
        x_GA_3, y_GA_3 = zip(*sub_lists_GA_3)

    Other_substitution_distance_list_5_end.sort()
    Other_substitution_distance_list_3_end.sort()

    counter_Other_5 = collections.Counter(Other_substitution_distance_list_5_end)
    counter_Other_3 = collections.Counter(Other_substitution_distance_list_3_end)

    counter_Other_5 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_5.items()}
    counter_Other_3 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_3.items()}

    sub_lists_Other_5 = sorted(counter_Other_5.items())
    x_Other_5, y_Other_5 = zip(*sub_lists_Other_5)

    sub_lists_Other_3 = sorted(counter_Other_3.items())
    x_Other_3, y_Other_3 = zip(*sub_lists_Other_3)

    A_freq_list_5_end.sort()
    A_freq_list_3_end.sort()
    T_freq_list_5_end.sort()
    T_freq_list_3_end.sort()
    C_freq_list_5_end.sort()
    C_freq_list_3_end.sort()
    G_freq_list_5_end.sort()
    G_freq_list_3_end.sort()
    Other_freq_list_5_end.sort()
    Other_freq_list_3_end.sort()

    counter_A_5 = collections.Counter(A_freq_list_5_end)
    counter_A_3 = collections.Counter(A_freq_list_3_end)
    counter_T_5 = collections.Counter(T_freq_list_5_end)
    counter_T_3 = collections.Counter(T_freq_list_3_end)
    counter_C_5 = collections.Counter(C_freq_list_5_end)
    counter_C_3 = collections.Counter(C_freq_list_3_end)
    counter_G_5 = collections.Counter(G_freq_list_5_end)
    counter_G_3 = collections.Counter(G_freq_list_3_end)
    counter_Other_5 = collections.Counter(Other_freq_list_5_end)
    counter_Other_3 = collections.Counter(Other_freq_list_3_end)

    counter_A_5 = {nt:nr/nt_sub_reads for nt, nr in counter_A_5.items()}
    counter_A_3 = {nt:nr/nt_sub_reads for nt, nr in counter_A_3.items()}
    counter_T_5 = {nt:nr/nt_sub_reads for nt, nr in counter_T_5.items()}
    counter_T_3 = {nt:nr/nt_sub_reads for nt, nr in counter_T_3.items()}
    counter_C_5 = {nt:nr/nt_sub_reads for nt, nr in counter_C_5.items()}
    counter_C_3 = {nt:nr/nt_sub_reads for nt, nr in counter_C_3.items()}
    counter_G_5 = {nt:nr/nt_sub_reads for nt, nr in counter_G_5.items()}
    counter_G_3 = {nt:nr/nt_sub_reads for nt, nr in counter_G_3.items()}
    counter_Other_5 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_5.items()}
    counter_Other_3 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_3.items()}

    sub_lists_A_5 = sorted(counter_A_5.items()) # sorted by key, return a list of tuples
    x_A_5, y_A_5 = zip(*sub_lists_A_5) # unpack a list of pairs into two tuples
    sub_lists_A_3 = sorted(counter_A_3.items())
    x_A_3, y_A_3 = zip(*sub_lists_A_3)

    sub_lists_T_5 = sorted(counter_T_5.items())
    x_T_5, y_T_5 = zip(*sub_lists_T_5)
    sub_lists_T_3 = sorted(counter_T_3.items())
    x_T_3, y_T_3 = zip(*sub_lists_T_3)

    sub_lists_C_5 = sorted(counter_C_5.items())
    x_C_5, y_C_5 = zip(*sub_lists_C_5)
    sub_lists_C_3 = sorted(counter_C_3.items())
    x_C_3, y_C_3 = zip(*sub_lists_C_3)

    sub_lists_G_5 = sorted(counter_G_5.items())
    x_G_5, y_G_5 = zip(*sub_lists_G_5)
    sub_lists_G_3 = sorted(counter_G_3.items())
    x_G_3, y_G_3 = zip(*sub_lists_G_3)

    sub_lists_Other_5 = sorted(counter_Other_5.items())
    x_Other_5, y_Other_5 = zip(*sub_lists_Other_5)
    sub_lists_Other_3 = sorted(counter_Other_3.items())
    x_Other_3, y_Other_3 = zip(*sub_lists_Other_3)

    ## Plotting block ##

    plt.subplot(2,2,1)
    plt.plot(X, y_CT_5, 'r', label='C -> T')
    if args.call_AG_substitution:
        plt.plot(X, y_GA_5, 'g', label='G -> A')
    plt.plot(X, y_Other_5, 'b', label='Other')
    plt.title("Number of C to T substitutions close to 5' and 3' ends")
    plt.xlabel("Distance to 5' end")
    plt.ylabel('Frequency of substitutions')
    plt.legend()
    plt.grid(True)

    plt.subplot(2,2,3)
    plt.plot(X, y_CT_3, 'r', label='C -> T')
    if args.call_AG_substitution:
        plt.plot(X, y_GA_3, 'g', label='G -> A')
    plt.plot(X, y_Other_3, 'b', label='Other')
    plt.xlabel("Distance to 3' end")
    plt.ylabel('Frequency of substitutions')
    plt.legend()
    plt.grid(True)

    plt.subplot(2,2,2)
    plt.plot(X, y_A_5, 'r', label='A')
    plt.plot(X, y_T_5, 'b', label='T')
    plt.plot(X, y_C_5, 'g', label='C')
    plt.plot(X, y_G_5, 'y', label='G')
    #plt.plot(X, y_Other_5, 'black', label='Undetermined')
    plt.title("Nucleotide frequency close to 5' and 3' ends")
    plt.xlabel("Distance to 5' end")
    plt.ylabel('Nucleotide Frequency')
    plt.legend()
    plt.grid(True)

    plt.subplot(2,2,4)
    plt.plot(X, y_A_3, 'r', label='A')
    plt.plot(X, y_T_3, 'b', label='T')
    plt.plot(X, y_C_3, 'g', label='C')
    plt.plot(X, y_G_3, 'y', label='G')
    #plt.plot(X, y_Other_3, 'black', label='Undetermined')
    plt.xlabel("Distance to 3' end")
    plt.ylabel('Nucleotide Frequency')
    plt.legend()
    plt.grid(True)

    plt.show()

elif args.call_CT_substitution or args.call_AG_substitution:

    X = list(range(1, X_axis_range+1))

    ## Setting up X and Y variables
    CT_substitution_distance_list_5_end.sort()
    CT_substitution_distance_list_3_end.sort()

    counter_CT_5 = collections.Counter(CT_substitution_distance_list_5_end)
    counter_CT_3 = collections.Counter(CT_substitution_distance_list_3_end)

    counter_CT_5 = {nt:nr/nt_sub_reads for nt, nr in counter_CT_5.items()}
    counter_CT_3 = {nt:nr/nt_sub_reads for nt, nr in counter_CT_3.items()}

    sub_lists_CT_5 = sorted(counter_CT_5.items()) # sorted by key, return a list of tuples
    x_CT_5, y_CT_5 = zip(*sub_lists_CT_5) # unpack a list of pairs into two tuples

    sub_lists_CT_3 = sorted(counter_CT_3.items())
    x_CT_3, y_CT_3 = zip(*sub_lists_CT_3)

    if args.call_AG_substitution:

        GA_substitution_distance_list_5_end.sort()
        GA_substitution_distance_list_3_end.sort()

        counter_GA_5 = collections.Counter(GA_substitution_distance_list_5_end)
        counter_GA_3 = collections.Counter(GA_substitution_distance_list_3_end)

        counter_GA_5 = {nt:nr/nt_sub_reads for nt, nr in counter_GA_5.items()}
        counter_GA_3 = {nt:nr/nt_sub_reads for nt, nr in counter_GA_3.items()}

        sub_lists_GA_5 = sorted(counter_GA_5.items())
        x_GA_5, y_GA_5 = zip(*sub_lists_GA_5)

        sub_lists_GA_3 = sorted(counter_GA_3.items())
        x_GA_3, y_GA_3 = zip(*sub_lists_GA_3)

    Other_substitution_distance_list_5_end.sort()
    Other_substitution_distance_list_3_end.sort()

    counter_Other_5 = collections.Counter(Other_substitution_distance_list_5_end)
    counter_Other_3 = collections.Counter(Other_substitution_distance_list_3_end)

    counter_Other_5 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_5.items()}
    counter_Other_3 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_3.items()}

    sub_lists_Other_5 = sorted(counter_Other_5.items())
    x_Other_5, y_Other_5 = zip(*sub_lists_Other_5)

    sub_lists_Other_3 = sorted(counter_Other_3.items())
    x_Other_3, y_Other_3 = zip(*sub_lists_Other_3)

    ## Plotting block ##

    plt.subplot(2,1,1)
    plt.plot(X, y_CT_5, 'r', label='C -> T')
    if args.call_AG_substitution:
        plt.plot(X, y_GA_5, 'g', label='G -> A')
    plt.plot(X, y_Other_5, 'b', label='Other')
    plt.title("Number of C to T substitutions close to 5' and 3' ends")
    plt.xlabel("Distance to 5' end")
    plt.ylabel('Frequency of substitutions')
    plt.legend()
    plt.grid(True)

    plt.subplot(2,1,2)
    plt.plot(X, y_CT_3, 'r', label='C -> T')
    if args.call_AG_substitution:
        plt.plot(X, y_GA_3, 'g', label='G -> A')
    plt.plot(X, y_Other_3, 'b', label='Other')
    plt.xlabel("Distance to 3' end")
    plt.ylabel('Frequency of substitutions')
    plt.legend()
    plt.grid(True)

    plt.show()

elif args.call_nt_freq:

    X = list(range(1, X_axis_range+1))

    A_freq_list_5_end.sort()
    A_freq_list_3_end.sort()
    T_freq_list_5_end.sort()
    T_freq_list_3_end.sort()
    C_freq_list_5_end.sort()
    C_freq_list_3_end.sort()
    G_freq_list_5_end.sort()
    G_freq_list_3_end.sort()
    Other_freq_list_5_end.sort()
    Other_freq_list_3_end.sort()

    counter_A_5 = collections.Counter(A_freq_list_5_end)
    counter_A_3 = collections.Counter(A_freq_list_3_end)
    counter_T_5 = collections.Counter(T_freq_list_5_end)
    counter_T_3 = collections.Counter(T_freq_list_3_end)
    counter_C_5 = collections.Counter(C_freq_list_5_end)
    counter_C_3 = collections.Counter(C_freq_list_3_end)
    counter_G_5 = collections.Counter(G_freq_list_5_end)
    counter_G_3 = collections.Counter(G_freq_list_3_end)
    counter_Other_5 = collections.Counter(Other_freq_list_5_end)
    counter_Other_3 = collections.Counter(Other_freq_list_3_end)

    counter_A_5 = {nt:nr/nt_sub_reads for nt, nr in counter_A_5.items()}
    counter_A_3 = {nt:nr/nt_sub_reads for nt, nr in counter_A_3.items()}
    counter_T_5 = {nt:nr/nt_sub_reads for nt, nr in counter_T_5.items()}
    counter_T_3 = {nt:nr/nt_sub_reads for nt, nr in counter_T_3.items()}
    counter_C_5 = {nt:nr/nt_sub_reads for nt, nr in counter_C_5.items()}
    counter_C_3 = {nt:nr/nt_sub_reads for nt, nr in counter_C_3.items()}
    counter_G_5 = {nt:nr/nt_sub_reads for nt, nr in counter_G_5.items()}
    counter_G_3 = {nt:nr/nt_sub_reads for nt, nr in counter_G_3.items()}
    counter_Other_5 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_5.items()}
    counter_Other_3 = {nt:nr/nt_sub_reads for nt, nr in counter_Other_3.items()}

    sub_lists_A_5 = sorted(counter_A_5.items()) # sorted by key, return a list of tuples
    x_A_5, y_A_5 = zip(*sub_lists_A_5) # unpack a list of pairs into two tuples
    sub_lists_A_3 = sorted(counter_A_3.items())
    x_A_3, y_A_3 = zip(*sub_lists_A_3)

    sub_lists_T_5 = sorted(counter_T_5.items())
    x_T_5, y_T_5 = zip(*sub_lists_T_5)
    sub_lists_T_3 = sorted(counter_T_3.items())
    x_T_3, y_T_3 = zip(*sub_lists_T_3)

    sub_lists_C_5 = sorted(counter_C_5.items())
    x_C_5, y_C_5 = zip(*sub_lists_C_5)
    sub_lists_C_3 = sorted(counter_C_3.items())
    x_C_3, y_C_3 = zip(*sub_lists_C_3)

    sub_lists_G_5 = sorted(counter_G_5.items())
    x_G_5, y_G_5 = zip(*sub_lists_G_5)
    sub_lists_G_3 = sorted(counter_G_3.items())
    x_G_3, y_G_3 = zip(*sub_lists_G_3)

    sub_lists_Other_5 = sorted(counter_Other_5.items())
    x_Other_5, y_Other_5 = zip(*sub_lists_Other_5)
    sub_lists_Other_3 = sorted(counter_Other_3.items())
    x_Other_3, y_Other_3 = zip(*sub_lists_Other_3)

    ## Plotting block ##

    plt.subplot(2,1,1)
    plt.plot(X, y_A_5, 'r', label='A')
    plt.plot(X, y_T_5, 'b', label='T')
    plt.plot(X, y_C_5, 'g', label='C')
    plt.plot(X, y_G_5, 'y', label='G')
    #plt.plot(X, y_Other_5, 'black', label='Undetermined')
    plt.title("Nucleotide frequency close to 5' and 3' ends")
    plt.xlabel("Distance to 5' end")
    plt.ylabel('Nucleotide Frequency')
    plt.legend()
    plt.grid(True)

    plt.subplot(2,1,2)
    plt.plot(X, y_A_3, 'r', label='A')
    plt.plot(X, y_T_3, 'b', label='T')
    plt.plot(X, y_C_3, 'g', label='C')
    plt.plot(X, y_G_3, 'y', label='G')
    #plt.plot(X, y_Other_3, 'black', label='Undetermined')
    plt.xlabel("Distance to 3' end")
    plt.ylabel('Nucleotide Frequency')
    plt.legend()
    plt.grid(True)

    plt.show()
