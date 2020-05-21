

             _____  _   _            _______          _ _               
            |  __ \| \ | |   /\     |__   __|        | | |              
        __ _| |  | |  \| |  /  \       | | ___   ___ | | |__   _____  __
       / _` | |  | | . ` | / /\ \      | |/ _ \ / _ \| | '_ \ / _ \ \/ /
      | (_| | |__| | |\  |/ ____ \     | | (_) | (_) | | |_) | (_) >  <
       \__,_|_____/|_| \_/_/    \_\    |_|\___/ \___/|_|_.__/ \___/_/\_\


                A proof-of-concept toolbox for ancient DNA

# README

Author: Gustaw Eriksson

Date: 06-06-2019

Version: 1.0

Contact: erikssongustaw@gmail.com

# Introduction
aDNA Toolbox is a command-line toolbox applied on BAM-files, to be used for
analyzing aDNA to illustrate concepts of C to T substitutions, patterns of
overlapping reads and nucleotide frequency at the ends of reads. Furthermore,
the toolbox allows the user to reconstruct aDNA molecules by reconstructing it
with the use of MD-tags and CIGAR-string. The program was developed using Python
version 3.7.1 and the Atom Editor version 1.29.0.

The program requires the use of samtools, which either can be downloaded locally
from samtools website http://www.htslib.org/download/ or using package manager.
aDNA Toolbox was written with samtools version 1.9 using htslib version 1.9.

aDNA Toolbox also uses several Python modules. The version used were
the following:

sys 3.7                              subprocess 3.7
collections 3.7                      matplot.pyplot 3.1.0
argparse 3.2                         numpy 1.3.0

# Usage:
The program utilizes argparse, allowing the user to call inputting a BAM-file
and call different analytic functions. Some of the functions allows the output
the data in an output file which the user has to supply. If an output is not
supplied, the output will be printed in the terminal. Other functions output
the data in plot format in the user window. For further information on
functions, flag and output format, please see below:

-h, --help                            Show help message

-v, --version                         Show program's version number

-b [BAM_FILE], --b [BAM_FILE]         User BAM input file

-rr, --rr                             Output reconstructed reference seq. If
                                      output file is not supplied, output is
                                      printed in terminal window

-c, --c                               Output consensus sequence of reads. If
                                      output file is not supplied, output is
                                      printed in terminal window

-ct, --ct                             Output plot illustrating C to T
                                      substitution frequeny close to read ends.
                                      Outputs plot in user window

-ag, --ag                             Output plot illustrating A to G
                                      substitution frequeny close to read ends.
                                      Outputs plot in user window

-freq, --freq                         Output plot illustrating nucleotide
                                      frequeny close to read ends. Outputs plot
                                      in user window

-frw, --frw                           Allow forward reads to construct
                                      consensus sequence

-acl, --acl                           Allow all consensus sequence lengths, i.e
                                      allow not fully covered consensus
                                      sequences

-rf, --rf                             Allow reads with insertions and deletions
                                      in analysis

-over, --over                         Output figure of 3'-5' distance of
                                      overlapping sequences. Outputs plot in
                                      user window

-r [X-AXIS_RANGE], --r [X-AXIS_RANGE] X-axis range of output 3'-5' distance of
                                      overlapping sequences histogram. Default
                                      is -100-100

-o [OUTPUT_FILE], --o [OUTPUT_FILE]   User output file

Flags can also be shown in the terminal by writing:

  $ ./aDNA_Toolbox -h

The program is runned by the following command:

  $ ./aDNA_Toolbox -b [BAM FILE] [-flag]

Besides calling the program and supplying a BAM-file, a flag with our without a
output file has to be provided.

# WARNING:
In its current state, the program will filter out reads containing soft (S in
cigar) and hard clips (H in cigar), spliced ends (N in cigar) as well as
padding (P in md) when calling the -ct, -ag and -freq flags.

The flags -frw and -acl which removes filters are not fully functional and
increases memory use, i.e. removal of the filters increase runtime.

The program is still in development and is to be used with caution
