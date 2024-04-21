#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Create a count matrix from a sorted_bam file.

Sort reads by name and count:
- if the 2 paired reads map the same target --> +1
- if the 2 paired reads map different targets --> +1,+1
- if one of the paired reads map a target --> +1

Usage:
    bam_map2count.py [--help] [--version] [--threads THREADS] --threshold <integer> REFERENCE IDENTIFIER REFFAI OUTCOUNTS OUTABUNDANCE OUTLEARN BAM ...

Options:
  --help
  --version
  -t THREADS, --threads THREADS   number of threads used [default: 1]

https://wabi-wiki.scilifelab.se/display/KB/Filter+uniquely+mapped+reads+from+a+BAM+file

XA: to report alternative sites and SA is a new flag to tag so-called split alignments for chimera reads
"""

__author__ = "Matthieu Pichaud"
__version__ = "1.1"
from optparse import OptionParser
import subprocess
import os
import string
import fileinput
from collections import Counter
import sys
from pprint import pprint

##
# Supporting functions
##


def run_command(command):
        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=True)
        return iter(p.stdout.readline, b'')


##
# Main
##

if __name__ == "__main__":
        
        parser = OptionParser()
        parser.add_option("-t", "--input.threads", dest="in_threads_int", help="in ref file", metavar="FILENAME")
        parser.add_option("-y", "--input.threshold", dest="in_threshold_int", help="in ref file", metavar="FILENAME")
        parser.add_option("-i", "--input.bam", dest="in_bam_file", help="in bam file", metavar="FILENAME")
        parser.add_option("-v", "--input.reffai", dest="in_reffai_file", help="in reffai file", metavar="FILENAME")
        #parser.add_option("-r", "--input.ref", dest="in_ref_file", help="in ref file", metavar="FILENAME")
        parser.add_option("-s", "--input.sample", dest="in_sample_name", help="in reffai file", metavar="FILENAME")
        parser.add_option("-c", "--output.counts", dest="out_counts_file", help="in reffai file", metavar="FILENAME")
        parser.add_option("-a", "--output.abundance", dest="out_abundance_file", help="in reffai file", metavar="FILENAME")
        parser.add_option("-l", "--output.learn", dest="out_learn_file", help="in reffai file", metavar="FILENAME")
        
        (options, args) = parser.parse_args()
        
        
        #arguments = docopt(__doc__, version="{} version:{}".format(os.path.basename(__file__), __version__))
        #identifier = arguments["IDENTIFIER"]
        #reference = arguments["REFERENCE"]
        #reffai = arguments["REFFAI"]
        #bam = arguments["BAM"]
        #Nthreads = arguments["--threads"]
        #read_threshold = int(arguments["--threshold"])
        #OUTCOUNTS=arguments["OUTCOUNTS"]
        #OUTABUNDANCE=arguments["OUTABUNDANCE"]
        #OUTLEARN=arguments["OUTLEARN"]

        print('Making count dictionary')
        # Initialize count_d using reference file
        count_d = {}
        with open(options.in_reffai_file, 'r') as fi:  # ADDED instead of file for Python 3.0
                for line in fi:
                        count_d[line.rstrip("\r\n").split("\t")[0]] = 0
        print('Making length dictionary')
        # Make a dictionary with the lengths of each gene for calculating abundance
        gene_lengths_d = {}
        with open(options.in_reffai_file, 'r') as fi:
                for line in fi:
                        gene_lengths_d[line.rstrip("\r\n").split("\t")[0]] = int(line.split('\t')[1])

        print('Starting counting')
        # Get counts
        neighbor_l = []
        i = 0
        for bam_ in [options.in_bam_file]:
                #cmd = "samtools sort -n -@{} {}|samtools view -@{} - | grep -v 'XA:Z' | grep -v 'SA:Z'".format(Nthreads, bam_, Nthreads)
                cmd = "samtools sort -n -@{} {}|samtools view -@{} - ".format(options.in_threads_int, bam_, options.in_threads_int)
                readid_prev = ""
                hit_prev = ""
                for i, line in enumerate(run_command(cmd)):
                        line = line.decode('UTF8')  # ADDED for Python 3.0
                        items = line.rstrip("\r\n").split("\t")
                        if len(items) >= 2:
                                readid = items[0]
                                hit = items[2]
                                if readid_prev == readid:
                                        # increment the read1 hit
                                        count_d[hit] += 1
                                        if hit_prev != hit:
                                                print("Found alternative hits")
                                                print(neighbor_l)
                                                # increment the read2 hi
                                                count_d[hit_prev] += 1
                                                # report neighbor
                                                neighbor_l.append((min(hit_prev, hit), max(hit_prev, hit)))
                                                # reset
                                        readid_prev = ""
                                        hit_prev = ""
                                elif readid_prev != "":
                                        # read is different from previous read, increment the previous read hit
                                        count_d[hit_prev] += 1
                                        readid_prev = readid
                                        hit_prev = hit
                                else:
                                        # new set of read - store
                                        readid_prev = readid
                                        hit_prev = hit

        # Report pairs of sequences mapped by reads from the same pair
        print(neighbor_l)
        c = Counter(neighbor_l)
        #neighborfile = os.path.join(os.getcwd(), outdir, identifier + 'x' + reference + '.count.learn.txt')
        with open(options.out_learn_file, 'w') as out:
                for k, v in c.items():
                        if v >= options.in_threshold_int:
                                out.write("# bam2counts - neighbor - {} {} {}\n".format(k[0], k[1], v))

        #countfile = os.path.join(outdir, i + 'x' + reference + '.count.txt')
        # Report counts to a file
        with open(options.out_counts_file, 'w') as out:
                out.write("\t{}\n".format(options.in_sample_name))
                for hit in sorted(count_d):
                        out.write("{}\t{}\n".format(hit, count_d[hit]))

        #abundancefile = os.path.join(outdir, identifier + 'x' + reference + '.abundance.txt')

        # Report abundance i.e. counts scale by gene length to file
        with open(options.out_abundance_file, 'w') as out:
                out.write("\t{}\n".format(options.in_sample_name))
                for hit in sorted(count_d):
                        genelength = float(gene_lengths_d[hit])
                        # Calculate Read Per Kilobase : RPK
                        out.write("{}\t{}\n".format(hit, float(count_d[hit]) / (genelength/float(1000))))



