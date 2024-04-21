#!/usr/bin/python
# -*- coding: utf-8 -*-

###################
# Filter bam file based on mapping identity
#
####################

import sys
import re
import pysam
import os
import os.path
import logging
from collections import deque

from optparse import OptionParser

# Logging file for debugging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def getIdentity(cigar_tuple_list, extra_info_tuple):  # return percentage ID
  cigar = {}
  for i in range(10):
    cigar[i] = 0
  for key, val in cigar_tuple_list:
    cigar[key] = cigar[key] + val

  edit_length = 0
  for key, val in extra_info_tuple:
    if re.search(r'NM', key):
      edit_length = val

  read_length = cigar[0] + cigar[1] + cigar[4] + cigar[5]  # CIGAR identifiers: cigar[0] - M, cigar[1] - I, cigar[4] - S, cigar[5] - H
  identity = (cigar[0] - edit_length) / float(read_length)
  # print 'Read length %d, matches is %d, edit length is %d, is Identity is %.2f' % (read_length, cigar[0], edit_length, identity)
  return identity


def main():

  parser = OptionParser()
  parser.add_option("-i", "--input.bam", dest="in_bam_file", help="in bam file", metavar="FILENAME")
  parser.add_option("-o", "--out.bam", dest="out_bam_file", help="out bam file", metavar="FILENAME")
  parser.add_option("-f", "--id.fraction", dest="id_fraction", help="identity level; recommended is 0.95", metavar="FLOAT")
  (options, args) = parser.parse_args()

  in_bam_file = options.in_bam_file
  out_bam_file = options.out_bam_file
  id_fraction = float(options.id_fraction)

  print("File to read: ", in_bam_file)
  print("Output file: ", out_bam_file)
  print("Identity level: ", str(id_fraction))

  ######################### Start parsing bam file ###################################
  print("Reading in bam file")
  samfileIN = pysam.AlignmentFile(in_bam_file, "rb")
  samfileOUT = pysam.AlignmentFile(out_bam_file, "wb", template=samfileIN)

  '''
  CAPPMANXX170326:7:1101:1121:63024       83      CSM5MCX3.gene.138947    624     40      101M    NM:i:0  MD:Z:101        AS:i:101        XS:i:101
  CAPPMANXX170326:7:1101:1121:63024       163     CSM5MCX3.gene.138947    366     48      101M    NM:i:0  MD:Z:101        AS:i:101        XS:i:96 XA:Z:PSMA269
  '''
  # Make a stack - this way, we always look at the next read mate in the deque. This assumes that we we'll never observe 3 of the same read in a row cause.

  print('Parsing bam file')
  read_stack = deque([])
  i = 0
  for alignedread in samfileIN.fetch(until_eof=True):
    i += 1

    read_stack.append(alignedread)

    # Get a read deque of at least 2 reads
    if len(read_stack) < 2:
      continue
    if i % 1000000 == 0:
      print('Parsed {} reads'.format(i))

    alignedread = read_stack[0]
    alignedread2 = read_stack[1]

    identity = getIdentity(alignedread.cigar, alignedread.tags)
    if identity >= id_fraction:
      AS = alignedread.get_tag('AS')
      XS = alignedread.get_tag('XS')
      if AS > XS:
        samfileOUT.write(alignedread)
        if alignedread.query_name == alignedread2.query_name and alignedread.reference_name == alignedread2.reference_name:
          identity = getIdentity(alignedread2.cigar, alignedread2.tags)
          if identity >= id_fraction:
              samfileOUT.write(alignedread2)
              read_stack.pop()  # Pop a read from the deque
              read_stack.pop()  # Pop a read from the deque
              continue
      else:
      	# The alignment quality of first read was poor, will the mate save it?
        if alignedread.query_name == alignedread2.query_name and alignedread.reference_name == alignedread2.reference_name:
          identity = getIdentity(alignedread2.cigar, alignedread2.tags)
          if identity >= id_fraction:
            AS = alignedread2.get_tag('AS')
            XS = alignedread2.get_tag('XS')
            if AS > XS:
              samfileOUT.write(alignedread)
              samfileOUT.write(alignedread2)
              read_stack.pop()  # Pop a read from the deque
              read_stack.pop()  # Pop a read from the deque
              continue
    read_stack.popleft()

  # We need to collect the reads in pairs, so if the 2nd one doesn't correspond to the first one, we only write on to file.
  # so we should check continuously and always consider two reads but only write out once at the time
  # prevAS is True, for each read, when the AS > XS else False

  print('Done')

if __name__ == '__main__':
  main()
