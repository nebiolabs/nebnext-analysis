# -*- sh -*-

# Tiny test files for testing adapter trimming using flexbar. Designed
# to cover major cases of adapter occurrences in NEB single-cell and
# low-input RNA-seq libraries using template-switching oligos.

# descriptions of individual tests:
# test_descriptions.tsv

# Each line describes a read pair such as the one below:

# ./in/r1_ins_1_28_tso_full_g5_rev_1.fastq
# ./in/r1_ins_1_28_tso_full_g5_rev_2.fastq

# Expected output after trimming:
# ./exp/r1_ins_1_28_tso_full_g5_rev_1.fastq
# ./exp/r1_ins_1_28_tso_full_g5_rev_2.fastq

# Adapter files:
# adapters/*.fasta

############################################################

# To run the tests:

# I used flexbar development branch, this version, running on ubuntu
# aws instance:

# commit 93bc8837ee335f5bd6e81266ae247f1e4a61b04c (HEAD -> develop, origin/develop)
# Date:   Thu May 10 17:19:20 2018 +0200

# Ensure that flexbar is in the path, then run test.sh, e.g.,

# /home/ubuntu/src/nebnext-analysis/flexbar/test/test.sh && echo 'ok' || echo 'not ok'
# ok
