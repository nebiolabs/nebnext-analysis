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
# adapters/ilmn_13.fasta
# adapters/tso.new.fasta

############################################################

# To run the tests:

# I used flexbar development branch, this version, running on ubuntu
# aws instance:

# commit e4076e5b781cfd84c375c90bc0f6f799b61ecc93 (HEAD -> develop, origin/develop)
# Date:   Mon Apr 9 16:30:10 2018 +0200

# launch aws instance, e.g.:
# AMI ID
# flexbar_dev_4_libs_1e5_1e6_pairs (ami-2707a658)

# IPv4 Public IP
# 54.211.238.166

# connect to aws from rl1 (note user name = ubuntu, not ec2-user):

ssh -AY -i /home/galaxy/.ssh/nextflow.pem.txt ubuntu@54.211.238.166

# on aws, mount storage if needed:

lsblk

sudo mount /dev/xvdb  /data

# on the cluster, rsync test files from cluster to aws instance (start
# instance and mount dest dir before rsync):

rsync -avR -e "ssh -i /home/galaxy/.ssh/nextflow.pem.txt" --include '*/' --include '*' --exclude '*' ${work_dir}/./test ubuntu@54.211.238.166:/data/flex_180404_1

# on aws, run tests from the perl wrapper:

export work_dir=/data/flex_180404_1

export adapter_dir=${work_dir}/test/fixtures/adapters

export fastq_dir=${work_dir}/test/fixtures/in

export exp_dir=${work_dir}/test/fixtures/exp

export obs_dir=${work_dir}/test/fixtures/obs

mkdir -p ${obs_dir}

# optimal run time, higher numbers have no effect:

export num_threads=4

cd ${work_dir}


# Run flexbar 2x: right trim illumina, then left trim tso adapters (as
# of this flexbar version, these two steps cannot be
# combined). Compare expected and observed output fastq files, write
# ok/not ok test results into tsv file, and diagnostic info into log
# file.

cat test_library_names.tsv | \
    perl -lne '
BEGIN {
        @fields = qw( library is_ok );
        print join "\t", @fields;

}
@cmds = ();
%val = ();
$status = 0;
$val{library} = $_;

$cmd = qq{flexbar --reads $ENV{fastq_dir}/$val{library}_1.fastq --reads2 $ENV{fastq_dir}/$val{library}_2.fastq --target $ENV{obs_dir}/$val{library}.end_right.adp_ilmn_13 --adapters $ENV{adapter_dir}/ilmn_13.fasta --adapter-trim-end RIGHT --threads $ENV{num_threads}};
print {*STDERR} "$cmd";
system $cmd;

$cmd = qq{flexbar --reads $ENV{obs_dir}/$val{library}.end_right.adp_ilmn_13_1.fastq --reads2 $ENV{obs_dir}/$val{library}.end_right.adp_ilmn_13_2.fastq --target $ENV{obs_dir}/$val{library} --adapters $ENV{adapter_dir}/tso.new.fasta --adapter-trim-end LEFT --threads $ENV{num_threads}};
print {*STDERR} "$cmd";
system $cmd;

for $mate (1, 2) {
        $cmd = qq{diff $ENV{exp_dir}/$val{library}_${mate}.fastq $ENV{obs_dir}/$val{library}_${mate}.fastq >/dev/null};
        print {*STDERR} "$cmd";
        system $cmd and $status++;
}
$val{is_ok} = $status ? q{not ok} : q{ok};
print join "\t", @val{ @fields };
' \
         1>test.end_right.adp_ilmn_13.end_left.adp_tso.new.tsv \
         2>test.end_right.adp_ilmn_13.end_left.adp_tso.new.log

cat test.end_right.adp_ilmn_13.end_left.adp_tso.new.tsv | column -t -s$'\t'

# Prints this table. See test_descriptions.tsv for details. Note that
# 2 tests fail, this is a known issue (discussing this with Johannes
# Roehr).

# library                            is_ok
# r1_ins_1_28_tso_full_g5_rev        not ok
# r1_ins_1_40_ilmn_13_fwd            ok
# r1_ins_1_40_tso_full_g5_rev        not ok
# r1_ins_1_59_ilmn_13_fwd            ok
# r1_tso_full_g4_fwd_1_43_ins_44_76  ok
# r1_tso_full_t29_fwd_1_66_ins       ok
# r1_tso_full_t39_fwd_1_76_no_ins    ok
# r2_tso_full_g4_fwd_1_43_ins_44_76  ok
# r2_tso_full_t39_fwd_1_76_no_ins    ok

