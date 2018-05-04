#!/usr/bin/env bash

read -d '' usage <<'EOF'

test.sh - runs tests for flexbar adapter trimming.

DESCRIPTION

flexbar must be installed and present in the path to run these
tests. Exit status = 0 is all tests succeeded, and non-zero
otherwise.

Writes files with test details in the current dir:
test_flexbar.tsv
test_flexbar.log

USAGE

test.sh

EXAMPLES OF USAGE

# Run tests and print ok if all tests succeeded:
test.sh && echo 'ok' || echo 'not ok'

# View the detailed results:

cat test_flexbar.tsv | column -t -s$'\t'

# Prints this table. See test_descriptions.tsv for details. 

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

EOF

# export usage='foo'

if [[ "$#" -ne "0" ]] ; then
    echo "${usage}"
    exit 1
fi

export script_dir=`dirname $0`

export adapter_dir=${script_dir}/fixtures/adapters

export fastq_dir=${script_dir}/fixtures/in

export exp_dir=${script_dir}/fixtures/exp

export obs_dir=${script_dir}/fixtures/obs

mkdir -p ${obs_dir}

# optimal run time, higher numbers have no effect:

export num_threads=4


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
         1>test_flexbar.tsv \
         2>test_flexbar.log


test_exit_status=`grep -c 'not ok' test_flexbar.tsv`

exit ${test_exit_status}
