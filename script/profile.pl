#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# DEFINE VARIABLES
my $root_node;

if (@ARGV >= 2 && $ARGV[0] eq '-z') {
    shift @ARGV; $root_node = shift @ARGV;
}

# CHECK THAT AT LEAST ONE ARGUMENT WAS PROVIDED
if (@ARGV < 1) {
    die "Usage: $0 [-z root_node] <program> [args...]\n";
}

# PROFILE THE PROGRAM
system('valgrind', '-q', '--tool=callgrind', '--callgrind-out-file=callgrind.data', @ARGV);

# GET THE '-z' OPTION IF PROVIDED
my $z_flag = defined $root_node ? "-z '$root_node'" : "";

# CREATE THE ANALYSIS PIPELINE
my $pipeline = "gprof2dot -e 1 -f callgrind -n 5 $z_flag callgrind.data | dot -T svg -o profile.svg";

# RUN THE PIPELINE
system($pipeline);

# CLEAN UP
if (-e 'perf.data') {
    unlink 'perf.data';
}
