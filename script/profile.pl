#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# CHECK THAT AT LEAST ONE ARGUMENT WAS PROVIDED
if (@ARGV < 1) {
    die "Usage: $0 <program> [args...]\n";
}

# PROFILE THE PROGRAM
system('perf', 'record', '-g', '--', @ARGV);

# CREATE THE ANALYSIS PIPELINE
my $pipeline = 'perf script | gprof2dot -e 1 -f perf -n 5 -z "main.main" | dot -T svg -o profile.svg';

# RUN THE PIPELINE
system($pipeline);

# CLEAN UP
if (-e 'perf.data') {
    unlink 'perf.data';
}
