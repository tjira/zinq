#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# IMPORT FUNCTIONS
use File::Path qw(make_path);

# DEFINE THE BASIS SETS
my @bases = (
    "sto-2g",
    "sto-3g",
    "sto-4g",
    "sto-5g",
    "sto-6g",
    "3-21g",
    "6-21g",
    "6-31g",
    "6-31g*",
    "6-31g**",
    "6-31+g",
    "6-31+g*",
    "6-31+g**",
    "6-31++g",
    "6-31++g*",
    "6-31++g**",
    "6-311g",
    "6-311g*",
    "6-311g**",
    "6-311+g",
    "6-311+g*",
    "6-311+g**",
    "6-311++g",
    "6-311++g*",
    "6-311++g**",
    "def2-svp",
    "def2-svpd",
    "def2-tzvp",
    "def2-tzvpd",
    "def2-tzvpp",
    "def2-tzvppd",
    "def2-qzvp",
    "def2-qzvpd",
    "def2-qzvpp",
    "def2-qzvppd",
    "cc-pvdz",
    "cc-pvtz",
    "cc-pvqz",
    "cc-pv5z",
    "aug-cc-pvdz",
    "aug-cc-pvtz",
    "aug-cc-pvqz",
    "aug-cc-pv5z",
);

# CREATE OUTPUT DIRECTORY
make_path("example/basis") unless -d "example/basis";

# RETRIEVE AND SAVE EACH BASIS SET
foreach my $basis (@bases) {

    # REPLACE SPECIAL CHARACTERS
    my $safe_basis = $basis;
    $safe_basis =~ s/\*/s/g;
    $safe_basis =~ s/\+/p/g;

    # DEFINE OUTPUT FILE NAME
    my $outfile = "example/basis/$safe_basis.g94";

    # CONSTRUCT THE COMMAND
    my $cmd = sprintf('bse get-basis %s gaussian94 > %s', quotemeta($basis), quotemeta($outfile));

    # EXECUTE THE COMMAND
    system($cmd);
}
