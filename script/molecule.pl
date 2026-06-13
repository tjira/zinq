#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# IMPORT FUNCTIONS
use File::Path qw(make_path rmtree);
use File::Copy qw(move            );

# DEFINE THE BASIS SET AND METHOD
my $method = $ARGV[0] ||    "RHF";
my $basis  = $ARGV[1] || "STO-3G";

# CLEAN UP TEMP DIRECTORY IF IT EXISTS
rmtree("temp") if -d "temp";

# MAKE THE MOLECULE DIRECTORY IF IT DOESN'T EXIST
make_path("example/molecule") unless -d "example/molecule";

# CREATE TEMP DIRECTORY
make_path("temp") unless -d "temp";

# DEFINE THE MOLECULES
my %molecules = (
    CO2 => [
        "3",
        "CO2",
        "C  0.0000 0.0000 0.0000",
        "O -1.0000 0.0000 0.0000",
        "O  1.0000 0.0000 0.0000"
    ],
    H2 => [
        "2",
        "H2",
        "H 0.0000 0.0000 0.0000",
        "H 1.0000 0.0000 0.0000"
    ],
    HCl => [
        "2",
        "HCl",
        "H  0.0000 0.0000 0.0000",
        "Cl 1.0000 0.0000 0.0000"
    ],
    HF => [
        "2",
        "HF",
        "H 0.0000 0.0000 0.0000",
        "F 1.0000 0.0000 0.0000"
    ],
    Ne2 => [
        "2",
        "Ne2",
        "Ne 0.0000 0.0000 0.0000",
        "Ne 1.0000 0.0000 0.0000"
    ],
    acetone => [
        "10",
        "acetone",
        "O -0.4797 -0.0446  1.2294",
        "C  1.4376  0.0818 -0.1875",
        "C -0.0377 -0.0035  0.0968",
        "C -0.9294 -0.0345 -1.1150",
        "H  2.0471  0.1029  0.7251",
        "H  1.6643  0.9894 -0.7620",
        "H  1.7616 -0.7796 -0.7860",
        "H -0.7844  0.8703 -1.7195",
        "H -0.6856 -0.8986 -1.7467",
        "H -1.9948 -0.0979 -0.8587"
    ],
    water => [
        "3",
        "water",
        "O -0.0669 -0.0232 -0.0267",
        "H  0.6848 -0.6120  0.5191",
        "H  0.3783  0.9803 -0.0938"
    ]
);

# LOOP THROUGH MOLECULES
foreach my $name (keys %molecules) {

    # VARIABLE FOR XYZ FILE
    my $xyz_file = "temp/$name.xyz";

    # CREATE XYZ FILE HANDLE
    open my $fh, '>', $xyz_file or die "CANNOT OPEN '$xyz_file' FOR WRITING: $!";

    # WRITE THE MOLECULE DATA TO THE XYZ FILE
    foreach my $line (@{$molecules{$name}}) {
        print $fh "$line\n";
    }

    # CLOSE THE FILE HANDLE
    close $fh;
}

# CHANGE DIRECTORY TO TEMP
chdir "temp" or die "CANNOT CHDIR TO TEMP: $!";

# PROCESS EACH XYZ FILE
foreach my $mol_file (glob("*.xyz")) {

    # EXTRACT THE BASE NAME (WITHOUT EXTENSION)
    my ($name) = $mol_file =~ /^(.*)\.xyz$/;

    # ORCA INPUT FILE NAME
    my $inp_file = "${name}.inp";

    # CREATE ORCA INPUT FILE HANDLE
    open my $fh_inp, '>', $inp_file or die "Cannot open $inp_file for writing: $!";

    # WRITE ORCA INPUT CONTENT
    print $fh_inp "! $method $basis OPT\n*xyzfile 0 1 $mol_file\n";

    # CLOSE THE ORCA INPUT FILE HANDLE
    close $fh_inp;

    # ORCA OUTPUT FILE NAME
    my $out_file = "${name}.out";

    # ORCA COMMAND
    my $cmd = sprintf('orca %s | tee %s', quotemeta($inp_file), quotemeta($out_file));

    # EXECUTE THE ORCA CALCULATION
    system($cmd);

    # PLACEHOLDER FOR ENERGY VALUE
    my $energy = "";

    # OPEN THE ORCA OUTPUT
    if (open my $fh_out, '<', $out_file) {

        # LOOP OVER LINES
        while (my $line = <$fh_out>) {

            # SEARCH FOR THE FINAL SINGLE POINT ENERGY AND EXTRACT IT
            if ($line =~ /FINAL SINGLE POINT ENERGY\s+(\S+)/) {
                $energy = $1;
            }
        }

        # CLOSE THE ORCA OUTPUT FILE HANDLE
        close $fh_out;
    }

    # TEMP FILE NAME FOR FORMATTING
    my $tmp_file = "${mol_file}.tmp";

    # OPEN THE ORIGINAL XYZ FILE AND THE TEMP FILE FOR WRITING
    open my $fh_in,  '<', $mol_file or die "CANNOT OPEN '$mol_file' FOR READING: $!";
    open my $fh_tmp, '>', $tmp_file or die "CANNOT OPEN '$tmp_file' FOR WRITING: $!";

    # LINE NUMBER TRACKER
    my $line_num = 0;

    # PROCESS EACH LINE OF THE XYZ FILE
    while (my $line = <$fh_in>) {

        # INCREMENT THE LINE NUMBER
        $line_num++;

        # STRIP LEADING AND TRAILING WHITESPACE
        chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;

        # SPLIT THE LINE INTO FIELDS
        my @fields = split /\s+/, $line;

        # IF THIS IS THE FIRST LINE (NUMBER OF ATOMS), PRINT IT AS IS
        if (@fields == 1) {
            print $fh_tmp "$fields[0]\n";
        }

        # IF THIS IS THE SECOND LINE (MOLECULE NAME), PRINT THE NAME, METHOD, BASIS, AND ENERGY
        if ($line_num == 2) {
            print $fh_tmp "${name}/${method}/${basis}/${energy}\n";
        }

        # IF THIS IS AN ATOM LINE (4 FIELDS), PRINT THE ATOM SYMBOL AND COORDINATES
        if (@fields == 4) {
            printf $fh_tmp "%2s % 3.8f % 3.8f % 3.8f\n", $fields[0], $fields[1], $fields[2], $fields[3];
        }
    }

    # CLOSE THE FILE HANDLES
    close $fh_in; close $fh_tmp;
    
    # MOVE THE TEMP FILE TO THE EXAMPLE DIRECTORY
    move($tmp_file, "../example/molecule/$mol_file") or die "CANNOT MOVE '$mol_file' TO EXAMPLE DIRECTORY: $!";
}

# GO BACK TO PARENT DIRECTORY
chdir ".." or die "CANNOT CHDIR TO PARENT DIRECTORY: $!";

# REMOVE TEMP DIRECTORY
rmtree("temp") if -d "temp";
