#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# IMPORT FUNCTIONS
use File::Basename qw(dirname   );
use Cwd            qw(abs_path  );
use File::Temp     qw(tempdir   );
use Getopt::Long   qw(GetOptions);

# DEFINE TEMPLATES FOR EACH METHOD
my %TEMPLATES = (
    rhf => {
        orca => "! NoRI RHF {{BASIS}} TightSCF NoAutoStart\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "hartree_fock": {
                "system": "{{XYZ}}",
                "basis": "builtin:{{BASIS}}",
                "integral_direct": true
            }
        }
    ]
}
EOF
    },
    mp2 => {
        orca => "! NoRI MP2 {{BASIS}} TightSCF NoFrozenCore NoAutoStart\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "moller_plesset": {
                "hartree_fock": {
                    "system": "{{XYZ}}",
                    "basis": "builtin:{{BASIS}}"
                },
                "order": 2
            }
        }
    ]
}
EOF
    },
    mp3 => {
        orca => "! NoRI MP3 {{BASIS}} TightSCF NoFrozenCore NoAutoStart\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "moller_plesset": {
                "hartree_fock": {
                    "system": "{{XYZ}}",
                    "basis": "builtin:{{BASIS}}"
                },
                "order": 3
            }
        }
    ]
}
EOF
    },
    cisd => {
        orca => "! NoRI CISD {{BASIS}} TightSCF NoFrozenCore NoAutoStart\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "configuration_interaction": {
                "hartree_fock": {
                    "system": "{{XYZ}}",
                    "basis": "builtin:{{BASIS}}"
                },
                "space": {
                    "excitations": [1, 2]
                },
                "nstate": 1
            }
        }
    ]
}
EOF
    },
    casscf => {
        orca => "! NoRI CASSCF {{BASIS}} TightSCF NoAutoStart\n%casscf\n  nel 2\n  norb 2\nend\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "casscf": {
                "hartree_fock": {
                    "system": "{{XYZ}}",
                    "basis": "builtin:{{BASIS}}"
                },
                "active_electrons": 2,
                "active_orbitals": 2,
                "nstate": 1
            }
        }
    ]
}
EOF
    },
    dft_b3lyp => {
        orca => "! NoRI B3LYP/G {{BASIS}} TightSCF NoAutoStart\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "hartree_fock": {
                "system": "{{XYZ}}",
                "basis": "builtin:{{BASIS}}",
                "dft": {
                    "exchange_correlation": "hyb_gga_xc_b3lyp"
                }
            }
        }
    ]
}
EOF
    },
    dft_pbe => {
        orca => "! NoRI PBE {{BASIS}} TightSCF NoAutoStart\n* xyzfile 0 1 {{XYZ}}\n",
        zinq => <<'EOF',
{
    "zinq": [
        {
            "hartree_fock": {
                "system": "{{XYZ}}",
                "basis": "builtin:{{BASIS}}",
                "dft": {
                    "exchange": "gga_x_pbe",
                    "correlation": "gga_c_pbe"
                }
            }
        }
    ]
}
EOF
    },
);

# DEFINE METHOD ORDER
my @methods = qw(rhf mp2 mp3 cisd casscf dft_b3lyp dft_pbe);

# DEFINE CLI OPTIONS
my $keep = 0; GetOptions("keep" => \$keep);

# DEFINE REPOSITORY ROOT
my $root = abs_path(dirname(abs_path($0)) . "/..");

# DEFINE EXECUTABLES
my $zinq = $ENV{ZINQ} || (-x "$root/zig-out/bin/zinq" ? "$root/zig-out/bin/zinq" : "zinq"); my $orca = $ENV{ORCA} || "orca";

# DEFINE TARGETS
my @targets = @ARGV ? @ARGV : map { "${_}_water_sto-3g" } @methods;

# CREATE WORK DIRECTORY
my $workdir = tempdir(CLEANUP => !$keep);

# INITIALIZE TABLE
my @table;

# EXECUTE BENCHMARKS
for my $t (@targets) {
    # PARSE TARGET STRING
    my ($method, $system, $basis) = ($t =~ /^(.+)_([^_]+)_([^_]+)$/);

    # CHECK IF TARGET FORMAT IS VALID
    die "INVALID TARGET FORMAT '$t'\n" unless defined $method && defined $system && defined $basis;

    # CHECK IF METHOD IS SUPPORTED
    die "UNKNOWN METHOD '$method'\n"   unless exists $TEMPLATES{$method};

    # FIND XYZ FILE
    my $xyz = find_xyz($root, $system);

    # CHECK IF XYZ FILE EXISTS
    die "XYZ FILE FOR '$system' NOT FOUND IN '$root/example/molecule'\n" unless $xyz;

    # PREPARE TEMPLATES
    my $zinq_cfg = $TEMPLATES{$method}{zinq};
    my $orca_cfg = $TEMPLATES{$method}{orca};

    # XYZ SUBSTITUTION
    $zinq_cfg =~ s/\{\{XYZ\}\}/$xyz/g;
    $orca_cfg =~ s/\{\{XYZ\}\}/$xyz/g;

    # BASIS SUBSTITUTION
    $zinq_cfg =~ s/\{\{BASIS\}\}/$basis/g;
    $orca_cfg =~ s/\{\{BASIS\}\}/$basis/g;

    # DEFINE ZINQ INPUT FILE
    my $zinq_in  = "$workdir/zinq_${t}.json";

    # DEFINE ZINQ OUTPUT FILE AND ORCA INPUT/OUTPUT FILES
    my $orca_in  = "$workdir/orca_${t}.inp";
    my $zinq_out = "$workdir/zinq_${t}.out";
    my $orca_out = "$workdir/orca_${t}.out";

    # WRITE INPUT FILES
    write_file($zinq_in, $zinq_cfg);
    write_file($orca_in, $orca_cfg);

    # PRINT PROGRESS
    print "Method: $method, System: $system, Basis: $basis\n";

    # RUN ZINQ
    system("$zinq \"$zinq_in\" > \"$zinq_out\" 2>&1");

    # RUN ORCA IN WORKDIR TO CONTAIN INTERMEDIATE FILES
    system("cd \"$workdir\" && $orca \"$orca_in\" > \"$orca_out\" 2>&1");

    # PARSE OUTPUTS
    my ($ez, $tz) = parse_zinq_output($zinq_out);
    my ($eo, $to) = parse_orca_output($orca_out);

    # COMPUTE ENERGY DIFFERENCE
    my $de = (defined $ez && defined $eo) ? ($ez - $eo) : 0;

    # PRINT LIVE RESULT ENERGIES
    printf "  ZINQ: %14.8f Eh (%s)\n", $ez // 0, $tz // "N/A";
    printf "  ORCA: %14.8f Eh (%s)\n", $eo // 0, $to // "N/A";

    # PRINT LIVE RESULT DIFFERENCE
    printf "  DE:   %14.8f Eh\n\n",   $de;

    # GENERATE FORMATTED TABLE ROW
    my $row = sprintf("%-12s %-8s %-10s %16.8f %-12s %16.8f %-10s %16.8f", $method, $system, $basis, $ez // 0, "($tz)", $eo // 0, "($to)", $de);

    # STORE ROW IN TABLE
    push @table, $row;
}

# PRINT TABLE HEADER
printf "%-12s %-8s %-10s %16s %-12s %16s %-10s %16s\n", "Method", "System", "Basis", "ZINQ (Eh)", "", "ORCA (Eh)", "", "DE (Eh)";

# PRINT TABLE SEPARATOR
print "-----------------------------------------------------------------------------------------------------------\n";

# PRINT TABLE ROWS
print "$_\n" for @table;

# HELPER FUNCTIONS =====================================================================================================

sub find_xyz {
    # EXTRACT ARGUMENTS
    my ($root, $system) = @_;

    # FIND CASE-INSENSITIVE MATCHING MOLECULE FILE
    my ($file) = grep { lc($_) eq lc("$root/example/molecule/$system.xyz") } glob("$root/example/molecule/*");

    # RETURN FOUND PATH
    return $file;
}

sub write_file {
    # EXTRACT ARGUMENTS
    my ($path, $content) = @_;

    # OPEN FILE FOR WRITING
    open my $fh, '>', $path or die "CANNOT WRITE '$path': $!\n";

    # WRITE CONTENT
    print $fh $content;

    # CLOSE FILE
    close $fh;
}

sub parse_zinq_output {
    # EXTRACT ARGUMENTS
    my ($path) = @_;

    # READ FILE CONTENT
    open my $fh, '<', $path or return (undef, undef); my $out = do { local $/; <$fh> }; close $fh;

    # EXTRACT LAST REPORTED ENERGY
    my ($ez) = $out =~ /.*(?:FINAL|MP0).*?\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+Eh/s;

    # EXTRACT EXECUTION TIME
    my ($tz) = $out =~ /TOTAL EXECUTION TIME:\s*(\S+)/;

    # RETURN PARSED VALUES
    return ($ez, $tz);
}

sub parse_orca_output {
    # EXTRACT ARGUMENTS
    my ($path) = @_;

    # READ FILE CONTENT
    open my $fh, '<', $path or return (undef, undef); my $out = do { local $/; <$fh> }; close $fh;

    # EXTRACT LAST SINGLE POINT ENERGY
    my ($eo) = $out =~ /.*FINAL SINGLE POINT ENERGY\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)/s;

    # EXTRACT RUN TIME BREAKDOWN
    my ($d, $h, $m, $s, $ms) = $out =~ /TOTAL RUN TIME:\s*(\d+) days\s*(\d+) hours\s*(\d+) minutes\s*(\d+) seconds\s*(\d+) msec/i;

    # CALCULATE TOTAL ELAPSED SECONDS
    my $sec = defined $ms ? ($d * 86400 + $h * 3600 + $m * 60 + $s) : 0;

    # FORMAT DURATION STRING
    my $to = defined $ms ? ($sec == 0 ? "${ms}ms" : sprintf("%.3fs", $sec + $ms / 1000)) : undef;

    # RETURN PARSED VALUES
    return ($eo, $to);
}
