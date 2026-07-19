#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# IMPORT FUNCTIONS
use File::Path   qw(make_path rmtree);
use Cwd          qw(getcwd          );
use File::Copy   qw(move copy       );
use Getopt::Long qw(GetOptions      );

# DEFINE LIBRARY NAMES
my @lib_names = qw(eigen libint libxc openblas fftw exprtk);

# DEFINE BUILD OPTIONS
my (%build, $generic);

# DEFINE DEFAULT TARGET AND HOST
my $target = qx(uname -m); chomp $target; $target .= "-$^O"; my $host = $target;

# GET THE NUMBER OF CPU CORES
my $cores = qx(nproc --all); chomp $cores; $cores ||= 1;

# PARSE COMMAND-LINE OPTIONS
GetOptions(
    (map {$_ => \$build{$_}} @lib_names),
    "target=s" => \$target,
    "cores|j=i"=> \$cores,
    "generic"  => \$generic,
);

# IF NO FLAGS ARE PROVIDED, BUILD EVERYTHING
if (!grep { $build{$_} } @lib_names) {
    $build{$_} = 1 for @lib_names;
}

# GET CWD AND PREFIX
my $pwd = getcwd(); my $prefix = "$pwd/external-$target";

# CLEAN PREVIOUS INSTALLATION ONLY IF BUILDING EVERYTHING
if ((grep {$build{$_}} @lib_names) == @lib_names) {
    rmtree($prefix) if -d $prefix;
}

# CREATE COMPILER WRAPPERS AND EXPORT ENVIRONMENT VARIABLES
create_compiler_wrappers($target, $pwd);

# CLEAN LIBRARIES
rmtree("lib") if -d "lib";

# COMPILE EACH PROGRAM
compile_exprtk  ($prefix,                $pwd                   ) if $build{exprtk};
compile_eigen   ($prefix, $cores,        $pwd                   ) if $build{eigen};
compile_libint  ($prefix, $cores,        $pwd                   ) if $build{libint};
compile_libxc   ($prefix, $cores, $host, $pwd                   ) if $build{libxc};
compile_openblas($prefix, $cores,        $pwd,          $generic) if $build{openblas};
compile_fftw    ($prefix, $cores, $host, $pwd,          $generic) if $build{fftw};

# REMOVE COMPILER WRAPPERS
clean_compiler_wrappers();

# CLEAN LIBRARIES
rmtree("lib") if -d "lib";

# LIBRARY DOWNLOAD AND EXTRACTION FUNCTION =============================================================================

sub download_library {
    # EXTRACT ARGUMENTS
    my ($url, $dest_name) = @_;

    # CREATE LIB DIRECTORY IF IT DOESN'T EXIST
    make_path("lib") if ! -d "lib";

    # GET EXISTING PATHS IN LIB TO DETECT NEWLY EXTRACTED CONTENT
    my %before = map {$_ => 1} glob("lib/*");

    # COMMAND TO DOWNLOAD AND EXTRACT THE ARCHIVE USING BSDTAR
    my $cmd = "curl -sL -A 'Mozilla/5.0' '$url' | bsdtar -x -f - -C lib";

    # EXECUTE THE COMMAND AND CHECK FOR SUCCESS
    system($cmd) == 0 or die "FAILED TO DOWNLOAD/EXTRACT '$dest_name' FROM '$url'";

    # FIND THE EXTRACTED PATH
    my @paths = grep {!$before{$_}} glob("lib/*");

    # CHECK IF THE EXTRACTED PATH WAS FOUND
    die "EXTRACTED PATH FOR '$dest_name' NOT FOUND" unless @paths;

    # REMOVE ANY EXISTING DIRECTORY WITH THE DESTINATION NAME
    rmtree("lib/$dest_name") if -d "lib/$dest_name"; 

    # RENAME THE EXTRACTED PATH TO THE DESTINATION NAME
    move($paths[0], "lib/$dest_name") or die "FAILED TO MOVE '$paths[0]' TO 'lib/$dest_name': $!";
}

sub create_compiler_wrappers {
    # EXTRACT ARGUMENTS
    my ($target, $pwd) = @_;

    # DEFINE WRAPPER CONTENTS
    my %wrappers = (
        zigar     => "#!/usr/bin/env bash\n\nzig ar                      \"\$@\"\n",
        zigcc     => "#!/usr/bin/env bash\n\nzig cc     --target=$target \"\$@\"\n",
        zigcpp    => "#!/usr/bin/env bash\n\nzig c++    --target=$target \"\$@\"\n",
        zigranlib => "#!/usr/bin/env bash\n\nzig ranlib                  \"\$@\"\n",
    );

    # LOOP OVER WRAPPERS
    while (my ($file, $content) = each %wrappers) {

        # OPEN FILE FOR WRITING
        open my $fh, '>', $file or die "CANNOT OPEN '$file' FOR WRITING: $!";

        # WRITE CONTENT TO FILE
        print $fh $content;

        # CLOSE FILE
        close $fh;

        # MAKE FILE EXECUTABLE
        chmod 0755, $file or die "CANNOT CHMOD '$file': $!";
    }

    # EXPORT ENVIRONMENT VARIABLES
    $ENV{CC}     =     "$pwd/zigcc";
    $ENV{CXX}    =    "$pwd/zigcpp";
    $ENV{AR}     =     "$pwd/zigar";
    $ENV{RANLIB} = "$pwd/zigranlib";
}

sub clean_compiler_wrappers {
    # DELETE WRAPPER FILES
    unlink qw(zigar zigcc zigcpp zigranlib);

    # UNSET ENVIRONMENT VARIABLES
    delete     $ENV{CC};
    delete    $ENV{CXX};
    delete     $ENV{AR};
    delete $ENV{RANLIB};
}

# LIBRARY COMPILATION FUNCTIONS ========================================================================================

sub compile_eigen {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $pwd) = @_;

    # DEFINE THE URL FOR THE EIGEN SOURCE ARCHIVE
    my $url = "https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "eigen");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/eigen" or die "CANNOT CHDIR TO 'lib/eigen': $!";
    
    # CONFIGURE COMMAND
    my @args = (
        "cmake",
        "-B", "build",
        "-DBUILD_SHARED_LIBS=False",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_INSTALL_PREFIX=$prefix",
        "-DCMAKE_PREFIX_PATH=$prefix",
        "-DEIGEN_BUILD_BLAS=False",
        "-DEIGEN_BUILD_LAPACK=False"
    );

    # RUN CONFIGURE
    system(@args) == 0 or die "EIGEN CONFIGURE FAILED";

    # INSTALL THE LIBRARY
    system("cmake", "--install", "build") == 0 or die "EIGEN INSTALL FAILED";

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_libint {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $pwd) = @_;

    # DEFINE THE URL FOR THE LIBINT SOURCE ARCHIVE
    my $url = "https://github.com/evaleev/libint/archive/refs/tags/v2.13.1.tar.gz";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "libint");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/libint" or die "CANNOT CHDIR TO 'lib/libint': $!";

    # SAVE CURRENT ENVIRONMENT VARIABLES TO RESTORE LATER
    my %saved_env; $saved_env{$_} = delete $ENV{$_} for grep { exists $ENV{$_} } qw(CC CXX AR RANLIB);

    # CONFIGURE LIBINT COMPILER
    my @compiler_args = (
        "cmake",
        "-B", "build",
        "-DLIBINT2_ENABLE_ONEBODY=1",
        "-DLIBINT2_ENABLE_ERI=1",
        "-DLIBINT2_MAX_AM=4"
    );
    system(@compiler_args) == 0 or die "LIBINT COMPILER CONFIGURE FAILED";

    # BUILD LIBINT COMPILER TARGET
    system("cmake", "--build", "build", "--parallel", $cores, "--target", "build_libint") == 0 or die "LIBINT COMPILER BUILD FAILED";

    # RUN EXPORT TARGET TO CREATE THE TGZ ARCHIVE
    system("cmake", "--build", "build", "--target", "export") == 0 or die "LIBINT COMPILER EXPORT FAILED";

    # RESTORE COMPILER ENVIRONMENT VARIABLES FOR TARGET BUILD
    for my $var (keys %saved_env) {
        $ENV{$var} = $saved_env{$var};
    }

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";

    # FIND THE EXPORTED TGZ ARCHIVE
    my ($tgz_file) = glob("lib/libint/build/*.tgz");

    # CHECK IF THE TGZ FILE WAS FOUND
    die "LIBINT EXPORTED ARCHIVE NOT FOUND" unless $tgz_file;

    # EXTRACT THE COMPILED LIBINT COMPILER PACKAGE
    system("tar", "-xzvf", $tgz_file) == 0 or die "FAILED TO EXTRACT LIBINT COMPILER EXPORT";

    # FIND THE EXTRACTED LIBINT DIRECTORY
    my ($extracted_dir) = glob("libint*");

    # CHECK IF THE EXTRACTED DIRECTORY WAS FOUND
    die "EXTRACTED LIBINT DIRECTORY NOT FOUND" unless $extracted_dir;

    # MOVE TO THE TARGET COMPILATION DIRECTORY
    move($extracted_dir, "lib/clibint") or die "FAILED TO MOVE '$extracted_dir' TO 'lib/clibint': $!";

    # CHANGE DIRECTORY TO THE TARGET BUILD PATH
    chdir "lib/clibint" or die "CANNOT CHDIR TO 'lib/clibint': $!";

    # CONFIGURE TARGET LIBINT LIBRARY
    my @args = (
        "cmake",
        "-B", "build",
        "-DCMAKE_INSTALL_PREFIX=$prefix",
        "-DBUILD_SHARED_LIBS=False",
        "-DCMAKE_DISABLE_FIND_PACKAGE_Boost=True"
    );

    # FIX COMPILATION ON WINDOWS
    push @args, "-DLIBINT2_ALIGN_SIZE=0" if $target =~ /windows/;

    # RUN CONFIGURE
    system(@args) == 0 or die "LIBINT CONFIGURE FAILED";

    # RUN BUILD
    system("cmake", "--build", "build", "--parallel", $cores) == 0 or die "LIBINT BUILD FAILED";

    # INSTALL THE LIBRARY
    system("cmake", "--install", "build") == 0 or die "LIBINT INSTALL FAILED";

    # PATCH BOYS.H TO USE WINDOWS ALIGNED MEMORY FUNCTIONS ON MINGW
    if ($target =~ /windows/) {
        system("perl", "-pi", "-e", "s/#ifdef _MSC_VER/#if defined(_MSC_VER) || defined(__MINGW32__)/g", "$prefix/include/libint2/boys.h");
    }

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_libxc {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $host, $pwd) = @_;

    # DEFINE THE URL FOR THE LIBXC SOURCE ARCHIVE
    my $url = "https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.bz2";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "libxc");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/libxc" or die "CANNOT CHDIR TO 'lib/libxc': $!";
    
    # RUN AUTORECONF
    system("autoreconf", "-i") == 0 or die "LIBXC AUTORECONF FAILED";

    # CONFIGURE COMMAND
    my @args = (
        "./configure",
        "--disable-fortran",
        "--host=$host",
        "--prefix=$prefix",
        "--enable-static",
        "--disable-shared"
    );

    # RUN CONFIGURE
    system(@args) == 0 or die "LIBXC CONFIGURE FAILED";

    # RUN MAKE
    system("make", "-j", $cores) == 0 or die "LIBXC MAKE FAILED";

    # INSTALL THE LIBRARY
    system("make", "install") == 0 or die "LIBXC INSTALL FAILED";

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_openblas {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $pwd, $generic) = @_;

    # DEFINE THE URL FOR THE OPENBLAS SOURCE ARCHIVE
    my $url = "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.33/OpenBLAS-0.3.33.tar.gz";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "openblas");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/openblas" or die "CANNOT CHDIR TO 'lib/openblas': $!";
    
    # DEFINE COMPILE ARGUMENTS
    my @args = (
        "HOSTCC=gcc",
        "NOFORTRAN=1",
        "NO_SHARED=1",
        "NUM_THREADS=128",
        "PREFIX=$prefix",
        $generic ? "TARGET=GENERIC" : "DYNAMIC_ARCH=1"
    );

    # MANUALLY PASS OS NAME IF WINDOWS IS SPECIFIED
    push @args, "OSNAME=WINNT" if $target =~ /windows/;

    # RUN MAKE
    system("make", @args, "-j", $cores, "libs", "shared") == 0 or die "OPENBLAS MAKE FAILED";

    # INSTALL THE LIBRARY
    system("make", @args, "install") == 0 or die "OPENBLAS INSTALL FAILED";

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_fftw {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $host, $pwd, $generic) = @_;

    # DEFINE THE URL FOR THE FFTW SOURCE ARCHIVE
    my $url = "https://www.fftw.org/fftw-3.3.11.tar.gz";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "fftw");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/fftw" or die "CANNOT CHDIR TO 'lib/fftw': $!";
    
    # CONFIGURE COMMAND
    my @args = (
        "cmake",
        "-B", "build",
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        "-DBUILD_SHARED_LIBS=OFF",
        "-DBUILD_TESTS=OFF",
        "-DDISABLE_FORTRAN=ON",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_INSTALL_PREFIX=$prefix"
    );

    # RUN CONFIGURE
    system(@args) == 0 or die "FFTW CONFIGURE FAILED";

    # BUILD THE LIBRARY
    system("cmake", "--build", "build", "--parallel", $cores) == 0 or die "FFTW BUILD FAILED";

    # INSTALL THE LIBRARY
    system("cmake", "--install", "build") == 0 or die "FFTW INSTALL FAILED";

    # PATCH FFTW3.H TO DISABLE __float128 ON WINDOWS
    if ($target =~ /windows/ and -f "$prefix/include/fftw3.h") {
        system("perl", "-pi", "-e", "s/#if \\(__GNUC__ > 4/#if 0 && (__GNUC__ > 4/g", "$prefix/include/fftw3.h");
    }

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_exprtk {
    # EXTRACT ARGUMENTS
    my ($prefix, $pwd) = @_;

    # DEFINE THE URL OF THE EXPRTK SOURCE ZIP FILE
    my $url = "https://www.partow.net/downloads/exprtk.zip";

    # CREATE LIB DIRECTORY IF IT DOESN'T EXIST
    make_path("lib") if ! -d "lib";

    # DEFINE THE PATH TO THE ZIP FILE
    my $zip_file = "lib/exprtk.zip";

    # DOWNLOAD THE ZIP FILE WITH A USER-AGENT HEADER
    system("curl", "-L", "-A", "Mozilla/5.0", "--retry", "3", "-o", $zip_file, $url) == 0 or die "FAILED TO DOWNLOAD EXPRTK";

    # UNPACK THE ZIP FILE
    system("unzip", "-q", "-d", "lib", $zip_file) == 0 or die "FAILED TO UNPACK EXPRTK";

    # REMOVE THE ZIP FILE
    unlink $zip_file or warn "COULD NOT REMOVE '$zip_file': $!";

    # CREATE INCLUDE DIRECTORY IF IT DOESN'T EXIST
    make_path("$prefix/include") if ! -d "$prefix/include";

    # DEFINE SOURCE AND DESTINATION PATHS FOR THE HEADER FILE
    my $src = "lib/exprtk/exprtk.hpp"; my $dst = "$prefix/include/exprtk.hpp";

    # RUN COPY
    copy($src, $dst) or die "FAILED TO COPY '$src' TO '$dst': $!";
}
