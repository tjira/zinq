#!/usr/bin/env perl

# ENABLE SAFE MODE
use strict; use warnings;

# IMPORT FUNCTIONS
use File::Path   qw(make_path rmtree);
use Cwd          qw(getcwd          );
use File::Copy   qw(move            );
use Getopt::Long qw(GetOptions      );

# DEFINE VARIABLES FOR BUILD OPTIONS
my ($build_eigen, $build_libint, $build_libxc, $build_openblas, $build_fftw);

# DEFINE DEFAULT TARGET
my $target = "x86_64-linux-musl";

# GET THE NUMBER OF CPU CORES
my $cores = qx(nproc --all); chomp $cores; $cores ||= 1;

# PARSE COMMAND-LINE OPTIONS
GetOptions(
    "eigen"    =>    \$build_eigen,
    "libint"   =>   \$build_libint,
    "libxc"    =>    \$build_libxc,
    "openblas" => \$build_openblas,
    "fftw"     =>     \$build_fftw,
    "target=s" =>         \$target,
    "cores|j=i"=>          \$cores,
);

# IF NO FLAGS ARE PROVIDED, BUILD EVERYTHING
if (!$build_eigen && !$build_libint && !$build_libxc && !$build_openblas && !$build_fftw) {
    $build_eigen = $build_libint = $build_libxc = $build_openblas = $build_fftw = 1;
}

# HOST TRIPLE FOR CONFIGURING CROSS-COMPILE
my $host = $target;

# GET CURRENT WORKING DIRECTORY
my $pwd = getcwd();

# EXTRACT OS AND ARCH FROM TARGET TO NAME INSTALLATION DIRECTORY
my $ext_dir = "external-" . ($target =~ /^([^-]+-[^-]+)/)[0];

# CLEAN PREVIOUS INSTALLATION ONLY IF BUILDING EVERYTHING
if ($build_eigen && $build_libint && $build_libxc && $build_openblas && $build_fftw) {
    rmtree($ext_dir) if -d $ext_dir;
}

# DEFINE INSTALL PREFIX
my $prefix = "$pwd/$ext_dir";

# CREATE COMPILER WRAPPERS AND EXPORT ENVIRONMENT VARIABLES
create_compiler_wrappers($target, $pwd);

# COMPILE EACH PROGRAM
compile_eigen   ($prefix, $cores,        $pwd) if    $build_eigen;
compile_libint  ($prefix, $cores,        $pwd) if   $build_libint;
compile_libxc   ($prefix, $cores, $host, $pwd) if    $build_libxc;
compile_openblas($prefix, $cores,        $pwd) if $build_openblas;
compile_fftw    ($prefix, $cores, $host, $pwd) if     $build_fftw;

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
    my %before = map { $_ => 1 } glob("lib/*");

    # DETERMINE TAR FLAGS BASED ON URL
    my $tar_flags = $url =~ /\.tar\.bz2$|\.tbz2$/ ? '-xj' : '-xz';

    # COMMAND TO DOWNLOAD AND EXTRACT THE ARCHIVE
    my $cmd = "curl -L '$url' | tar $tar_flags -C lib";

    # EXECUTE THE COMMAND AND CHECK FOR SUCCESS
    if (system($cmd) != 0) {
        die "FAILED TO DOWNLOAD/EXTRACT '$dest_name' from '$url'";
    }

    # FIND THE EXTRACTED PATH
    my @paths = grep { !$before{$_} } glob("lib/*");

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
    system("cmake", "--install", "build", "--parallel", $cores, "--verbose") == 0 or die "EIGEN INSTALL FAILED";

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_libint {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $pwd) = @_;

    # DEFINE THE URL FOR THE LIBINT SOURCE ARCHIVE
    my $url = "https://github.com/evaleev/libint/releases/download/v2.13.1/libint-2.13.1-mpqc4.tgz";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "libint");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/libint" or die "CANNOT CHDIR TO 'lib/libint': $!";
    
    # CONFIGURE COMMAND
    my @args = (
        "cmake",
        "-B", "build",
        "-DCMAKE_INSTALL_PREFIX=$prefix"
    );

    # RUN CONFIGURE
    system(@args) == 0 or die "LIBINT CONFIGURE FAILED";

    # RUN BUILD
    system("cmake", "--build", "build", "--parallel", $cores, "--verbose") == 0 or die "LIBINT BUILD FAILED";

    # INSTALL THE LIBRARY
    system("cmake", "--install", "build") == 0 or die "LIBINT INSTALL FAILED";

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
        "--prefix=$prefix"
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
    my ($prefix, $cores, $pwd) = @_;

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
        "PREFIX=$prefix"
    );

    # RUN MAKE
    system("make", @args, "-j", $cores, "libs", "shared") == 0 or die "OPENBLAS MAKE FAILED";

    # INSTALL THE LIBRARY
    system("make", @args, "install") == 0 or die "OPENBLAS INSTALL FAILED";

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}

sub compile_fftw {
    # EXTRACT ARGUMENTS
    my ($prefix, $cores, $host, $pwd) = @_;

    # DEFINE THE URL FOR THE FFTW SOURCE ARCHIVE
    my $url = "https://www.fftw.org/fftw-3.3.11.tar.gz";

    # DOWNLOAD AND EXTRACT THE LIBRARY
    download_library($url, "fftw");

    # CHANGE DIRECTORY TO THE EXTRACTED LIBRARY
    chdir "lib/fftw" or die "CANNOT CHDIR TO 'lib/fftw': $!";
    
    # CONFIGURE COMMAND
    my @args = (
        "./configure",
        "--disable-fortran",
        "--host=$host",
        "--prefix=$prefix"
    );

    # RUN CONFIGURE
    system(@args) == 0 or die "FFTW CONFIGURE FAILED";

    # RUN MAKE
    system("make", "-j", $cores) == 0 or die "FFTW MAKE FAILED";

    # INSTALL THE LIBRARY
    system("make", "install") == 0 or die "FFTW INSTALL FAILED";

    # CHANGE BACK TO ORIGINAL DIRECTORY
    chdir $pwd or die "CANNOT CHDIR TO '$pwd': $!";
}
