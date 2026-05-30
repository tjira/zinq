#!/bin/bash

# ERROR HANDLING
set -euo pipefail

# TARGETS AND HOST
TARGETS=(x86_64-linux-musl) && HOST="x86_64-linux-musl"

# CLEAN PREVIOUS LIBRARIES
rm -rf lib && mkdir lib 

# CLONE LIBRARIES
curl -L https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz | tar -xz -C lib && mv lib/eigen* lib/eigen
curl -L https://github.com/evaleev/libint/releases/download/v2.13.1/libint-2.13.1-mpqc4.tgz | tar -xz -C lib && mv lib/libint* lib/libint
curl -L https://gitlab.com/libxc/libxc/-/archive/7.0.0/libxc-7.0.0.tar.bz2 | tar -xj -C lib && mv lib/libxc* lib/libxc
curl -L https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.33/OpenBLAS-0.3.33.tar.gz | tar -xz -C lib && mv lib/OpenBLAS* lib/openblas
curl -L https://www.fftw.org/fftw-3.3.11.tar.gz | tar -xz -C lib && mv lib/fftw* lib/fftw

# LOOP OVER TARGETS
for TARGET in "${TARGETS[@]}"; do

    # CLEAN PREVIOUS INSTALLATIONS
    rm -rf external-$(echo $TARGET | cut -d- -f1-2)

    # CORES AND PREFIX
    CORES=$(nproc --all) && PREFIX=$PWD/external-$(echo $TARGET | cut -d- -f1-2)

    # EIGEN COMPILATION OPTIONS
    MAKE_EIGEN=(
        -B build
        -DBUILD_SHARED_LIBS=False
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX=$PREFIX
        -DCMAKE_PREFIX_PATH=$PREFIX
        -DEIGEN_BUILD_BLAS=False
        -DEIGEN_BUILD_LAPACK=False
    )

    # LIBINT COMPILATION OPTIONS
    MAKE_LIBINT=(
        -B build
        -DCMAKE_INSTALL_PREFIX=$PREFIX
    )

    # LIBXC COMPILATION OPTIONS
    MAKE_LIBXC=(
        --disable-fortran
        --host=$HOST
        --prefix=$PREFIX
    )

    # OPENBLAS COMPILATION OPTIONS
    MAKE_OPENBLAS=(
        HOSTCC=gcc
        NOFORTRAN=1
        NO_SHARED=1
        NUM_THREADS=128
        PREFIX=$PREFIX
    )

    # CREATE COMPILER WRAPPERS
    echo -e "#!/usr/bin/env bash\n\nzig ar                      \"\$@\"" > zigar     && chmod +x     zigar
    echo -e "#!/usr/bin/env bash\n\nzig cc     --target=$TARGET \"\$@\"" > zigcc     && chmod +x     zigcc
    echo -e "#!/usr/bin/env bash\n\nzig c++    --target=$TARGET \"\$@\"" > zigcpp    && chmod +x    zigcpp
    echo -e "#!/usr/bin/env bash\n\nzig ranlib                  \"\$@\"" > zigranlib && chmod +x zigranlib

    # EXPORT COMPILERS
    export CC="$PWD/zigcc"; export CXX="$PWD/zigcpp"; export AR="$PWD/zigar"; export RANLIB="$PWD/zigranlib"

    # COMPILE EIGEN
    cd lib/eigen && cmake "${MAKE_EIGEN[@]}" && cmake --install build && cd ../..

    # COMPILE LIBINT
    cd lib/libint && cmake "${MAKE_LIBINT[@]}" && cmake --build build --parallel $CORES --verbose && cmake --install build && cd ../..

    # COMPILE LIBXC
    cd lib/libxc && autoreconf -i && ./configure "${MAKE_LIBXC[@]}" && make -j $CORES && make install && cd ../..

    # COMPILE OPENBLAS
    cd lib/openblas && make "${MAKE_OPENBLAS[@]}" -j $CORES libs shared && make "${MAKE_OPENBLAS[@]}" install && cd ../..

    # COMPILE FFTW
    cd lib/fftw && ./configure --host=$HOST --prefix=$PREFIX --disable-fortran && make -j $CORES && make install && cd ../..

    # REMOVE COMPILER WRAPPERS
    rm -rf zigar zigcc zigcpp zigranlib
done
