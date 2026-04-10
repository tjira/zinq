#!/bin/bash

rm -rf external lib && mkdir lib

CORES=$(nproc --all) && PREFIX=$PWD/external

MAKE_EIGEN=(
    -B build
    -DBUILD_SHARED_LIBS=False
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=$PREFIX
    -DCMAKE_PREFIX_PATH=$PREFIX
    -DEIGEN_BUILD_BLAS=False
    -DEIGEN_BUILD_LAPACK=False
)

MAKE_LIBINT_COMPILER=(
    -B build
)

MAKE_LIBINT=(
    -B build
    -DCMAKE_INSTALL_PREFIX=$PREFIX
)

MAKE_LIBXC=(
    --prefix=$PREFIX
)

MAKE_OPENBLAS=(
    HOSTCC=gcc
    NOFORTRAN=1
    NO_SHARED=1
    NUM_THREADS=128
    PREFIX=$PREFIX
)

# CLONE LIBRARIES
git clone https://gitlab.com/libeigen/eigen.git       lib/eigen
git clone https://github.com/evaleev/libint.git       lib/libint
git clone https://gitlab.com/libxc/libxc.git          lib/libxc
git clone https://github.com/OpenMathLib/OpenBLAS.git lib/openblas

# COMPILE LIBINT COMPILER
cd lib/libint && cmake "${MAKE_LIBINT_COMPILER[@]}" && cmake --build build --parallel $CORES --target build_libint --verbose && cmake --build build --target export && cd ../..

# EXTRACT COMPILED LIBINT COMPILER
tar -xzvf lib/libint/build/*.tgz && mv libint* lib/clibint

# CREATE COMPILER WRAPPERS
echo -e "#!/usr/bin/env bash\n\nzig ar                                \"\$@\"" > zigar     && chmod +x     zigar
echo -e "#!/usr/bin/env bash\n\nzig cc     --target=x86_64-linux-musl \"\$@\"" > zigcc     && chmod +x     zigcc
echo -e "#!/usr/bin/env bash\n\nzig c++    --target=x86_64-linux-musl \"\$@\"" > zigcpp    && chmod +x    zigcpp
echo -e "#!/usr/bin/env bash\n\nzig ranlib                            \"\$@\"" > zigranlib && chmod +x zigranlib

# EXPORT COMPILERS
export CC="$PWD/zigcc"; export CXX="$PWD/zigcpp"; export AR="$PWD/zigar"; export RANLIB="$PWD/zigranlib"

# COMPILE EIGEN
cd lib/eigen && cmake "${MAKE_EIGEN[@]}" && cmake --install build && cd ../..

# COMPILE LIBINT
cd lib/clibint && cmake "${MAKE_LIBINT[@]}" && cmake --build build --parallel $CORES --verbose && cmake --install build && cd ../..

# COMPILE LIBXC
cd lib/libxc && autoreconf -i && ./configure "${MAKE_LIBXC[@]}" && make -j $CORES && make install && cd ../..

# COMPILE OPENBLAS
cd lib/openblas && make "${MAKE_OPENBLAS[@]}" -j $CORES libs shared && make "${MAKE_OPENBLAS[@]}" install && cd ../..

# REMOVE COMPILER WRAPPERS
rm -rf zigar zigcc zigcpp zigranlib
