#!/bin/bash

rm -rf lib && mkdir lib

CORES=$(nproc --all) && PREFIX=$PWD/external

git clone --depth 1 https://gitlab.com/libxc/libxc.git lib/libxc

echo -e "#!/usr/bin/env bash\n\nzig ar     \"\$@\"" > zigar     && chmod +x     zigar
echo -e "#!/usr/bin/env bash\n\nzig cc     \"\$@\"" > zigcc     && chmod +x     zigcc
echo -e "#!/usr/bin/env bash\n\nzig c++    \"\$@\"" > zigcpp    && chmod +x    zigcpp
echo -e "#!/usr/bin/env bash\n\nzig ranlib \"\$@\"" > zigranlib && chmod +x zigranlib

export CC="$PWD/zigcc"; export CXX="$PWD/zigcpp"; export AR="$PWD/zigar"; export RANLIB="$PWD/zigranlib"

cd lib/libxc && autoreconf -i && ./configure --prefix=$PREFIX && make -j $CORES && make install && cd ../..

rm -rf zigar zigcc zigcpp zigranlib
