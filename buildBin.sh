#! /bin/sh
arch=`uname -m`
version=1.4.x-130524

rm CMakeCache.txt 
ccmake -D CMAKE_INSTALL_PREFIX=$HOME/.local -DBUILD_STATIC=YES ..

make 

strip MafFilter/maffilter
tar cvzf maffilter-${arch}-bin-static-${version}.tar.gz MafFilter/maffilter

