#! /bin/sh
arch=`uname -m`
version=1.1.2-1

strip MafFilter/maffilter
tar cvzf maffilter-${arch}-bin-static-${version}.tar.gz MafFilter/maffilter

