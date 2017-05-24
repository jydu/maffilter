%define _basename maffilter
%define _version 1.2.0
%define _release 1
%define _prefix /usr

URL: http://bioweb.me/maffilter

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/repos/sources/maffilter/%{_basename}-%{_version}.tar.gz
Summary: The Multiple Alignment Format file processor
Group: Productivity/Scientific/Other

Requires: libbpp-phyl-omics1 = 2.3.0
Requires: libbpp-seq-omics1 = 2.3.0
Requires: libbpp-phyl9 = 2.3.0
Requires: libbpp-seq9 = 2.3.0
Requires: libbpp-core2 = 2.3.0
Requires: zlib
Requires: libbz2

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.3.0
BuildRequires: libbpp-core-devel = 2.3.0
BuildRequires: libbpp-seq9 = 2.3.0
BuildRequires: libbpp-seq-devel = 2.3.0
BuildRequires: libbpp-phyl9 = 2.3.0
BuildRequires: libbpp-phyl-devel = 2.3.0
BuildRequires: libbpp-seq-omics1 = 2.3.0
BuildRequires: libbpp-seq-omics-devel = 2.3.0
BuildRequires: libbpp-phyl-omics1 = 2.3.0
BuildRequires: libbpp-phyl-omics-devel = 2.3.0
BuildRequires: zlib-devel
BuildRequires: libbz2-devel


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion}
%if 0%{?mdkversion} >= 201100
BuildRequires: xz
%define zipext xz
%else
BuildRequires: lzma
%define zipext lzma
%endif
%else
BuildRequires: gzip
%define zipext gz
%endif

#Distribution specific requirements:

%if 0%{?suse_version} >= 1200 && 0%{?suse_version} < 1300 
Requires: libboost_iostreams1_49_0 
BuildRequires: boost-devel
%else
%if 0%{?suse_version} >= 1300
Requires: libboost_iostreams1_53_0 
BuildRequires: boost-devel
%else
Requires: libboost-iostreams 
BuildRequires: boost-devel
%endif
%endif


%description
MafFilter is a program for processing and analysing files in the Multiple Alignment Format (MAF).
A description of MAF files can be found on the UCSC genome browser (http://genome.ucsc.edu/FAQ/FAQformat.html#format5).
MafFilter can be used to design a pipeline as a series of consecutive filters, each performing a dedicated analysis.
Many filters are available, from alignment cleaning to phylogeny reconstruction and population genetics analysis.

%prep
%setup -q

%build
CFLAGS="-I%{_prefix}/include $RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix}"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
if [ %{zipext} == 'lzma' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=lzma -DDOC_COMPRESS_EXT=lzma"
fi
if [ %{zipext} == 'xz' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DDOC_COMPRESS=xz -DDOC_COMPRESS_EXT=xz"
fi

cmake $CMAKE_FLAGS .
make
make info

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/maffilter
%{_prefix}/share/man/man1/maffilter.1.%{zipext}
%{_prefix}/share/info/maffilter.info.%{zipext}

%changelog
* Wed May 24 2017 Julien Dutheil <dutheil@evolbio.mpg.de 1.2.0-1
* Fri Sep 26 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.1.0-1
- Initial spec file.
