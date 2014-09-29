%define _basename maffilter
%define _version 1.1.0
%define _release 1
%define _prefix /usr

URL: http://home.gna.org/bppsuite/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/repos/sources/%{_basename}-%{_version}.tar.gz
Summary: The Bio++ Program Suite
Group: Productivity/Scientific/Other

Requires: libbpp-phyl-omics1 = 2.2.0
Requires: libbpp-seq-omics1 = 2.2.0
Requires: libbpp-phyl9 = 2.2.0
Requires: libbpp-seq9 = 2.2.0
Requires: libbpp-core2 = 2.2.0

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.2.0
BuildRequires: libbpp-core-devel = 2.2.0
BuildRequires: libbpp-seq9 = 2.2.0
BuildRequires: libbpp-seq-devel = 2.2.0
BuildRequires: libbpp-phyl9 = 2.2.0
BuildRequires: libbpp-phyl-devel = 2.2.0
BuildRequires: libbpp-seq-omics1 = 2.2.0
BuildRequires: libbpp-seq-omics-devel = 2.2.0
BuildRequires: libbpp-phyl-omics1 = 2.2.0
BuildRequires: libbpp-phyl-omics-devel = 2.2.0


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

%description
MafFilter is a program for processing and analysing files in the Multiple Alignment Format (MAF).
A description of MAF files can be found on the UCSC genome browser (http://genome.ucsc.edu/FAQ/FAQformat.html#format5).
MafFilter can be used to design a pipeline as a serie of consecutive filters, performing each a dedicated analysis.
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

%changelog
* Fri Sep 26 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.1.0-1
- Initial spec file.
