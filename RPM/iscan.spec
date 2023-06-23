Summary: Twinscan/N-SCAN Gene Predictor 
Name: iscan
Version: 4.1.2
Release: 1
License: GPL
Group: Bioinformatics
URL: http://mblab.wustl.edu/software/twinscan
Source0: %{name}-%{version}.tar.gz
Buildroot: %{_tmppath}/%{name}-%{version}-root
Provides: pairagon

%description
Twinscan/N-SCAN is a suite of software for gene-structure prediction. Twinscan is currently available for Mammals, Caenorhabditis (worm), Dicot plants, and Cryptococci. N-SCAN is available for human and Drosophila (fruitfly).  This package also contains Pairagon, a PairHMM-based cDNA to genome alignment program.  See http://mblab.wustl.edu for more information.

%prep
%setup -q 

%ifarch x86_64
    %define target linux64
%else
    %define target linux
%endif


%build
make %{target} pairagon-linux

%install
rm -rf $RPM_BUILD_ROOT
make install install-pairagon DESTDIR=$RPM_BUILD_ROOT/usr

%clean
rm -rf $RPM_BUILD_ROOT

%files 
%defattr(-,root,root)
%doc README*
/usr/share/iscan
/usr/lib
/usr/bin

%changelog
* Tue Jul 24 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 4.1.2-1
- Fixed makefile to include BPdeluxe.pm
* Tue Jul 24 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 4.1.1-1
- Adding pairagon README, parameters, examples
- Adding BPdeluxe.pm
* Tue Jun 12 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 4.1-1
- Added summary
- Split off pairagon and iscan scripts in install and install-pairagon targets
- Removed native perl extensions
- Updated README for removal of native perl extensions
- Bumped version for release
- Removed FAlite.pm from the bin directory
- Removed GTF_Parser_SU.pm
- Removed split_zhmm.pl, which is only used for parameter estimation

* Tue Jun 12 2007 Bob Zimmermann <rpz@cse.wustl.edu> - 4.0-1
- Initial build
