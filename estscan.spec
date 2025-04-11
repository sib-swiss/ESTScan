# $Id: estscan.spec,v 1.6 2007/02/01 15:18:17 c4chris Exp $
Name:           estscan
Version:        3.0.1
Release:        0
Summary:        Detect coding regions in EST sequences

Group:          Applications/Engineering
License:        ESTScan
URL:            http://estscan.sourceforge.net
Source0:        http://dl.sf.net/estscan/%{name}-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

%description
ESTScan is a program that can detect coding regions in DNA sequences, even if
they are of low quality.  ESTScan will also detect and correct sequencing
errors that lead to frameshifts.

ESTScan is not a gene prediction program , nor is it an open reading frame
detector.  In fact, its strength lies in the fact that it does not require an
open reading frame to detect a coding region.  As a result, the program may
miss a few translated amino acids at either the N or the C terminus, but will
detect coding regions with high selectivity and sensitivity.


%package devel
Summary:	Development tools to create matrices for estscan
Group:		Applications/Engineering
Requires:	%{name} = %{version}-%{release}
Provides:	perl(build_model_utils.pl)

%description devel
The estscan-devel package contains various tools to develop and evaluate your
own score matrices for use with estscan.


%prep
%setup -q
sed -i 's+/usr/molbio/share/ESTScan+%{_sysconfdir}/%{name}+' estscan.c
# Help RPM depsolver find the requirements
sed -i 's+/usr/bin/env perl+%{_bindir}/perl+' build_model build_model_utils.pl evaluate_model extract_EST extract_mRNA extract_UG_EST prepare_data


%build
make CFLAGS="-std=gnu99 $RPM_OPT_FLAGS" %{?_smp_mflags} estscan maskred makesmat


%install
rm -rf $RPM_BUILD_ROOT
mkdir -p ${RPM_BUILD_ROOT}%{_bindir}
install -m755 estscan ${RPM_BUILD_ROOT}%{_bindir}
install -m755 maskred ${RPM_BUILD_ROOT}%{_bindir}
install -m755 makesmat ${RPM_BUILD_ROOT}%{_bindir}
install -m755 build_model ${RPM_BUILD_ROOT}%{_bindir}
install -m755 evaluate_model ${RPM_BUILD_ROOT}%{_bindir}
install -m755 extract_EST ${RPM_BUILD_ROOT}%{_bindir}
install -m755 extract_mRNA ${RPM_BUILD_ROOT}%{_bindir}
install -m755 extract_UG_EST ${RPM_BUILD_ROOT}%{_bindir}
install -m755 prepare_data ${RPM_BUILD_ROOT}%{_bindir}

mkdir -p ${RPM_BUILD_ROOT}%{perl_vendorarch}
install -m644 build_model_utils.pl ${RPM_BUILD_ROOT}%{perl_vendorarch}

mkdir -p ${RPM_BUILD_ROOT}%{_sysconfdir}/%{name}


%check


%clean
rm -rf $RPM_BUILD_ROOT


%files
%defattr(-,root,root,-)
%doc COPYRIGHT
%dir %{_sysconfdir}/%{name}/
%{_bindir}/estscan


%files devel
%defattr(-,root,root,-)
%{_bindir}/maskred
%{_bindir}/makesmat
%{_bindir}/build_model
%{_bindir}/evaluate_model
%{_bindir}/extract_EST
%{_bindir}/extract_mRNA
%{_bindir}/extract_UG_EST
%{_bindir}/prepare_data
%{perl_vendorarch}/build_model_utils.pl


%changelog
* Thu Feb  1 2007 Christian Iseli <Christian.Iseli@licr.org> - 3.0.1-0
- version 3.0.1
- 2007-02-01 16:15  c4chris
	* estscan.c, estscan.spec: Bump to version 3.0.1.
- 2007-01-25 15:25  c4chris
	* extract_mRNA: Make use of new BTLib version 0.16 (can now parse
	  general GenBank format hopefully).
- 2007-01-25 14:39  c4chris
	* prepare_data: Properly count nt (skip newlines).

* Tue Dec 19 2006 Christian Iseli <Christian.Iseli@licr.org> - 3.0-0
- created
