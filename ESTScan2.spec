# $Id: ESTScan2.spec,v 1.1.1.1 2006/12/13 23:03:45 c4chris Exp $
Name:           ESTScan2
Version:        2.1
Release:        0
Summary:        Detect coding regions in EST sequences

Group:          Applications/Engineering
License:        ESTScan
URL:            http://estscan.sourceforge.net
Source0:        http://dl.sf.net/estscan/%{name}-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)

BuildRequires:  perl
Requires:       perl(:MODULE_COMPAT_%(eval "`%{__perl} -V:version`"; echo $version))

%description
ESTScan is a program that can detect coding regions in DNA sequences, even if
they are of low quality.  ESTScan will also detect and correct sequencing
errors that lead to frameshifts.

ESTScan is not a gene prediction program , nor is it an open reading frame
detector.  In fact, its strength lies in the fact that it does not require an
open reading frame to detect a coding region.  As a result, the program may
miss a few translated amino acids at either the N or the C terminus, but will
detect coding regions with high selectivity and sensitivity.


%prep
%setup -q
sed -i 's+/usr/molbio/share/ESTScan+%{_sysconfdir}/%{name}+' ESTScan
# Help RPM depsolver find the requirements
sed -i 's+/usr/bin/env perl+%{_bindir}/perl+' ESTScan


%build
%{__perl} Makefile.PL INSTALLDIRS=vendor OPTIMIZE="$RPM_OPT_FLAGS"
make %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
make pure_install PERL_INSTALL_ROOT=$RPM_BUILD_ROOT
find $RPM_BUILD_ROOT -type f -name .packlist -exec rm -f {} ';'
find $RPM_BUILD_ROOT -type f -name '*.bs' -a -size 0 -exec rm -f {} ';'
find $RPM_BUILD_ROOT -type d -depth -exec rmdir {} 2>/dev/null ';'
chmod -R u+w $RPM_BUILD_ROOT/*

mkdir -p ${RPM_BUILD_ROOT}%{_sysconfdir}/%{name}
install -m644 Hs.smat ${RPM_BUILD_ROOT}%{_sysconfdir}/%{name}/Hs.smat


%check || :
make test


%clean
rm -rf $RPM_BUILD_ROOT


%files
%defattr(-,root,root,-)
#%doc Changes README
%dir %{_sysconfdir}/%{name}/
%config(noreplace) %{_sysconfdir}/%{name}/Hs.smat
%{_bindir}/ESTScan
%{perl_vendorarch}/*.pm
%{perl_vendorarch}/auto/ESTScan/
%{_mandir}/man3/*.3*


%changelog
* Wed Dec 13 2006 Christian Iseli <Christian.Iseli@licr.org> - 2.1-0
- created
