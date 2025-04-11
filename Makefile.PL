# $Id: Makefile.PL,v 1.3 2006/12/14 23:24:45 c4chris Exp $
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=> 'ESTScan',
    'VERSION_FROM' => 'ESTScan.pm', # finds $VERSION
    'OBJECT'    => '$(O_FILES)',
    'AUTHOR'    => 'Christian Iseli <Christian.Iseli@licr.org>',
    'EXE_FILES' => [ 'ESTScan' ]
);

sub MY::postamble {
        '
# Make sure we do not overwrite the ESTScan script due to default
# rules of the make command...
ESTScan: ESTScan.o
	touch ESTScan

';
}
