# $Id: Makefile.PL,v 1.1.1.1 2006/12/18 15:44:47 c4chris Exp $
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=> 'ESTScan1',
    'VERSION_FROM' => 'ESTScan1.pm', # finds $VERSION
    'OBJECT'    => '$(O_FILES)',
    'AUTHOR'    => 'Christian Iseli <Christian.Iseli@licr.org>',
    'EXE_FILES' => [ 'ESTScan1' ]
);

sub MY::postamble {
        '
# Make sure we do not overwrite the ESTScan script due to default
# rules of the make command...
ESTScan1: ESTScan1.o
	touch ESTScan1

';
}
