# $Id: Makefile.PL,v 1.3 2000/08/22 13:21:35 clottaz Exp $
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
    'LIBS'	=> [''],   # e.g., '-lm' 
    'DEFINE'	=> '',     # e.g., '-DHAVE_SOMETHING' 
    'INC'	=> '',     # e.g., '-I/usr/include/other' 
    'OBJECT'    => '$(O_FILES)',
    'MYEXTLIB'  => 'libestscan$(LIB_EXT)',
);

sub MY::postamble {
        '
all :: static

static :: libestscan$(LIB_EXT)

libestscan$(LIB_EXT): $(O_FILES)
	$(AR) cr libestscan$(LIB_EXT) $(O_FILES)
	$(RANLIB) libestscan$(LIB_EXT)

';
}
