package ESTScan;
#
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

require Exporter;
require DynaLoader;

@ISA = qw(Exporter DynaLoader);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(
	
);
$VERSION = '2.1';
bootstrap ESTScan $VERSION;

sub LoadMatrix {
    my ($file, $type, $correction) = @_;
    local *FH;
    my @res;
    open FH, $file;
    while ( <FH> ) {
	if (/^FORMAT: \S+ (\S+)/ && $type eq $1) {
	    while (1) {
		my @a;
		my($name, $fType, $mType, $order, $period, $offset, $CGmin, $CGmax) =
		    m/^FORMAT: (\S+) (\S+) (\S+) (\d+) (\d+) (\d+) s C\+G: (\S+) (\S+)/;
		while ( <FH> ) {
		    if (/^(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/) {push @a,$1,$2,$3,$4;} 
		    else { last; }
		}
		last if $#a < 0;
		if ($CGmin > 0.0) { $CGmin += $correction; }
		if ($CGmin < 0.0) { $CGmin = 0.0; }
		if ($CGmax < 100.0) { $CGmax += $correction; }
		if ($CGmax > 100.0) { $CGmax = 100.0; }
		push @res, [$name, $fType, $mType, $order, $period, $offset, 
			    $CGmin, $CGmax, \@a, 0];
		if (/^FORMAT: \S+ (\S+)/ && $type eq $1) { next; }
		while ( <FH> ) {
		    if (/^FORMAT: \S+ (\S+)/ && $type eq $1) { last; }
		}
	    }
	    last;
	}
    }
    close FH;
    return @res;
}

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is the stub of documentation for the module.

=head1 NAME

ESTScan - Perl extension for coding region prediction is EST sequences.

=head1 SYNOPSIS

  use ESTScan;

=head1 DESCRIPTION

Better read the ISMB 99 paper...

=head1 AUTHOR

Christian Iseli, LICR ITO, Christian.Iseli@licr.org

Copyright (c) 1999 Swiss Institute of Bioinformatics. All rights reserved.

=head1 SEE ALSO

perl(1).

=cut
