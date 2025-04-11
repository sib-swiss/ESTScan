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
$VERSION = '0.01';
bootstrap ESTScan $VERSION;

sub LoadMatrix {
    my ($file, $type, $correction) = @_;
    local *FH;
    my @res;
    open FH, $file;
    my ($name, $fType, $mType, $order, $period, $offset, $CGmin, $CGmax);
    while ( <FH> ) {
	if (/^FORMAT: \S+ (\S+)/ && $type eq $1) {
	    ($name, $fType, $mType, $order, $period, $offset, $CGmin, $CGmax)
		= $_ =~ /^FORMAT: (\S+) (\S+) (\S+) (\d+) (\d+) (\d+) s C\+G: (\S+) (\S+)/;
	    while (1) {
		my @a;
		while ( <FH> ) {
		    if (/^(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)/) {
			push @a, $1, $2, $3, $4;
		    } else {
			last;
		    }
		}
		last if $#a < 0;
		if ($CGmin > 0.0) {
		    $CGmin += $correction;
		    if ($CGmin < 0.0) {
			$CGmin = 0.0;
		    }
		}
		if ($CGmax < 100.0) {
		    $CGmax += $correction;
		    if ($CGmax > 100.0) {
			$CGmax = 100.0;
		    }
		}
		push @res, [$name, $fType, $mType, $order, $period, $offset, $CGmin, $CGmax, \@a, 0];
		if (/^FORMAT: \S+ (\S+)/ && $type eq $1) {
		    ($name, $fType, $mType, $order, $period, $offset, $CGmin, $CGmax)
			= $_ =~ /^FORMAT: (\S+) (\S+) (\S+) (\d+) (\d+) (\d+) s C\+G: (\S+) (\S+)/;
		    next;
		}
		while ( <FH> ) {
		    if (/^FORMAT: \S+ (\S+)/ && $type eq $1) {
			($name, $fType, $mType, $order, $period, $offset, $CGmin, $CGmax)
			    = $_ =~ /^FORMAT: (\S+) (\S+) (\S+) (\d+) (\d+) (\d+) s C\+G: (\S+) (\S+)/;
			last;
		    }
		}
	    }
	    last;
	}
    }
    close FH;
    return @res;
}

sub LoadFPRMatrices {
    my ($FPRMat, $startFlag) = @_;
    my @FPRs;
    local *FH;
    open FH, $FPRMat or die "Couldn't open file $FPRMat: $!";
    if ($startFlag == 1) {
	while (<FH>) { if (m/\-start/) { last; } }
    }
    while ( <FH> ) {
	if (m"^Iso/Length") {
	    s,^Iso/Length\s+\|\s+,,;
	    @FPRs = split;
	    last;
	}
    }
    my $curMin = -1.0;
    my $curMax = -1.0;
    my $curMat;
    while ( <FH> ) {
	if (m"^Iso/Length") { last; }
	if (s/^([\d.]+)-([\d.]+)\s+(\d+)\s+\|\s+//) {
	    my $min = $1;
	    my $max = $2;
	    my $len = $3;
	    my @a = split;
	    if ($min != $curMin || $max != $curMax) {
		# Create new matrix
		my %mat;
		foreach my $a (@FPRs) {
		    $mat{$a} = [];
		}
		$curMat = \%mat;
		push @main::FPRMatrices, [$min, $max, $curMat];
		$curMin = $min;
		$curMax = $max;
	    }
	    foreach my $a (@FPRs) {
		my $val = shift @a;
		die "Not enough values" unless defined $val;
		push @ {$$curMat{$a}}, [$len, $val];
	    }
	    die "Too many values" unless $#a == -1;
	}
    }
    close FH;
}

sub ComputeCutOff {
    my ($percent, $length, $FPR) = @_;
    my $matrix;
    # Determine correct matrix
    foreach my $a (@main::FPRMatrices) {
	if ($percent >= $$a[0] && $percent <= $$a[1]) {
	    $matrix = $$a[2];
	    last;
	}
    }
    die "No false positive rate matrix found for $percent %CG"
	unless defined $matrix;
    my @rates = sort {$a <=> $b} keys %$matrix;
    my ($r0, $r1);
    if ($FPR <= $rates[0]) {
	$r0 = $rates[0];
	$r1 = $rates[1];
    } elsif ($FPR >= $rates[$#rates]) {
	$r0 = $rates[$#rates - 1];
	$r1 = $rates[$#rates];
    } else {
	foreach my $k (@rates) {
	    if ($k <= $FPR) {
		$r0 = $k;
	    }
	    if ($k >= $FPR) {
		$r1 = $k;
		last;
	    }
	}
    }
    my $c0 = $$matrix{$r0};
    my $c1 = $$matrix{$r1};
    my ($l0, $l1);
    my ($s00, $s01, $s10, $s11);
    my $l = $ {$$c0[0]}[0];
    my $maxLIdx = $#$c0;
    if ($length <= $ {$$c0[0]}[0]) {
	($l0, $s00) = @ {$$c0[0]};
	($l1, $s10) = @ {$$c0[1]};
	($l0, $s01) = @ {$$c1[0]};
	($l1, $s11) = @ {$$c1[1]};
    } elsif ($length >= $ {$$c0[$maxLIdx]}[0]) {
	($l0, $s00) = @ {$$c0[$maxLIdx - 1]};
	($l1, $s10) = @ {$$c0[$maxLIdx]};
	($l0, $s01) = @ {$$c1[$maxLIdx - 1]};
	($l1, $s11) = @ {$$c1[$maxLIdx]};
    } else {
	for my $i (0 .. $maxLIdx) {
	    if ($length >= $ {$$c0[$i]}[0]) {
		next;
	    }
	    ($l0, $s00) = @ {$$c0[$i - 1]};
	    ($l1, $s10) = @ {$$c0[$i]};
	    ($l0, $s01) = @ {$$c1[$i - 1]};
	    ($l1, $s11) = @ {$$c1[$i]};
	    last;
	}
    }
#    print "($percent, $length, $FPR): "
#	. "$l0, $l1, $r0, $r1, $s00, $s10, $s01, $s11\n";
    # Do a linear interpolation on the length
    my $lDiff = $l1 - $l0;
    if ($lDiff != 0) {
	$s00 = (($s10 - $s00) / $lDiff) * ($length - $l0) + $s00;
	$s01 = (($s11 - $s01) / $lDiff) * ($length - $l0) + $s01;
    }
#    print "($percent, $length, $FPR): $r0, $r1, $s00, $s01\n";
    # Do a log interpolation on the FPR.
    my $logFPR = log $FPR;
    my $logr0 = log $r0;
    my $logDiff = log($r1) - $logr0;
    if ($logDiff != 0.0) {
	$s00 = (($s01 - $s00) / $logDiff) * ($logFPR - $logr0) + $s00;
    }
#    print "($percent, $length, $FPR): $s00, $logFPR, $logr0, $logDiff\n";
    return int $s00;
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
