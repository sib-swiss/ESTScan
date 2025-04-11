# $Id: build_model_utils.pl,v 1.10 2006/12/19 16:01:57 c4chris Exp $

################################################################################
#
# build_model_utils
# -----------------
#
# Claudio Lottaz, SIB-ISREC, Claudio.Lottaz@isb-sib.ch
# Christian Iseli, LICR ITO, Christian.Iseli@licr.org
#
# Copyright (c) 1999-2002, 2006 Swiss Institute of Bioinformatics.
# All rights reserved.
#
################################################################################

use strict;
use Symbol;
use POSIX ();

my $reportfh = gensym;
my $local_datadir = '';
my $local_filestem = '';
my $local_verbose = 0;

return 1;

################################################################################
#
#   Define and show paramters
#

sub readConfig {
    # read the configuration specified in the given file
    my($parFileName, $forcedtuplesize, $forcedminmask, $forcedpseudocounts,
       $forcedminscore, $forcedstartlength, $forcedstartpreroll,
       $forcedstoplength, $forcedstoppreroll, $verbose) = @_;

    my($organism, $hightaxo, $dbfiles, $ugdata, $estdata, $datadir, $filestem,
       $rnafile, $estfile, $estcdsfile, $estutrfile,  $trainingfile,
       $testfile, $utrfile, $cdsfile, $tuplesize, $minmask, $pseudocounts,
       $minscore, $startlength, $startpreroll, $stoplength, $stoppreroll,
       $smatfile, $estscanparams, $nb_isochores);

    # default values
    $tuplesize = 6;
    $minmask = 30;
    $pseudocounts = 1;
    $minscore = -100;
    $startpreroll = 2;
    $stoppreroll = 2;
    $nb_isochores = 0;
    $estscanparams = "-m -100 -d -50 -i -50 -N 0";
    my(@isochore_borders) = (0.0, 43.0, 47.0, 51.0, 100.0);

    # read parameter file
    unless (-s $parFileName) {
	warn "The parameter-file $parFileName does not exist, skipped";
	next;
    }
    if (!eval `cat $parFileName`) { die "Error in '$parFileName': $@"; }
    if (!defined($organism)) { die '$organism not specified in '.$parFileName; }
    if (!defined($datadir)) { die '$datadir not specified in '.$parFileName; }

    # set up directories
    if (!(-e $datadir)) { mkdir($datadir, 0775); }
    if (!(-e "$datadir/Report")) { mkdir("$datadir/Report", 0775); }
    if (!(-e "$datadir/Matrices")) { mkdir("$datadir/Matrices", 0775); }
    if (!(-e "$datadir/Isochores")) { mkdir("$datadir/Isochores", 0775); }
    if (!(-e "$datadir/Shuffled")) { mkdir("$datadir/Shuffled", 0775); }
    if (!(-e "$datadir/Evaluate")) { mkdir("$datadir/Evaluate", 0755); }

    # compute further default values and forced values
    if (defined($forcedtuplesize))    { $tuplesize    = $forcedtuplesize;    }
    if (defined($forcedminmask))      { $minmask      = $forcedminmask;      }
    if (defined($forcedpseudocounts)) { $pseudocounts = $forcedpseudocounts; }
    if (defined($forcedminscore))     { $minscore     = $forcedminscore;     }
    if (defined($forcedstartlength))  { $startlength  = $forcedstartlength;  }
    if (defined($forcedstartpreroll)) { $startpreroll = $forcedstartpreroll; }
    if (defined($forcedstoplength))   { $stoplength   = $forcedstoplength;   }
    if (defined($forcedstoppreroll))  { $stoppreroll  = $forcedstoppreroll;  }
    if (!defined($startlength)) {
      $startlength = $startpreroll + POSIX::ceil($tuplesize / 3);
    }
    if (!defined($stoplength)) {
      $stoplength  = $stoppreroll  + POSIX::ceil($tuplesize / 3);
    }
    $filestem =sprintf("%01d_%05d\_%07d\_%1d%1d%1d%1d", $tuplesize, $minmask,
		       $pseudocounts, $startlength, $startpreroll, $stoplength,
		       $stoppreroll);
    $local_datadir = $datadir;
    $local_filestem = $filestem;
    $local_verbose = $verbose;

    if (!defined($rnafile)) {
      $rnafile      = "$datadir/mrna.seq";
    }
    if (!defined($estfile)) {
      $estfile      = "$datadir/ests.seq";
    }
    if (!defined($estcdsfile)) {
      $estcdsfile   = "$datadir/Evaluate/estcds.seq";
    }
    if (!defined($estutrfile)) {
      $estutrfile   = "$datadir/Evaluate/estutr.seq";
    }
    if (!defined($trainingfile)) {
      $trainingfile = "$datadir/training.seq";
    }
    if (!defined($testfile)) {
      $testfile     = "$datadir/test.seq";
    }
    if (!defined($utrfile)) {
      $utrfile      = "$datadir/Evaluate/rnautr.seq";
    }
    if (!defined($cdsfile)) {
      $cdsfile      = "$datadir/Evaluate/rnacds.seq";
    }
    if (!defined($smatfile)) {
      $smatfile     = "$datadir/Matrices/$filestem.smat";
    }

    # check if parameters are plausible
    if ($startpreroll >= $startlength) {
	die("build_model: preroll in start profile too large " .
	    "($startpreroll >= $startlength)");
    }
    if ($stoppreroll > $stoplength) {
	die("build_model: preroll in stop profile too large " .
	    "($stoppreroll >= $stoplength)");
    }
    if ($tuplesize > 11) {
      die "build_model: analysed tuples too big ($tuplesize)";
    }

    return($organism, $hightaxo, $dbfiles, $ugdata, $estdata, $datadir, $filestem,
	   $rnafile, $estfile, $estcdsfile, $estutrfile, $trainingfile,
	   $testfile, $utrfile, $cdsfile, $tuplesize, $minmask, $pseudocounts,
	   $minscore, $startlength, $startpreroll, $stoplength, $stoppreroll,
	   $smatfile, $estscanparams, $nb_isochores, \@isochore_borders);
}

sub showConfig {
    # Show chosen paramters
    my($parFileName, $organism, $hightaxo, $dbfiles, $ugdata, $estdata, $datadir,
       $rnafile, $estfile, $estcdsfile, $estutrfile, $trainingfile, $testfile,
       $utrfile, $cdsfile, $tuplesize, $minmask, $pseudocounts, $minscore,
       $startlength, $startpreroll, $stoplength, $stoppreroll, $smatfile,
       $estscanparams, $nb_isochores, $isochore_borders_ref) = @_;
    my @isochore_borders = @{$isochore_borders_ref};

    my($i);
    log_print("Build ESTScan Tables for $parFileName");
    $parFileName =~ s/\S/\-/g;
    log_print("-------------------------" . $parFileName);
    log_print("\nCurrent parameters:");
    if ($hightaxo eq "") {
      log_print(" - organism:                     $organism");
    } else {
      log_print(" - organism:                     $organism");
      log_print(" - taxonomic level:              $hightaxo");
    }
    log_print(" - database files are:           $dbfiles");
    log_print(" - UniGene data is in:           $ugdata");
    log_print(" - ESTs for testing:             $estdata\n");

    log_print(" - data directory:               $datadir");
    log_print(" - mRNA file is:                 $rnafile");
    log_print(" - EST file is:                  $estfile");
    log_print(" - ESTs with coding:             $estcdsfile");
    log_print(" - ESTs without coding:          $estutrfile");
    log_print(" - training file is:             $trainingfile");
    log_print(" - test file is:                 $testfile");
    log_print(" - clean UTR file is:            $utrfile");
    log_print(" - clean CDS file is:            $cdsfile");
    log_print(" - HMM paramters file:           $smatfile\n");

    log_print(" - tuple size:                   $tuplesize");
    log_print(" - min redundancy mask:          $minmask");
    log_print(" - added pseudocounts:           $pseudocounts");
    log_print(" - minimum score:                $minscore");
    log_print(" - start profile length/preroll: $startlength/$startpreroll");
    log_print(" - stop profile length/preroll:  $stoplength/$stoppreroll");
    if ($nb_isochores>0) {
      log_print(" - nb of isochores:              $nb_isochores");
    } else {
      my @isos;
      for ($i = 0; $i < $#isochore_borders; $i++) {
	push @isos, $isochore_borders[$i] . "-" . $isochore_borders[$i+1];;
      }
      wrapped_log_print(" - Isochores:                    ", 80, @isos);
    }
    log_print(" - options passed scan program:  $estscanparams");
}

sub wrapped_log_print {
    # prints the array given comma separated after its name, wrapping
    # it nicely
    my ($name, $width, @array) = @_;
    my $filler = $name; $filler =~ s/./ /g;

    my $buffer = $name;
    my $linelen = length($name);
    foreach (@array) {
      $_ .= ', ';
      if (($linelen + length($_)) < $width) {
	$buffer .= $_;
	$linelen += length($_);
      }	else {
	$buffer .= "\n$filler$_";
	$linelen = length("$filler$_");
      }
    }
    chop($buffer);chop($buffer);
    log_print($buffer);
}

################################################################################
#
#   Odds and Ends
#

sub log_open {
    my($fname) = "$local_datadir/Report/$local_filestem\_$_[0]";
    open($reportfh, ">$fname");
}

sub log_print {
    my(@stuff) = @_;
    if ($local_verbose) { print @stuff, "\n"; }
    print $reportfh @stuff, "\n";
}

sub log_close {
    close($reportfh);
}

#
#   End of file
#
################################################################################
