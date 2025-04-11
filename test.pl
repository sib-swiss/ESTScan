# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..26\n"; }
END {print "not ok 1\n" unless $loaded;}
use ESTScan;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

my $cdsIndex = 0;
my @mat = ESTScan::LoadMatrix("Hs.smat", "CODING", 4.0);
my $utrIndex = $#mat + 1;
push @mat, ESTScan::LoadMatrix("Hs.smat", "UNTRANSLATED", $main::percent);
my $startIndex = $#mat + 1;
push @mat, ESTScan::LoadMatrix("Hs.smat", "START", $main::percent);
my $stopIndex = $#mat + 1;
push @mat, ESTScan::LoadMatrix("Hs.smat", "STOP", $main::percent);
print $#mat == 15 ? "ok 2\n" : "not ok 2\n";
my @EM = (
  [6, 3, 1],
  [6, 3, 1],
  [6, 3, 1],
  [6, 3, 1],
  [6, 1, 1],
  [6, 1, 1],
  [6, 1, 1],
  [6, 1, 1],
  [1, 12, 7],
  [1, 12, 7],
  [1, 12, 7],
  [1, 12, 7],
  [1, 12, 6],
  [1, 12, 6],
  [1, 12, 6],
  [1, 12, 6]
);
for my $i (0 .. $#mat) {
    my $matType = -1;
    if ($ {$mat[$i]}[1] eq 'CODING')       { $matType = 0; }
    if ($ {$mat[$i]}[1] eq 'UNTRANSLATED') { $matType = 1; }
    if ($ {$mat[$i]}[1] eq 'START')        { $matType = 2; }
    if ($ {$mat[$i]}[1] eq 'STOP')         { $matType = 3; }
    my $order = $ {$mat[$i]}[3];
    my $length = $ {$mat[$i]}[4];
    my $offset = $ {$mat[$i]}[5];
    my $CGmin = int($ {$mat[$i]}[6]);
    my $CGmax = int($ {$mat[$i]}[7]);
    my $ar = $EM[$i];
    my $aRef = $ {$mat[$i]}[8];
    my $ref = ESTScan::CreateMatrix($matType, $order, $length, $offset,
				    $CGmin, $CGmax, $aRef, -20, 0);
    $ {$mat[$i]}[9] = $ref;
    if ($ref == $i
	&& $order == $$ar[0]
	&& $length == $$ar[1]
	&& $offset == $$ar[2]) {
      print "ok ", $i + 3, "\n";
    } else {
      print "not ok ", $i + 3, "\n";
    }
}
ESTScan::StoreTransits(-10, -10, -5, -80, -40, -80, -40, -20);

my @res;
my $seq = "TTTTGGTGGTTAGCTCCTCTTGCCAAACCAGCCATGAGCTCCCAGATTCG
TCAGAATTATTTCACCGACGTGGAGGCAGCCGTCAACAGCCTGGTCAATT
TGTACCTGCAGGCCTCCTACACCTACCTCTCTCTGGGCTTCTATTTCGAC
CGCGATGATGTGGCTCTGGAAGGCGTGAGCCACTTCTTCCGCGAACTGGC
CGAGGAGAAGCGCGAGGGCTACGAGCGTCTCCTGAAGATGCAAAACCAGC
GTGGGGCCGCGCTCTCTTCCAGGACATCAAGAAGCCAGCTGAAGATGAGT
GGGTAAAACCCCAGACGCCATGAAAGCTGCCATGGCCCTGGAGAAAAAGC
TGAACCAGGCCCTTTTGGATCTTCATGCCCTGGGTTCTGCCCGCACGGAC
CCCCATCTCTGTGACTTCCTGGAGACTCACTTCCTAGATGAGGAAGTGAA
GCTTATCAAGAAGATGGGTGACCACCTGACCAACCTCCACAGGCTGGGTG
GCCCGGAGGCTGGGCTGGGCGAGTATCTCTTCGAAAGGCTCACTCTCAAG
CACGACTAAGAGCCTTCTGAGCCCAGCGACTTCTGAAGGGCCCCTTGCAA
AGTAATAGGCTTCTGCCTAAGCCTCTCCCTCCGGCCAATAGGCAGCTTTC
TTAACTATCCTAACAAGCCTTGGACCAAATGGAA";
$seq =~ s/\s//g;
ESTScan::Compute("GCTACGAGCGTTA", -10, -10, -100, \@res, 
			$cdsIndex, $utrIndex, $startIndex, $stopIndex, 50);
ESTScan::Compute("GCTACGAGCGTTAGCTACGAGCGTTA", -10, -10, -100, \@res, 
			$cdsIndex, $utrIndex, $startIndex, $stopIndex, 60);

@res = ();
my $bigMax = ESTScan::Compute($seq, -150, -150, -100, \@res, 
		 $cdsIndex, $utrIndex, $startIndex, $stopIndex, 50);
if (length($seq) == 684 && $bigMax == 570 && $#res == 0) {
  print "ok 19\n";
} else {
  print "not ok 19\n";
}
my $r = $res[0];
if ($$r[0] == 570 && $$r[1] == 0 && $$r[2] == 558) {
  print "ok 20\n";
} else {
  print "not ok 20\n";
}
$seq .= $seq;
@res = ();
$bigMax = ESTScan::Compute($seq, -150, -150, -100, \@res, 
		 $cdsIndex, $utrIndex, $startIndex, $stopIndex, 50);
if (length($seq) == 1368 && $bigMax == 995 && $#res == 0) {
  print "ok 21\n";
} else {
  print "not ok 21\n";
}
$r = $res[0];
if ($$r[0] == 995 && $$r[1] == 0 && $$r[2] == 1242) {
  print "ok 22\n";
} else {
  print "not ok 22\n";
}
$seq = "TTCCATTTGGTCCAAGGCTTGTTAGGATAGTTAAGAAAGCTGCCTATTGGCCGGAGGGAG
AGGCTTAGGCAGAAGCCTATTACTTTGCAAGGGGCCCTTCAGAAGTCGCTGGGCTCAGAA
GGCTCTTAGTCGTGCTTGAGAGTGAGCCTTTCGAAGAGATACTCGCCCAGCCCAGCCTCC
GGGCCACCCAGCCTGTGGAGGTTGGTCAGGTGGTCACCCATCTTCTTGATAAGCTTCACT
TCCTCATCTAGGAAGTGAGTCTCCAGGAAGTCACAGAGATGGGGGTCCGTGCGGGCAGAA
CCCAGGGCATGAAGATCCAAAAGGGCCTGGTTCAGCTTTTTCTCCAGGGCCATGGCAGCT
TTCATGGCGTCTGGGGTTTTACCCACTCATCTTCAGCTGGCTTCTTGATGTCCTGGAAGA
GAGCGCGGCCCCACGCTGGTTTTGCATCTTCAGGAGACGCTCGTAGCCCTCGCGCTTCTC
CTCGGCCAGTTCGCGGAAGAAGTGGCTCACGCCTTCCAGAGCCACATCATCGCGGTCGAA
ATAGAAGCCCAGAGAGAGGTAGGTGTAGGAGGCCTGCAGGTACAAATTGACCAGGCTGTT
GACGGCTGCCTCCACGTCGGTGAAATAATTCTGACGAATCTGGGAGCTCATGGCTGGTTT
GGCAAGAGGAGCTAACCACCAAAA";
$seq =~ s/\s//g;
@res = ();
$bigMax = ESTScan::Compute($seq, -150, -150, -100, \@res, 
		   	   $cdsIndex, $utrIndex, $startIndex, $stopIndex, 50);
if (length($seq) == 684 && $bigMax == 268 && $#res == 0) {
  print "ok 23\n";
} else {
  print "not ok 23\n";
}
$r = $res[0];
if ($$r[0] == 268 && $$r[1] == 351 && $$r[2] == 674) {
  print "ok 24\n";
} else {
  print "not ok 24\n";
}
$seq = "AGTTGTTGCTTATGATGTGTGAGTGAACATATGCCATGCCTGGCCTTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACCA
TGAGCTCCCAGATTCGTCAGAATTATTCCACCGACGTGGAGGCAGCCGTCAACAGCCTGGTCAATTTGTACCTGCAGGCC
TCCTACACCTACCTCTCTCTGGGCTTCTATTTCGACCGCGATGATGTGGCTCTGGAAGGCGTGAGCCACTTCTTCCGCGA
ACTGGCCGAGGAGAAGCGCGAGGGCTACGAGCGTCTCCTGAAGATGCAAAACCAGCGTGGCGGCCGCGCTCTCTTCCAGG
ACATCAAGAAGCCAGCTGAAGATGAGTGGGGTAAAACCCCAGACGCCATGAAAGCTGCCATGACCCTGGAGAAAAAGCT
AACCAGGCCCTTTTGGATCTTCATGCCCTGGGTTCTGCCCGCACGGACCCCCATCTCTGTGACTTCCTGGAGACTCACTT
CCTAGATGAGGAAGTGAAGCTTATCAAGAAGATGGGTGACCACCTGACCAACCTCCACAGGCTGGGTGGCCCGGAGGCTG
GGCTGGGCGAGTATCTCTTCGAAAGGCTCACTCTCAAGCACGACTAAGAGCCTTCTGAGCCCAGCGACTTCTGAAGGGCC
CCTTGCAAAGTAATAGGGCTTCTGCCTAAGCCTCTCCCTCCAGCCAATAGGCAGCTTTCTTAACTATCCTAACAAGCCTT
GGA";
$seq =~ s/\s//g;
@res = ();
$bigMax = ESTScan::Compute($seq, -150, -150, -100, \@res, 
			   $cdsIndex, $utrIndex, $startIndex, $stopIndex, 50);
if (length($seq) == 722 && $bigMax == 594 && $#res == 0) {
  print "ok 25\n";
} else {
  print "not ok 25\n";
}
$r = $res[0];
if ($$r[0] == 594 && $$r[1] == 79 && $$r[2] == 605) {
  print "ok 26\n";
} else {
  print "not ok 26\n";
}
exit 0;
