# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..13\n"; }
END {print "not ok 1\n" unless $loaded;}
use ESTScan1;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

my @mat = ESTScan1::LoadMatrix("HumanIso.smat", "CODING", 4.0);
print $#mat == 1 ? "ok 2\n" : "not ok 2\n";
for my $i (0 .. $#mat) {
    my $order = $ {$mat[$i]}[3];
    my $frames = $ {$mat[$i]}[4];
    my $offset = $ {$mat[$i]}[5];
    my $aRef = $ {$mat[$i]}[8];
    my $ref = ESTScan1::CreateMatrix($order, $frames, $offset, $aRef, -20, 1);
    $ {$mat[$i]}[9] = $ref;
    if ($ref == $i && $order == 6 && $frames == 3 && $offset == 1) {
      print "ok ", $i + 3, "\n";
    } else {
      print "not ok ", $i + 3, "\n";
    }
}
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
ESTScan1::Compute("GCTACGAGCGTTA", -10, -10, 20, 40, \@res, 0, -1, 0, -1, 0, "");
ESTScan1::Compute("GCTACGAGCGTTAGCTACGAGCGTTA", -10, -10, 20, 40, \@res, 0, -1, 0, -1, 0, "");
@res = ();
my $bigMax = ESTScan1::Compute($seq, -150, -150, 100, 400, \@res, 0, -1, 0, -1, 0, "");
if (length($seq) == 684 && $bigMax == 644 && $#res == 0) {
  print "ok 5\n";
} else {
  print "not ok 5\n";
}
my $r = $res[0];
if ($$r[0] == 644 && $$r[1] == 15 && $$r[2] == 582) {
  print "ok 6\n";
} else {
  print "not ok 6\n";
}
$seq .= $seq;
@res = ();
$bigMax = ESTScan1::Compute($seq, -150, -150, 100, 200, \@res, 0, -1, 0, -1, 0, "");
if (length($seq) == 1368 && $bigMax == 648 && $#res == 1) {
  print "ok 7\n";
} else {
  print "not ok 7\n";
}
$r = $res[0];
if ($$r[0] == 644 && $$r[1] == 15 && $$r[2] == 582) {
  print "ok 8\n";
} else {
  print "not ok 8\n";
}
$r = $res[1];
if ($$r[0] == 648 && $$r[1] == 660 && $$r[2] == 1545) {
  print "ok 9\n";
} else {
  print "not ok 9\n";
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
$bigMax = ESTScan1::Compute($seq, -150, -150, 100, 400, \@res, 0, -1, 0, -1, 0, "");
if (length($seq) == 684 && $bigMax == 229 && $#res == 0) {
  print "ok 10\n";
} else {
  print "not ok 10\n";
}
$r = $res[0];
if ($$r[0] == 229 && $$r[1] == 405 && $$r[2] == 671) {
  print "ok 11\n";
} else {
  print "not ok 11\n";
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
$bigMax = ESTScan1::Compute($seq, -150, -150, 100, 400, \@res, 0, -1, 0, -1, 0, "");
if (length($seq) == 722 && $bigMax == 699 && $#res == 0) {
  print "ok 12\n";
} else {
  print "not ok 12\n";
}
$r = $res[0];
if ($$r[0] == 699 && $$r[1] == 70 && $$r[2] == 629) {
  print "ok 13\n";
} else {
  print "not ok 13\n";
}
exit 0;
