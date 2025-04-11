# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..1\n"; }
END {print "not ok 1\n" unless $loaded;}
use ESTScan;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

my @mat = ESTScan::LoadMatrix("/usr/molbio/share/ESTScan/HumanIso.smat", "CODING", 4.0);
for my $i (0 .. $#mat) {
    my $order = $ {$mat[$i]}[3];
    my $frames = $ {$mat[$i]}[4];
    my $offset = $ {$mat[$i]}[5];
    print("$i, Read matrix, $order, $frames, $offset\n");
    my $aRef = $ {$mat[$i]}[8];
    my $ref = ESTScan::CreateMatrix($order, $frames, $offset, $aRef, -20, 1);
    print "CreateMatrix returned $ref\n";
    $ {$mat[$i]}[9] = $ref;
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
ESTScan::Compute("GCTACGAGCGTTA", -10, -10, 20, 40, \@res, 0, -1, 0, -1, 0, "");
ESTScan::Compute("GCTACGAGCGTTAGCTACGAGCGTTA", -10, -10, 20, 40, \@res, 0, -1, 0, -1, 0, "");
@res = ();
my $bigMax = ESTScan::Compute($seq, -150, -150, 100, 400, \@res, 0, -1, 0, -1, 0, "");
print "Sequence length: ", length($seq), ", bigMax: $bigMax\n";
print "We have ", $#res + 1, " results.\n";
foreach my $r (@res) {
    print "Maximum: $$r[0] at position $$r[2], starting from $$r[1].\n";
    my $s = $$r[3];
    print ">fake\n";
    $s =~ s/(.{70})/$1\n/g;
    $s =~ s/\s+$//; # remove a trailing newline, since we add one below.
    print "$s\n";
}
$seq .= $seq;
@res = ();
$bigMax = ESTScan::Compute($seq, -150, -150, 100, 200, \@res, 0, -1, 0, -1, 0, "");
print "Sequence length: ", length($seq), ", bigMax: $bigMax\n";
print "We have ", $#res + 1, " results.\n";
foreach my $r (@res) {
    print "Maximum: $$r[0] at position $$r[2], starting from $$r[1].\n";
    my $s = $$r[3];
    print ">fake\n";
    $s =~ s/(.{70})/$1\n/g;
    $s =~ s/\s+$//; # remove a trailing newline, since we add one below.
    print "$s\n";
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
#ESTScan::Compute($seq, -150, -150, 100, 400, \@res, 0, -1, 0, -1, 0, "");
#print "Sequence length: ", length($seq), "\n";
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
ESTScan::Compute($seq, -150, -150, 100, 400, \@res, 0, -1, 0, -1, 0, "");
print "Sequence length: ", length($seq), "\n";
print "We have ", $#res + 1, " results.\n";
foreach my $r (@res) {
    print "Maximum: $$r[0] at position $$r[2], starting from $$r[1].\n";
    my $s = $$r[3];
    print ">fake\n";
    $s =~ s/(.{70})/$1\n/g;
    $s =~ s/\s+$//; # remove a trailing newline, since we add one below.
    print "$s\n";
}
