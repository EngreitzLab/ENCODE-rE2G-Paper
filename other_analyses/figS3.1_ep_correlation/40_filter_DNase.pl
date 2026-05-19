#!/usr/bin/perl -w
use strict;

open FILE_IN,"unpigz -p3 -c input/DNAse_seq.matrix.gz |";
my $string = <FILE_IN>;

open FILE_OUT,">work/DNase_seq.example.matrix";
print FILE_OUT $string;
while($string = <FILE_IN>)
{
    my @Data = split "\t", $string;
    print FILE_OUT $string if $Data[0] eq "HDAC6" or $Data[0] eq "CD69";
}
close FILE_IN;

open FILE_IN,"unpigz -p3 -c input/GRCh38-cCREs.V4.sorted.100bp.allSamples.counts.normalized.tsv.gz |";
$string = <FILE_IN>;
while($string = <FILE_IN>)
{
    my @Data = split "\t", $string;
    print FILE_OUT $string if $Data[0] eq "chrX:48800681-48801014" or $Data[0] eq "chr12:9764781-9765131";
}
close FILE_OUT;
close FILE_IN;
print "Done DNase filter\n";
