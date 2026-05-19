#!/usr/bin/perl -w
use strict;

open FILE_OUT,"| LC_ALL=C sort --parallel=3 --buffer-size=4G -k 1,1n - --temporary-directory='./temp/' | pigz -p4 > work/glsCoefficient.data.gz";

open FILE_IN,"unpigz -p3 -c input/correlation.DHS_tagAlign.100bp.tsv.gz |";
my $string = <FILE_IN>;
my $total = 0;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    print FILE_OUT "$Data[6]\t$Data[9]\n";
    $total++;
    print "Parced $total\n" if $total%500000 == 0;
}
close FILE_IN;
close FILE_OUT;
print "Done. Total = $total\n";
