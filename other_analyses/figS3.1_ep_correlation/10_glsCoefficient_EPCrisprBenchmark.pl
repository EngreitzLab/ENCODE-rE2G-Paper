#!/usr/bin/perl -w
use strict;

open FILE_IN,"input/EPCrisprBenchmark_ensemble_data_GRCh38.K562_AllFeatures_withNA.FullModel.tsv";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;

open FILE_POS,">work/EPCrisprBenchmark_Pos_glsCoefficient.data";
print FILE_POS "pair\tdistance\tglsCoefficient\n";

open FILE_NEG,">work/EPCrisprBenchmark_Neg_glsCoefficient.data";
print FILE_NEG "pair\tdistance\tglsCoefficient\n";

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[22] eq "";

    if($Data[16] eq "TRUE")
    {
	print FILE_POS "$Data[4]\t$Data[24]\t$Data[22]\n";
    }
    elsif($Data[16] eq "FALSE")
    {
	print FILE_NEG "$Data[4]\t$Data[24]\t$Data[22]\n";
    }
    else
    {
	map { print $_,"\t",$Header[$_],"\t",$Data[$_],"\n" } (0..$#Data);exit;
    }
}
close FILE_POS;
close FILE_NEG;
close FILE_IN;
print "Done EPCrisprBenchmark glsCoefficient\n";
