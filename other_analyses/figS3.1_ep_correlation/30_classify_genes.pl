#!/usr/bin/perl -w
use strict;

open FILE_IN,"input/EPCrisprBenchmark_ensemble_data_GRCh38.K562_AllFeatures_withNA.FullModel.tsv";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
my %Symbols = ();

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    $Symbols{$Data[9]}++;
}
close FILE_IN;

open FILE_OUT,">work/Genes.classification";
print FILE_OUT "symbol\tclass\n";
my ($broadly,$varably) = (0,0);
open FILE_IN,"unpigz -p3 -c input/DNAse_seq.matrix.gz |";
$string = <FILE_IN>;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $gene_id = shift @Data;
    next unless defined $Symbols{$gene_id};
    my @Above = grep { $_ > 2 } @Data;
    if($#Above == 88)
    {
	print FILE_OUT "$gene_id\tbroad\n";
	$broadly++;
    }
    else
    {
	print FILE_OUT "$gene_id\tvariable\n";
	$varably++;
    }
}
close FILE_IN;
close FILE_OUT;
print "b = $broadly\nv = $varably\n";
