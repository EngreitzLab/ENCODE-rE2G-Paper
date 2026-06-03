#!/usr/bin/perl -w
use strict;

my $Genes = load_classes("Genes");
my $Enhancers = load_classes("Enhancers");

open FILE_IN,"input/EPCrisprBenchmark_ensemble_data_GRCh38.K562_AllFeatures_withNA.FullModel.tsv";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
my %Cross_classes = ();

open FILE_OUT,">work/EPCrisprBenchmark.classification";
print FILE_OUT "dataset\tenhancer\tgene\tdistance\te_class\tg_class\tRegulated\tglsCoefficient\n";

while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next unless defined $Genes->{$Data[9]};
    next unless defined $Enhancers->{"$Data[1]:$Data[2]-$Data[3]"};
    next if $Data[22] eq "";
    $Data[24] =~ s/\.5$//;
    die $Data[24] if $Data[24] eq "";
    die $Data[24] if $Data[24] =~ /\D/;
    print FILE_OUT "$Data[0]\t$Data[1]:$Data[2]-$Data[3]\t$Data[9]\t$Data[24]\t",$Enhancers->{"$Data[1]:$Data[2]-$Data[3]"},"\t",$Genes->{$Data[9]},"\t$Data[16]\t$Data[22]\n";
}
close FILE_IN;


sub load_classes
{
    my $prefix = $_[0];
    open FILE_IN,"work/$prefix.classification";
    my $string = <FILE_IN>;
    my %ret_h = ();
    while(my $string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	$ret_h{$Data[0]} = $Data[1];
    }
    close FILE_IN;
    return \%ret_h;
}
