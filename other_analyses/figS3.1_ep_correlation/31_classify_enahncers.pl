#!/usr/bin/perl -w
use strict;

open FILE_IN,"input/EPCrisprBenchmark_ensemble_data_GRCh38.K562_AllFeatures_withNA.FullModel.tsv";
my $string = <FILE_IN>;
my %Regions = ();
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    $Regions{"$Data[1]\t$Data[2]\t$Data[3]\t$Data[1]:$Data[2]-$Data[3]"}++;
}
close FILE_IN;

open FILE_OUT,"| bedSort stdin work/EPCrisprBenchmark.regions.bed";
map { print FILE_OUT $_,"\n" } keys %Regions;
close FILE_OUT;

open FILE_IN,"bedtools intersect -a work/EPCrisprBenchmark.regions.bed -b input/DHS_tagAlign.bed.gz -wao -sorted | ";
my %Mapping = ();
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @{$Mapping{$Data[3]}},\@Data;
}
close FILE_IN;

open FILE_OUT,">work/Enhancers.mapping";
foreach my $enhancer (sort keys %Mapping)
{
    my @Sorted = sort { $b->[-1] <=> $a->[-1] } @{$Mapping{$enhancer}};
    next if $Sorted[0]->[4] eq ".";
    print FILE_OUT $enhancer,"\t",$Sorted[0]->[4],":",$Sorted[0]->[5],"-",$Sorted[0]->[6],"\t",$Sorted[0]->[7],"\n";
}
close FILE_OUT;
print "Done enhancer mapping\n";

%Mapping = ();
open FILE_IN,"work/Enhancers.mapping";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    if(defined $Mapping{$Data[1]})
    {
	print "CCRE $Data[1] overlaps enhancers: ",$Mapping{$Data[1]}," and $Data[0]\n";
	$Mapping{$Data[1]} .=",$Data[0]";
    }
    else
    {
	$Mapping{$Data[1]} = $Data[0];
    }
}
close FILE_IN;

open FILE_OUT,">work/Enhancers.classification";
print FILE_OUT "region\tclass\n";
open FILE_IN,"unpigz -p3 -c input/GRCh38-cCREs.V4.sorted.100bp.allSamples.counts.normalized.tsv.gz |";
my ($broad,$specific,$zero) = (0,0,0);
$string = <FILE_IN>;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $region = shift @Data;
    next unless defined $Mapping{$region};
    my @Above = grep { $_ > 1.5 } @Data;
    if($#Above == -1)
    {
	foreach my $enhancer (split ",", $Mapping{$region})
	{
	    print FILE_OUT "$enhancer\tzero\n";
	    $zero++;
	}
    }
    elsif($#Above >= 3)
    {
	foreach my $enhancer (split ",", $Mapping{$region})
	{
	    print FILE_OUT "$enhancer\tbroad\n";
	    $broad++;
	}
    }
    else
    {
	foreach my $enhancer (split ",", $Mapping{$region})
	{
	    print FILE_OUT "$enhancer\tspecific\n";
	    $specific++;
	}
    }
}
close FILE_OUT;
close FILE_IN;
print "b = $broad\ns = $specific\nz = $zero\n";
