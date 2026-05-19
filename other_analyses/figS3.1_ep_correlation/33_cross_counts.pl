#!/usr/bin/perl -w
use strict;

open FILE_IN,"work/EPCrisprBenchmark.classification";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
my %Cross_classes = ();
my $total = 0;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    $Cross_classes{"$Data[4]<>$Data[5]"}++;
    $total++;
    if($Data[6] eq "TRUE")
    {
        $Cross_classes{"$Data[4]<>$Data[5]<>TRUE"}++;
    }
}
close FILE_IN;
map { print $_,"\t",$Cross_classes{$_},"\t",sprintf("(%.1f%s)",100*$Cross_classes{$_}/$total,"%"),"\n"; } grep { $_ !~ /TRUE/ } sort keys %Cross_classes;
print "\n";
foreach my $enh ("broad","specific","zero")
{
    print "enhancer $enh active\t",$Cross_classes{"$enh<>broad"}," (",$Cross_classes{"$enh<>broad<>TRUE"},")\t",$Cross_classes{"$enh<>variable"}," (",$Cross_classes{"$enh<>variable<>TRUE"},")\n";
}

open FILE_IN,"work/EPCrisprBenchmark.classification";
$string = <FILE_IN>;

open FILE_ONE,">work/Classes.broad_broad.data";
print FILE_ONE $string;

open FILE_TWO,">work/Classes.variable_broad.data";
print FILE_TWO $string;

open FILE_THREE,">work/Classes.broad_specific.data";
print FILE_THREE $string;

open FILE_FOUR,">work/Classes.variable_specific.data";
print FILE_FOUR $string;

while(my $string = <FILE_IN>)
{
    my @Data = split "\t", $string;
    next if $Data[4] eq "zero";
    if("$Data[4]_$Data[5]" eq "broad_broad")
    {
        print FILE_ONE $string;
    }
    elsif("$Data[4]_$Data[5]" eq "broad_variable")
    {
        print FILE_TWO $string;
    }
    elsif("$Data[4]_$Data[5]" eq "specific_broad")
    {
        print FILE_THREE $string;
    }
    elsif("$Data[4]_$Data[5]" eq "specific_variable")
    {
        print FILE_FOUR $string;
    }
    else
    {
        die "$Data[4]_$Data[5]";
    }
}
