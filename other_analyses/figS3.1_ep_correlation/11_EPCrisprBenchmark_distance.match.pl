#!/usr/bin/perl -w
use strict;

my @Observations = ();
my $seed = 1234;

open FILE_IN,"work/EPCrisprBenchmark_Pos_glsCoefficient.data";
my $string = <FILE_IN>;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @Data,"positive";
    push @Observations,\@Data;
}
close FILE_IN;

open FILE_IN,"work/EPCrisprBenchmark_Neg_glsCoefficient.data";
$string = <FILE_IN>;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @Data,"negative";
    push @Observations,\@Data;
}
close FILE_IN;

@Observations = sort { $a->[1] <=> $b->[1] } @Observations;

my @Positives = grep { $Observations[$_]->[3] eq "positive" } (0..$#Observations);
srand($seed);
my %Random = map { $_ , rand(1) } @Positives;

@Positives = sort { $Random{$a} <=> $Random{$b} } @Positives;
my $max_match = 300;
my $match_now = 0;
foreach my $index (@Positives)
{
    next if $index == 0;
    next unless $Observations[$index  ]->[3] eq "positive";
    if($Observations[$index-1]->[3] eq "negative")
    {
	$Observations[$index-1]->[3] = "selected";
	$match_now++;
    }
    elsif($index > 2 and $Observations[$index-2]->[3] eq "negative")
    {
	$Observations[$index-2]->[3] = "selected";
	$match_now++;
    }
    elsif($index > 3 and $Observations[$index-3]->[3] eq "negative")
    {
	$Observations[$index-3]->[3] = "selected";
	$match_now++;
    }
    elsif($index > 4 and $Observations[$index-4]->[3] eq "negative")
    {
	$Observations[$index-4]->[3] = "selected";
	$match_now++;
    }
    elsif($index > 5 and $Observations[$index-5]->[3] eq "negative")
    {
	$Observations[$index-5]->[3] = "selected";
	$match_now++;
    }
    last if $match_now >= $max_match;
}

open FILE_OUT,">work/EPCrisprBenchmark_Match_glsCoefficient.data";
print FILE_OUT "pair\tdistance\tglsCoefficient\n";
map { print FILE_OUT join "\t", @{$_}[0..2];print FILE_OUT "\n" } grep { $_->[3] eq "selected" } @Observations;
close FILE_OUT;

my ($distance,$total,$value) = ();
map { $distance+=$_->[1];$value+=$_->[2];$total++ } grep { $_->[3] eq "positive" } @Observations;
print "$total  positive have ",sprintf("average distance =  %.0f average glsCoefficient = %.4f",$distance/$total,$value/$total),"\n";

($distance,$total,$value) = ();
map { $distance+=$_->[1];$value+=$_->[2];$total++ } grep { $_->[3] ne "positive" } @Observations;
print "$total negative have ",sprintf("average distance = %.0f average glsCoefficient = %.4f",$distance/$total,$value/$total),"\n";

($distance,$total,$value) = ();
map { $distance+=$_->[1];$value+=$_->[2];$total++ } grep { $_->[3] eq "selected" } @Observations;
print "$total  selected have ",sprintf("average distance =  %.0f average glsCoefficient = %.4f",$distance/$total,$value/$total),"\n";
print "Seed = $seed\n";
