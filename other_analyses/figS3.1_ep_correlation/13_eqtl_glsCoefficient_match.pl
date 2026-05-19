#!/usr/bin/perl -w
use strict;

my @Observations = ();
my $seed = 1234;

open FILE_IN,"unpigz -p3 -c work/GTEx_ge_LCL.eQTLs.positive.data.gz |";
my $string = <FILE_IN>;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @Observations,\@Data;
}
close FILE_IN;

open FILE_IN,"unpigz -p3 -c work/GTEx_ge_LCL.eQTLs.negative.data.gz |";
$string = <FILE_IN>;
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    push @Observations,\@Data;
}
close FILE_IN;

@Observations = sort { $a->[7] <=> $b->[7] } @Observations;

my ($distance,$total,$value) = (0,0,0);
map { $distance+=$_->[7];$value+=$_->[-3];$total++ } grep { $_->[-1] eq "positive" } @Observations;
print "$total    positive have ",sprintf("average distance =  %.0f average glsCoefficient = %.4f",$distance/$total,$value/$total),"\n";

($distance,$total,$value) = (0,0,0);
map { $distance+=$_->[7];$value+=$_->[-3];$total++ } grep { $_->[-1] ne "positive" } @Observations;
print "$total negative have ",sprintf("average distance = %.0f average glsCoefficient = %.4f",$distance/$total,$value/$total),"\n";

my @Match = ();
foreach my $index (0..$#Observations)
{
    next unless $Observations[$index]->[-1] eq "positive";
    my $nearest_negative = get_nearest($index);
    $Observations[$nearest_negative]->[-1] = "selected 0";
    my @Replace = map { "$_" } @{$Observations[$nearest_negative]};
    push @Match,\@Replace;
}

foreach my $index (0..$#Observations)
{
    next unless $Observations[$index]->[-1] eq "positive";
    my $nearest_negative = get_nearest($index);
    $Observations[$nearest_negative]->[-1] = "selected 1";
    my @Replace = map { "$_" } @{$Observations[$nearest_negative]};
    push @Match,\@Replace;
}

foreach my $index (0..$#Observations)
{
    next unless $Observations[$index]->[-1] eq "positive";
    my $nearest_negative = get_nearest($index);
    $Observations[$nearest_negative]->[-1] = "selected 3";
    my @Replace = map { "$_" } @{$Observations[$nearest_negative]};
    push @Match,\@Replace;
}

foreach my $index (0..$#Observations)
{
    next unless $Observations[$index]->[-1] eq "positive";
    my $nearest_negative = get_nearest($index);
    $Observations[$nearest_negative]->[-1] = "selected 4";
    my @Replace = map { "$_" } @{$Observations[$nearest_negative]};
    push @Match,\@Replace;
}

#### I filter @Match for very low distances that are smaller than smallest positive
@Match = grep { $_->[7] > 986 } @Match;

open MATCH,"| pigz -p3 > work/GTEx_ge_LCL.eQTLs.match.data.gz";
print MATCH "chrom\tsnp_pos\tpip\te_start\te_end\tgene_id\tsymbol\tdistance\tglsCoefficient\tactivity\tstatus\n";
map { print MATCH join "\t", @{$_};print MATCH "\n" } @Match;
close MATCH;

($distance,$total,$value) = (0,0,0);
map { $distance+=$_->[7];$value+=$_->[-3];$total++ } @Match;
print "$total   selected have ",sprintf("average distance =  %.0f average glsCoefficient = %.4f",$distance/$total,$value/$total),"\n";

sub get_nearest
{
    my $index = $_[0];
    my $before = "NA";
    my $local = "$index";
    while($local > 0)
    {
	if($Observations[$local]->[-1] eq "negative")
	{
	    $before = $local;
	    last;
	}
	$local--;
    }
    $local = "$index";
    my $after = "NA";
    while($local <= $#Observations)
    {
	if($Observations[$local]->[-1] eq "negative")
	{
	    $after = $local;
	    last;
	}
	$local++;
    }
    if($before ne "NA" and $after ne "NA")
    {
	return $before if abs($Observations[$before]->[7] - $Observations[$index]->[7]) <= abs($Observations[$after]->[7] - $Observations[$index]->[7]);
	return $after;
    }
    return $after if $before eq "NA";

    print "index = $index. before = $before. after = $after\n";
    print "before: ",(join "\t", @{$Observations[$before]}),"\n";
    print "curent: ",(join "\t", @{$Observations[$index]}),"\n";
    print "after:  ",(join "\t", @{$Observations[$after]}),"\n";
    exit;
}
