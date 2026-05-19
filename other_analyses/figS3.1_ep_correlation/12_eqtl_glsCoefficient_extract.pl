#!/usr/bin/perl -w
use strict;

print "It takes around 30 seconds\n";
open FILE_IN,"unpigz -p3 -c input/GM12878.activity.gz |";
my %Activity = ();
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    $Activity{$Data[0]} = $Data[1];
}
close FILE_IN;

open FILE_POS,"| pigz -p3 > work/GTEx_ge_LCL.eQTLs.positive.data.gz";
print FILE_POS "chrom\tsnp_pos\tpip\te_start\te_end\tgene_id\tsymbol\tdistance\tglsCoefficient\tactivity\tstatus\n";
open FILE_NEG,"| pigz -p3 > work/GTEx_ge_LCL.eQTLs.negative.data.gz";
print FILE_NEG "chrom\tsnp_pos\tpip\te_start\te_end\tgene_id\tsymbol\tdistance\tglsCoefficient\tactivity\tstatus\n";
my %Tested = ();
my ($positive,$negative) = ();
open FILE_IN,"bedtools intersect -a input/GTEx_ge_LCL.eQTLs.bed.gz -b input/correlation.DHS_tagAlign.100bp.tsv.gz -sorted -wao |";
while(my $string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    next if $Data[6] eq ".";
    my $region = "$Data[6]:$Data[7]-$Data[8]";
    next unless defined $Activity{$region};
    next if $Activity{$region} < 0.1;
    $Data[17] = $Activity{$region};
    next if defined $Tested{"$Data[0]\t$Data[1]\t$Data[9]"};
    $Tested{"$Data[0]\t$Data[1]\t$Data[9]"} = "yes";
    my $distance = abs($Data[11]-$Data[1]);
    $Data[12] = $distance;
    if($Data[3] eq $Data[9])
    {
	if($Data[5] >= 0.05)
	{
	    print FILE_POS join "\t", @Data[0,1,5,7,8,9,10,12,15,17];
	    print FILE_POS "\tpositive\n";
	    $positive++;
	}
    }
    else
    {
	$Data[5] = "set-to-zero";
	print FILE_NEG join "\t", @Data[0,1,5,7,8,9,10,12,15,17];
	print FILE_NEG "\tnegative\n";
	$negative++;
    }
}
close FILE_POS;
close FILE_NEG;
close FILE_IN;
print "Done bedtools of eQTLs and glsCoefficient\n";
print "positive = $positive\nnegative = $negative\n";
