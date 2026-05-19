#!/usr/bin/perl -w
use strict;

open FILE_IN,"input/EPCrisprBenchmark_ensemble_data_GRCh38.K562_AllFeatures_NAfilled_withPreds.tsv";
my $string = <FILE_IN>;
chomp $string;
my @Header = split "\t", $string;
my @Pairs  = ();

while($string = <FILE_IN>)
{
    chomp $string;
    my @Data = split "\t", $string;
    my $pair = $Data[4];
    my $regulated = $Data[16];
    my $baseline  = sprintf("%.7f",$Data[69]);
    my $baseline_plus = sprintf("%.7f",$Data[70]);
    my $ABC = sprintf("%.7f",$Data[35]);
    my $ABC_plus = sprintf("%.7f",$Data[71]);
    my $Encode_E2G = sprintf("%.7f",$Data[67]);
    my $Encode_E2G_minus = sprintf("%.7f",$Data[68]);
    my $glsCoefficient = sprintf("%.7f",$Data[22]);
    my $pearsonCorrelation = sprintf("%.7f",$Data[72]);
    my $RNAglsCoefficient = sprintf("%.7f",$Data[73]);
    push @Pairs,[$pair,$regulated,$baseline,$baseline_plus,$ABC,$ABC_plus,$Encode_E2G,$Encode_E2G_minus,$glsCoefficient,$pearsonCorrelation,$RNAglsCoefficient];
}
close FILE_IN;

&PR_curve(2, "Baseline", "Baseline model");
&PR_curve(3, "Baseline_plus", "Baseline with glsCoefficient");
&PR_curve(4, "ABC", "ABC model");
&PR_curve(5, "ABC_plus", "ABC with glsCoefficient");
&PR_curve(6, "Encode_E2G", "Encode_E2G model");
&PR_curve(7, "Encode_E2G_minus", "Encode_E2G without glsCoefficient");
&PR_curve(8, "glsCoefficient", "glsCoefficient");
&PR_curve(9, "pearsonCorrelation", "pearsonCorrelation");
&PR_curve(10, "RNAglsCoefficient", "RNAglsCoefficient");

sub PR_curve
{
    my ($model,$label,$report) = @_;
    my (%Scores,@AAA) = ();
    @Pairs = sort { $b->[$model] <=> $a->[$model] } @Pairs;

    map { $Scores{$_->[$model]}++ } grep { $_->[1] eq "True" } @Pairs;
    my $highest_pair = $Pairs[0];
    $Scores{$highest_pair->[$model]}++;

    open FILE_OUT,">work/PRC.$label.data";
    print FILE_OUT "threshold\tprecision\trecall\tgroup\tTP\tFN\tFP\n";
    my @List = sort { $b <=> $a } keys %Scores;shift @List;

    foreach my $threshold (@List)
    {
	@AAA = grep { $_->[$model] >= $threshold and $_->[1] eq "True" }  @Pairs;
	my $true_positive  = 1 + $#AAA;

	@AAA = grep { $_->[$model]  < $threshold and $_->[1] eq "True" }  @Pairs;
	my $false_negative = 1 + $#AAA;

	@AAA = grep { $_->[$model] >= $threshold and $_->[1] eq "False" } @Pairs;
	my $false_positive = 1 + $#AAA;

	my $current_precision = sprintf("%.7f",$true_positive/($true_positive+$false_positive));
	my $current_recall    = sprintf("%.7f",$true_positive/($true_positive+$false_negative));

	print FILE_OUT "$threshold\t$current_precision\t$current_recall\t$report\t$true_positive\t$false_negative\t$false_positive\n";
    }
    close FILE_IN;
    close FILE_OUT;
    print "Done $report\n";
}
