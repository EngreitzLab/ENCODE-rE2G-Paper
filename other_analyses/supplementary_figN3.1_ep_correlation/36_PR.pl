#!/usr/bin/perl -w
use strict;

map { PR_curve($_) } ("broad_broad","broad_specific","variable_broad","variable_specific");

sub PR_curve
{
    my $label = $_[0];
    open FILE_IN,"work/Classes.$label.data";
    my $string = <FILE_IN>;
    chomp $string;
    my @Header = split "\t", $string;
    my @Pairs  = ();
    while($string = <FILE_IN>)
    {
	chomp $string;
	my @Data = split "\t", $string;
	my $pair = $Data[1];
	my $regulated = $Data[6];
	my $score = sprintf("%.7f",$Data[7]);
	push @Pairs,["$pair","$regulated","$score"];
    }
    close FILE_IN;

    my $model = 2;
    my (%Scores,@AAA) = ();
    @Pairs = sort { $b->[$model] <=> $a->[$model] } @Pairs;

    map { $Scores{$_->[$model]}++ } grep { $_->[1] eq "TRUE" } @Pairs;
    my $highest_pair = $Pairs[0];
    $Scores{$highest_pair->[$model]}++;

    open FILE_OUT,">work/Classes.$label.PRC.data";
    print FILE_OUT "threshold\tprecision\trecall\tgroup\tTP\tFN\tFP\n";
    my @List = sort { $b <=> $a } keys %Scores;shift @List;

    foreach my $threshold (@List)
    {
	@AAA = grep { $_->[$model] >= $threshold and $_->[1] eq "TRUE" }  @Pairs;
	my $true_positive  = 1 + $#AAA;

	@AAA = grep { $_->[$model]  < $threshold and $_->[1] eq "TRUE" }  @Pairs;
	my $false_negative = 1 + $#AAA;

	@AAA = grep { $_->[$model] >= $threshold and $_->[1] ne "TRUE" } @Pairs;
	my $false_positive = 1 + $#AAA;

	my $current_precision = sprintf("%.7f",$true_positive/($true_positive+$false_positive));
	my $current_recall    = sprintf("%.7f",$true_positive/($true_positive+$false_negative));

	print FILE_OUT "$threshold\t$current_precision\t$current_recall\t$label\t$true_positive\t$false_negative\t$false_positive\n";
    }
    close FILE_IN;
    close FILE_OUT;
    print "Done $label\n";
}

__END__
