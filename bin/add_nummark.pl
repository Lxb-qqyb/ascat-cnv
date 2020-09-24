use strict;
my ($segfile,$ploidyfile,$baffile,$outdir) = @ARGV;
my %hash;
open (IN,"$segfile");

while (<IN>){
	chomp;
	next if ($_ =~ /sample/);
	$_ =~ s/\"//g;
	my @all = split /\t+|\s+/,$_;
	my $samplename = $all[1];
	my $chr = $all[2];
	my $start = $all[3];
	my $end = $all[4];
	push @{$hash{$samplename}{$chr}},[$start,$end,$all[5],$all[6],0];
}
close IN;
foreach my $sample (sort keys %hash){
	open (IN,"$baffile");
	<IN>;
	while (<IN>){
		chomp;
		next if $_ =~ /^\s/;
		my @all = split /\t+|\s+/,$_;
		my @info = @{$hash{$sample}{$all[1]}};
		foreach my $info (@info){
			my $start = $$info[0];
			my $end = $$info[1];
			my $num = $$info[4];
			if ($all[2] >= $start && $all[2] <= $end){
				$num++;
			}
			$$info[4] = $num;
			next;
		}
	}
	close IN;

}
foreach my $sample (sort keys %hash){
	open (OUT,">$outdir/$sample.seg.txt");
	print OUT "chr\tstart\tend\tnMajor\tnMinor\tnum.mark\n";
	foreach my $chr (sort_chrs (keys %{$hash{$sample}})){
		foreach my $info (@{$hash{$sample}{$chr}}){
			print OUT "$chr\t";
			print OUT join "\t",@$info;
			print OUT "\n";
		}
	}
	close OUT;
}
open (IN,"$ploidyfile");
open (OUT,">$outdir/ploidy.tsv");
while (<IN>){
	chomp;
	next if ($_ =~ /^\"x\"/);
	$_ =~ s/\"//g;
	print OUT "$_\n";
}
close OUT;
system ("mv $outdir/ploidy.tsv $outdir/ploidy.txt");

sub sort_chrs {
	my @unsort = @_;
	my @presort1;my @presort2;
	foreach my $chr (@unsort){
                if($chr =~ /^\d/ || $chr =~ /^chr\d/){
                        push @presort1,[$chr];
                }else{
                        push @presort2,[$chr];
                }
        }
        my @sort1 = map{$_}sort{$a->[0] <=> $b->[0]}@presort1;
        my @sort2 = map{$_}sort{$a->[0] cmp $b->[0]}@presort2;
        my @finalout;
        foreach my $sort1 (@sort1){
                my $sortchr = $$sort1[0];
                push @finalout,"$sortchr";
        }
        foreach my $sort2 (@sort2){
                my $sortchr = $$sort2[0];
                push @finalout,"$sortchr";
        }
        return @finalout;
}

