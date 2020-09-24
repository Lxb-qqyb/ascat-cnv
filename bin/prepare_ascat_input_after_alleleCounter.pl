#!/usr/bin/perl
use strict;
use List::Util qw(max min sum maxstr minstr shuffle);
my ($tumorcount,$normalcount,$outdir) = @ARGV;
my @tumorcount = split /,/,$tumorcount;
my %tumorinfo;my %normalinfo;
open (LOG,">$outdir/Sample.info.txt");
my $i = 1;
foreach my $tumorfile (@tumorcount){
	my $samplename = "S$i";
	my @samplename = split /\//,$tumorfile;
        my $filename = pop @samplename;
        my ($sample) = split /\./,$filename;
	print LOG "Tumor\t$sample\tS$i\n";
	open (IN,"$tumorfile");
	while (<IN>){
		chomp;
		next if ($_ =~ /^#CHR/);
		my @all = split /\t/,$_;
		my $site = "$all[0]\t$all[1]";
		my @atcg = @all[2,3,4,5];
		my @atcg_total = @all[2,3,4,5,6];
		my $check2 = check_cov (@atcg_total);
		if($check2){
			my $maxcount = max (@atcg);
			if($maxcount <= $all[6]){
				my $baf;
				if($all[6] > 0){
					$baf = $maxcount/$all[6];
				}else{
					$baf = 0;
				}
				$tumorinfo{$site}{$samplename} = "$baf\t$all[6]";
			}
		}
	}
	close IN;
	$i++;
}
print "Read tumor alleleCounter files finished!\n";
my @normalcount = split /,/,$normalcount;
$i = 1;
foreach my $normalfile (@normalcount){
        my $samplename = "S$i";
	my @samplename = split /\//,$normalfile;
        my $filename = pop @samplename;
        my ($sample) = split /\./,$filename;
	print LOG "Normal\t$sample\tS$i\n";
        open (IN,"$normalfile");
        while (<IN>){
                chomp;
                next if ($_ =~ /^#CHR/);
                my @all = split /\t+/,$_;
		next if ($all[6] < 10);
                my $site = "$all[0]\t$all[1]";
                my @atcg = @all[2,3,4,5];
                my @atcg_total = @all[2,3,4,5,6];
		my $check1 = check_hem (@atcg);
                my $check2 = check_cov (@atcg_total);
                if($check1 && $check2){
                        my $maxcount = max (@atcg);
                        if($maxcount <= $all[6]){
                                my $baf = $maxcount/$all[6];
                                $normalinfo{$site}{$samplename} = "$baf\t$all[6]";
                        }
                }
        }
        close IN;
        $i++;
}
close LOG;
print "Read normal alleleCounter files finished!\n";

my %sites;
foreach my $site (keys %normalinfo){
	my @normalinfotmp = (keys %{$normalinfo{$site}});
	next if (@normalinfotmp < @normalcount);
	my @tumorinfotmp = (keys %{$tumorinfo{$site}});
	next if (@tumorinfotmp < @tumorcount);
	$sites{$site} = 1;
}
print "Select SNP sites finished!\n";

open (OUT1,">$outdir/Tumor_BAF.tmp.txt");
open (OUT2,">$outdir/Tumor_LogR.tmp.txt");
open (OUT3,">$outdir/Germline_BAF.tmp.txt");
open (OUT4,">$outdir/Germline_LogR.tmp.txt");
print OUT1 "\tchrs\tpos";
print OUT2 "\tchrs\tpos";
print OUT3 "\tchrs\tpos";
print OUT4 "\tchrs\tpos";
my @siteskeys = keys %sites;
my $firstsite = shift @siteskeys;
foreach my $samplename (sort keys %{$normalinfo{$firstsite}}){
	print OUT2 "\t$samplename";
	print OUT4 "\t$samplename";
}
foreach my $samplename (sort keys %{$tumorinfo{$firstsite}}){
        print OUT1 "\t$samplename";
        print OUT3 "\t$samplename";
}
print OUT1 "\n";
print OUT2 "\n";
print OUT3 "\n";
print OUT4 "\n";
my $flag = 1;
my %tumorratio;
foreach my $site (sort_tab_key(keys %sites)){
	print OUT1 "SNP$flag\t$site";
	print OUT3 "SNP$flag\t$site";
	foreach my $samplename (sort keys %{$tumorinfo{$site}}){
		my @t_baf_depth = split /\t/,$tumorinfo{$site}{$samplename};
		print OUT1 "\t$t_baf_depth[0]";
		my $t_depth = $t_baf_depth[1];
		my @n_baf_depth = split /\t/,$normalinfo{$site}{$samplename};
		print OUT3 "\t$n_baf_depth[0]";
		my $n_depth = $n_baf_depth[1];
		my $tumorratio = $t_depth/$n_depth;  #$n_depth never be 0;
		push @{$tumorratio{$samplename}},$tumorratio;
	}
	$flag++;
	print OUT1 "\n";
	print OUT3 "\n";
}
close OUT1;
close OUT3;
print "Calcuate tumor SNP's logr and write tumor and normal's baf finished!\n";

my %ratio_ref;
foreach my $sample (keys %tumorratio){
	my @ratiotmp = @{$tumorratio{$sample}};
	my $t_mean_ratio = average (@ratiotmp);
	$ratio_ref{$sample} = $t_mean_ratio;
}
print "Calcuate tumor SNP's mean logr finished!\n";

$flag = 1;
foreach my $site (sort_tab_key(keys %sites)){
        print OUT2 "SNP$flag\t$site";
        print OUT4 "SNP$flag\t$site";
        foreach my $samplename (sort keys %{$tumorinfo{$site}}){
                my @t_baf_depth = split /\t/,$tumorinfo{$site}{$samplename};
                my $t_depth = $t_baf_depth[1];
                my @n_baf_depth = split /\t/,$normalinfo{$site}{$samplename};
                my $n_depth = $n_baf_depth[1];
                my $tumorratio = $t_depth/$n_depth;  #$n_depth never be 0;
                my $log2ratio = log2($tumorratio/$ratio_ref{$samplename});
		print OUT2 "\t$log2ratio";
		print OUT4 "\t0";
        }
	$flag++;
	print OUT2 "\n";
	print OUT4 "\n";
}
close OUT2;
close OUT4;
print "Write tumor and normal's final logr finished!\n";

my %info;
open (IN,"$outdir/Sample.info.txt");
while (<IN>){
	chomp;
	my ($type,$name,$Sinfo) = split /\t/,$_;
	$info{$Sinfo}{$type} = $name;
}
close IN;

open (IN,"$outdir/Tumor_BAF.tmp.txt");
open (OUTS,">$outdir/Tumor_BAF.txt");
while (<IN>){
	chomp;
	if ($. == 1){
		my @line = split /\t+|\s+/,$_;
		foreach my $i (0..$#line){
			if(exists $info{$line[$i]}){
				$line[$i] = $info{$line[$i]}{Tumor};
			}
		}
		print OUTS join "\t",@line;
		print OUTS "\n";
		
	}else{
		print OUTS "$_\n";
	}
}
close OUTS;
close IN;

open (IN,"$outdir/Tumor_LogR.tmp.txt");
open (OUTS,">$outdir/Tumor_LogR.txt");
while (<IN>){
        chomp;
        if ($. == 1){
                my @line = split /\t+|\s+/,$_;
                foreach my $i (0..$#line){
                        if(exists $info{$line[$i]}){
                                $line[$i] = $info{$line[$i]}{Tumor};
                        }
                }
                print OUTS join "\t",@line;
                print OUTS "\n";

        }else{
                print OUTS "$_\n";
        }
}
close OUTS;
close IN;

open (IN,"$outdir/Germline_BAF.tmp.txt");
open (OUTS,">$outdir/Germline_BAF.txt");
while (<IN>){
        chomp;
        if ($. == 1){
                my @line = split /\t+|\s+/,$_;
                foreach my $i (0..$#line){
                        if(exists $info{$line[$i]}){
                                $line[$i] = $info{$line[$i]}{Normal};
                        }
                }
                print OUTS join "\t",@line;
                print OUTS "\n";

        }else{
                print OUTS "$_\n";
        }
}
close OUTS;
close IN;

open (IN,"$outdir/Germline_LogR.tmp.txt");
open (OUTS,">$outdir/Germline_LogR.txt");
while (<IN>){
        chomp;
        if ($. == 1){
                my @line = split /\t+|\s+/,$_;
                foreach my $i (0..$#line){
                        if(exists $info{$line[$i]}){
                                $line[$i] = $info{$line[$i]}{Normal};
                        }
                }
                print OUTS join "\t",@line;
                print OUTS "\n";

        }else{
                print OUTS "$_\n";
        }
}
close OUTS;
system ("rm $outdir/Germline_LogR.tmp.txt $outdir/Germline_BAF.tmp.txt $outdir/Tumor_BAF.tmp.txt $outdir/Tumor_LogR.tmp.txt $outdir/Sample.info.txt");

sub log2 {
	my $n = shift;
	return log($n)/log(2);
}

sub sort_tab_key {
	my @unsort = @_;
	my @presort1;my @presort2;
	foreach my $unsort (@unsort){
		my ($chr,$pos) = split /\t/,$unsort;
		if($chr =~ /^\d/ || $chr =~ /^chr\d/){
			push @presort1,[$chr,$pos];
		}else{
			push @presort2,[$chr,$pos];
		}
	}
	my @sort1 = map{$_}sort{$a->[0] <=> $b->[0] || $a->[1]<=>$b->[1]}@presort1;
	my @sort2 = map{$_}sort{$a->[0] cmp $b->[0] || $a->[1]<=>$b->[1]}@presort2;
	my @finalout;
	foreach my $sort1 (@sort1){
		my $sortchr = $$sort1[0];
		my $sortpos = $$sort1[1];
		push @finalout,"$sortchr\t$sortpos";
	}
	foreach my $sort2 (@sort2){
                my $sortchr = $$sort2[0];
                my $sortpos = $$sort2[1];
                push @finalout,"$sortchr\t$sortpos";
        }
	return @finalout;
}

sub check_hem {
	my (@basecount) = @_;
	my $zerocount = grep { $_ eq "0" } @basecount;
	if($zerocount > 2){
		return 0;
	}else{
		return 1;
	}
} 
sub check_cov {
	my (@basecount) = @_;
        my $zerocount = grep { $_ eq "0" } @basecount;
        if($zerocount >= 4){
                return 0;
        }else{
                return 1;
        }
}

sub average {
	my (@num) = @_;
	my $num = scalar @num;
	my $total;
	foreach (0..$#num) {
		$total += $num[$_];
	}
	return ($total/$num);
}
