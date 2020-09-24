#!/usr/bin/perl -w
use lib '/GPFS01/home/zhouyh/perl5/lib/perl5/';
use strict;
use Cwd;
use Cwd 'abs_path';
use DateTime;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin '$Bin';
sub usage
{
        print STDERR <<USAGE;
==============================================================================
Description     
Options
	-dir <s>: dir of bam files
	-outdir <s>: pathway of outdir
	-info <s>: pair information of tumor-normal, for example "normal tumor patiant",first cloumn is normal,second column is tumor,third column is patient
	-alleleCounter <s>: path of alleleCounter
	-multi: run ASCAT with multi-sample segmentation (when shared breakpoints are expected)
	-cpu <i>: cpu used
	-node <s>: nodes used
        -h|?|help : Show this help
==============================================================================
USAGE
}
my ($help,$dir,$outdir,$info,$alleleCounter,$cpu,$nodes,$multi);
GetOptions(
        "h|?|help"=>\$help,
        "outdir=s"=>\$outdir,
        "dir=s"=>\$dir,
	"info=s"=>\$info,
	"alleleCounter=s"=>\$alleleCounter,
	"multi!"=>\$multi,
	"cpu=i"=>\$cpu,
	"node=s"=>\$nodes,
);
if(!defined ($info) || !defined($dir) || defined($help)){
        &usage;
        exit 0;
}

$cpu ||= 5;
$nodes ||= "ibm";
$outdir ||= getcwd();
$alleleCounter ||= "/GPFS01/home/liuxb/anaconda3/envs/cgp2/bin/alleleCounter";
my $snpbed = "$Bin/bed/Select.1000g.snp.bed";
my $ab_outdir = abs_path ($outdir);
my $ab_dir = abs_path ($dir);
open (IN,"$info");
system ("mkdir -p $ab_outdir/script") unless (-e "$ab_outdir/script");
system ("mkdir -p $ab_outdir/alleleCounter") unless (-e "$ab_outdir/alleleCounter");
open (OUT,">$ab_outdir/script/run_alleleCounter.sh");
my %hash;my %flag;
while (<IN>){
	chomp;
	my ($tumor,$normal,$patient) = split /\t+|\s+/,$_;
	system ("mkdir -p $ab_outdir/alleleCounter/$patient") unless (-e "$ab_outdir/alleleCounter/$patient");
	push @{$hash{$patient}},[$tumor,$normal];
	if(!exists $flag{$patient}{$tumor}){
		print OUT "$alleleCounter -l $snpbed -b $ab_dir/$tumor.sorted.rmdup.realigned.recal.bam -o $ab_outdir/alleleCounter/$patient/$tumor.alleleCounter.txt\n";
		$flag{$patient}{$tumor} = 1;
	}
	if(!exists $flag{$patient}{$normal}){
		print OUT "$alleleCounter -l $snpbed -b $ab_dir/$normal.sorted.rmdup.realigned.recal.bam -o $ab_outdir/alleleCounter/$patient/$normal.alleleCounter.txt\n";
		$flag{$patient}{$normal} = 1;
	}
}
close IN;
close OUT;

chdir ("$ab_outdir/script");
#pid_stat("run_alleleCounter.sh");

open (OUT,">$ab_outdir/script/run_prepare_ascat.sh");
foreach my $patient (keys %hash){
	system ("mkdir -p $ab_outdir/Ascat_Result/$patient") unless (-e "$ab_outdir/Ascat_Result/$patient");
	my $tumorfile;my $normalfile;
	foreach my $sample (@{$hash{$patient}}){
		$tumorfile .= "$ab_outdir/alleleCounter/$patient/$$sample[0].alleleCounter.txt,";
		$normalfile .= "$ab_outdir/alleleCounter/$patient/$$sample[1].alleleCounter.txt,";
	}
	$tumorfile =~ s/,$//;
	$normalfile =~ s/,$//;
	if($multi){
		print OUT "perl $Bin/bin/prepare_ascat_input_after_alleleCounter.pl $tumorfile $normalfile $ab_outdir/Ascat_Result/$patient && Rscript $Bin/bin/asmultipcf.R $ab_outdir/Ascat_Result/$patient && perl $Bin/bin/add_nummark.pl $ab_outdir/Ascat_Result/$patient/segment.txt $ab_outdir/Ascat_Result/$patient/ploidy.txt $ab_outdir/Ascat_Result/$patient/Tumor_BAF.txt $ab_outdir/Ascat_Result/$patient\n";
	}else{
		print OUT "perl $Bin/bin/prepare_ascat_input_after_alleleCounter.pl $tumorfile $normalfile $ab_outdir/Ascat_Result/$patient && Rscript $Bin/bin/aspcf.R $ab_outdir/Ascat_Result/$patient && perl $Bin/bin/add_nummark.pl $ab_outdir/Ascat_Result/$patient/segment.txt $ab_outdir/Ascat_Result/$patient/ploidy.txt $ab_outdir/Ascat_Result/$patient/Tumor_BAF.txt $ab_outdir/Ascat_Result/$patient\n";
	}
}
close OUT;
pid_stat("run_prepare_ascat.sh");


sub pid_stat{
        my $qsub_cmd = shift;
        my @cmd_id;
        if (-e "lsf_$qsub_cmd"){
                unlink glob "lsf_$qsub_cmd/* lsf_$qsub_cmd/.*";
                rmdir "lsf_$qsub_cmd" or die "cannot rmdir lsf_$qsub_cmd";
        }
        system ("/GPFS02/gubh/script/bsubjobs.py -c $cpu -o $nodes $qsub_cmd");
        open CMD,'<',"./lsf_$qsub_cmd/sucesslist.csv" or die "Can't open the file lsf_$qsub_cmd/sucesslist.csv";
        while (<CMD>){
                if (/^(\d+),/){
                push @cmd_id, $1;
                }
        }
        my $a = 0 ;
        while (1){
               my $tmp1 = `bjobs`;
                my @all_id = $tmp1 =~ /^(\d+)/mg;
                my $i = 0;
                foreach my $id (@cmd_id) {
                        if (@all_id ~~ /$id$/){$i++;}
                 }
                if ($a != $i){
                        &log("$i job(s) running!");
                        $a = $i;
                }
                last  if ($i == 0);
                sleep 5;
        }
}
sub log(){
        my $message=shift;
        my $dt=DateTime->now()->set_time_zone('Asia/Chongqing');
        $message.="\t".$dt->ymd()." ".$dt->hms() ;
        print "$message\n";
}


