#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);
use FindBin qw($Bin);

die "perl <input list> <outputdir> <prefix>  " unless (@ARGV==3);

my $outdir= $ARGV[1];
my $pre= $ARGV[2];
$outdir ||= "./";
system "mkdir $outdir " unless(-d "$outdir");
system "mkdir $outdir/lis " unless(-d "$outdir/lis");

open IN, "$ARGV[0]" || die $!;
open OUT, ">$outdir/lis/$pre.2lis" || die $!;
open LOG, ">$outdir/lis/$pre.2log" || die $!;
open IOG, ">$outdir/lis/$pre.2idlog" || die $!;
my $cut=0.1;

my %mitocondon = (
        "TTT" => "F",   "TTC" => "F",   "TTA" => "L",
        "TTG" => "L",   "TCT" => "S",   "TCC" => "S",
        "TCA" => "S",   "TCG" => "S",   "TAT" => "Y",
        "TAC" => "Y",   "TGT" => "C",   "TGC" => "C",
        "TGG" => "W",   "CTT" => "L",   "CTC" => "L",
        "CTA" => "L",   "CTG" => "L",   "CCT" => "P",
        "CCC" => "P",   "CCA" => "P",   "CCG" => "P",
        "CAT" => "H",   "CAC" => "H",   "CAA" => "Q",
        "CAG" => "Q",   "CGT" => "R",   "CGC" => "R",
        "CGA" => "R",   "CGG" => "R",   "ATT" => "I",
        "ATC" => "I",   "ATA" => "M",   "ATG" => "M",
        "ACT" => "T",   "ACC" => "T",   "ACA" => "T",
        "ACG" => "T",   "AAT" => "N",   "AAC" => "N",
        "AAA" => "K",   "AAG" => "K",   "AGT" => "S",
        "AGC" => "S",   "AGA" => "S",   "AGG" => "S",
        "GTT" => "V",   "GTC" => "V",   "GTA" => "V",
        "GTG" => "V",   "GCT" => "A",   "GCC" => "A",
        "GCA" => "A",   "GCG" => "A",   "GAT" => "D",
        "GAC" => "D",   "GAA" => "E",   "GAG" => "E",
        "GGT" => "G",   "GGC" => "G",   "GGA" => "G",
        "GGG" => "G",	"TGA" => "W",
);

while(<IN>){
	my %hash;
	my %count;
	chomp;
	my $ori;
	my $nam=$_;
	my $na=$nam;
	$na=~s/^.*_(For)/$1\t/;
	$na=~s/^.*_(Rev)/$1\t/;
	if ($nam=~/For/){$ori=1}elsif($nam=~/Rev/){$ori=2}else{die "$nam \t file name cannot judge direction\n"}
	open FH1, "$nam" || die $!;
	while (my $t=<FH1>){
		chomp($t);
		chomp(my $seq=<FH1>);
		#seq number
		if (exists $hash{$seq}){ 
			$hash{$seq}.=$t;
			$count{$seq}++;	
		}else{
			$hash{$seq}=$t;
			$count{$seq}=1;
		}
	}
	close FH1;

	if ($ori==1){
		TNT:for my $key (keys %hash){
			my $len=length $key;
			#codon
			for (my $i=1;$i<=$len-3;$i+=3){
				my $con=substr ($key,$i,3);
				unless (exists $mitocondon{$con}){
					delete $hash{$key};
					next TNT;
				}
			}
		}
	}else{
		TNT:for my $key (keys %hash){
			my $len=length $key;  #hifi02 Rev length 119
			my $keyr = $key;
			$keyr = reverse $keyr;
			$keyr =~ tr/ATCG/TAGC/;
			my $sta=($len % 3);  #658bp -- 119 --- the 3rd
			for (my $i=$sta;$i<=$len-3;$i+=3){
				my $con=substr ($keyr,$i,3);
				unless (exists $mitocondon{$con}){
					delete $hash{$key};
					next TNT;
				}
			}
		}
		
	}
	##usearch delet identity>=98% 
	my $tout="temp.fa";
	my $teuc="temp.uc";
	open (IFA,">",$tout);
	print IOG "$na";
	for my $keyi  (sort {$count{$b} <=> $count{$a}} keys %hash){
		print IFA ">$keyi\n$keyi\n";
	}
	system "$Bin/vsearch --cluster_smallmem temp.fa --uc temp.uc --id 0.98";
	open (UC,"<",$teuc);
	while (<UC>){
		chomp;
		next unless ($_=~/^H/);
		my ($ide,$ki)=(split (/\t/,$_))[3,8];
		print IOG "\t$ki\t$count{$ki}\t$ide";
		delete $hash{$ki};
	}
	close UC;
	print IOG "\n";
	close IFA;
	unlink $tout;
	unlink $teuc;

	my $topseq;
	my @title;
	my $acum=0;
	my $max;
	TTT:for my $key (sort {$count{$b} <=> $count{$a}} keys %hash){
		$acum++;
		if ($acum==1){
			$max=$count{$key};
			print LOG "\n$na";
		}
		my @title=split /\>/,$hash{$key};
		my $rate=$count{$key}/$max;
		last TTT if ($rate<$cut);
		last TTT if ($acum>=3 && $count{$key}<10);
		print OUT "$nam\t$count{$key}\t$key\n@title\n";
		print LOG "\t$count{$key}";
	}
}
close IN;
close LOG;
close IOG;
close OUT;
