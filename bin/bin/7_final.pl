#! usr/bin/perl -w
use strict;

die "perl <6outdir> <pre> <odir> " unless (@ARGV==3);

system "mkdir $ARGV[2] " unless(-d "$ARGV[2]");

my @files = `find $ARGV[0] -name "*F"`;
my %files;
foreach my $file (@files){
	my $num=$file;
	$num=~s/.contig.*//g;
	$num=~s/.*bar//g;
	$files{$num}=$file;
}

my %out;
my %cons;
open (CN,">","$ARGV[2]/nontarget.fa");
open (CNLOG,">","$ARGV[2]/$ARGV[1].note.txt");
print CNLOG "This is a ID list in the case that the most abundant 'end-sequences' are NOT 10 times higher to the second most ones:\n";
foreach my $key (sort {$a<=>$b} keys %files){
	open (FA,"<$files{$key}") || die $!;
	$/="\>";<FA>;$/="\n";
	while (<FA>){
		chomp (my $title=$_);
		my ($val,$kmer)=(split /\;/,$title)[1,4];
		$kmer=~s/k=//g;
		$/="\>";
		chomp(my $seq = <FA>);
		$/="\n";
		$seq=~s/\n//g;
		my $len=length $seq;
		my $id=(split /\_/,$title)[0];
		
		if (exists $out{$id}{seq}){
			if  ($out{$id}{val}<$val && $out{$id}{len}>680 && $out{$id}{len}<725){
				next if exists $cons{$title};
				if ($out{$id}{val}==1 && $val==2){
					print CN ">$title\;len=$len\n$seq\n";
					my $logname=$title;
					$logname=~s/_.*//g;
					print CNLOG "$logname\n";
					$cons{$title}=1;
				}
				next;}
		}
		
		$out{$id}{seq}=$seq;
		$out{$id}{name}=$title;
		$out{$id}{val}=$val;
		$out{$id}{kmer}=$kmer;
		$out{$id}{len}=$len;
		
	}
	close FA;
}
close CN;
close CNLOG;

open OT ,">$ARGV[2]/$ARGV[1].final.fa" || die $!;
for my $id (sort {$a<=>$b} keys %out){
	my $len=length $out{$id}{seq};
	print OT ">$out{$id}{name}\;l=$len\n$out{$id}{seq}\n";
}
close OT;
