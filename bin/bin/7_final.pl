#! usr/bin/perl -w
use strict;

die "perl <6outdir> <pre> <odir> " unless (@ARGV==3);

system "mkdir $ARGV[2] " unless(-d "$ARGV[2]");

my @files = `find $ARGV[0] -name "*F"`;
my %files;
my %adds;
foreach my $file (@files){
	my $num=$file;
	$num=~s/.contig.*//g;
	$num=~s/.*bar//g;
	$files{$num}=$file;
	$adds{$num}="$file.add";
	$adds{$num}=~s/\n//;
}

my %out;
my %cons;
open (CN,">","$ARGV[2]/nontarget.fa");
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

my %addfa;
foreach my $key (sort {$a<=>$b} keys %adds){
	open (AD,"<$adds{$key}") || die $!;
    $/="\>";<AD>;$/="\n";
    while (<AD>){
        chomp (my $title=$_);
        my $addnum=$title;
        $addnum=~s/.*AddSeq//g;
        $title=~s/;AddSeq.*//g;
		$/="\>";
		chomp (my $seq=<AD>);
		$/="\n";
        $addfa{$title}{$addnum}=$seq;
    }
    close AD;
}

open OT ,">$ARGV[2]/$ARGV[1].final.fa" || die $!;
open ADD,">$ARGV[2]/$ARGV[1].additional.fa" || die $!;

for my $id (sort {$a<=>$b} keys %out){
	my $len=length $out{$id}{seq};
	my $title=$out{$id}{name};
	print OT ">$title\;l=$len\n$out{$id}{seq}\n";
	next unless exists $addfa{$title};
	foreach my $addnum (sort {$a <=> $b} keys %{$addfa{$title}}){
		print ADD ">$title;AddSeq$addnum\n$addfa{$title}{$addnum}\n";
	}			

}
close OT;
close ADD;

open LOG,">$ARGV[2]/$ARGV[1].note.txt" || die $!;
print LOG "##### NOTE #####
The proposed method assemble short-read Illumina sequences 
based on k-mer sequence matches, and such misassembly was 
not observed in our published real data, but it's still 
possible theoretically. 
Thus, we produce additional sequences (.additional.fa)
which have similar or same scores comparing to their 
best alternative.\n"; 
close LOG;
