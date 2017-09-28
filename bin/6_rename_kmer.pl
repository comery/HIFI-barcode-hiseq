#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);

die "perl <assemble.contig> <end.fasta> <outdir>" unless (@ARGV==3);
open ASS, "$ARGV[0]" || die $!;
open OLD, "$ARGV[1]" || die $!;

system "mkdir $ARGV[2]" unless -e $ARGV[2];
my $base=basename $ARGV[0];
open ASF, ">$ARGV[2]\/$base.F" || die $!;
open POI, ">$ARGV[2]\/$base.F.add" || die $!;

my %name;
my $count=0;
while(<OLD>){
	chomp;
	next unless (s/^>//);
	$count++;
	$name{$count}=$_;	
}
my %aha;
$/="\>"; <ASS>; $/="\n";
while(<ASS>){
        chomp;
		my $lk=$_;
		$lk=~s/>//;
		$lk=~s/_.*//g;
        $/="\>";
        chomp(my $seq = <ASS>);
        $/="\n";
        $seq=~s/\n//g;
        my $len=length $seq;
        my @a=split /\t/;
        my @b=split /\_/, $a[0];
        $aha{$b[1]}=0 unless(exists $aha{$b[1]});
        if ($b[2]==1){
                $aha{$b[1]}++;
                print ASF ">$name{$b[1]};k=$lk\n$seq\n";
        }elsif($aha{$b[1]}==0){
                $aha{$b[1]}++;
                print ASF ">$name{$b[1]};k=$lk\n$seq\n";
        }else{	
				my $length=length $seq;
				next unless ($length>680 && $length<725);
				my $outnum=$b[2]-1;
                print POI ">$name{$b[1]};k=$lk;AddSeq$outnum\n$seq\n";
        }
}
close ASS;
close ASF;
close POI;
