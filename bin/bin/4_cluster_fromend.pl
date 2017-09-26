#! usr/bin/perl -w
use strict;

die "perl <fasta list> <outputdir> " unless (@ARGV==2);
my %cluf;
my %clur;
my $odir=$ARGV[1];

system "mkdir $odir" unless (-d $odir);
system "mkdir $odir/ends" unless (-d "$odir/ends");

open LI, "$ARGV[0]" || die $!;
open OUT, ">$odir/ends/rawends.fas" || die $!;
while(my $ff=<LI>){
	chomp($ff);
	my ($tcov,$tseq) = &cluster ($ff);
	if (-z $ff){$tcov=0;$tseq="N";}
	my $id = $ff;
	$id=~s/.*\///g;
	$id=~s/\..*$//g;
	if ($id =~ s/^For//){
		$cluf{$id}[0]=$tcov;
		$cluf{$id}[1]=$tseq;
		#	print $tcov,$tseq;
	}elsif($id =~ s/^Rev//){
		$clur{$id}[0]=$tcov;
		$tseq=reverse $tseq;
		$tseq=~tr/ATCG/TAGC/;
		$clur{$id}[1]=$tseq;
	}
#	print "$tcov,$id\n";
}

warn "not equ in For and Rev" unless (keys %cluf == keys %clur);
my %forkey= (keys %cluf >= keys %clur)?  %cluf :  %clur;
for my $key (keys %forkey){
#	print "$key\n";
	if (!exists$clur{$key}[0]){
		$clur{$key}[0]="0";
		$clur{$key}[1]="O";
	}
	if (!exists $cluf{$key}[0]){
		$cluf{$key}[0]="0";
		$cluf{$key}[1]="O";
	}
	print OUT ">$key;$cluf{$key}[0]\_$clur{$key}[0]\n$cluf{$key}[1]\n$clur{$key}[1]\n"
}
close OUT;
close LI;
#####sub####
sub cluster{
	my $file = $_[0];
	open IN, "$file" || die $!;
	my %hash;
	my $maxlen=0;
	$/="\>";
	while (my $rr=<IN>){
		chomp($rr);
		next if ($rr eq "");
		my @a=split /\n/,$rr,2;
		$a[1]=~s/\n//g;
		my @s=split /\s*/,$a[1];
		$maxlen = ($#s+1) if ($maxlen < ($#s+1));
		for my $i (0..$#s){
			if (exists $hash{$i}{$s[$i]}){
				$hash{$i}{$s[$i]}++;     #numbers of ATCG locate in $i
			}else{
				$hash{$i}{$s[$i]}=1;
			}
		}
	}
	close IN;
	$/="\n";
	my $eve;
	my $fin=0;
	my $seq;
	TTT:for my $m (0..($maxlen-1)){
		my $n=($maxlen-1-$m);
		my ($cov,$max) = &total (\%{$hash{$n}});
		if ($max >= 5){   #ends trimmed of < 5 read coverag
			$fin=$n+1;
			$maxlen=$n+1;
			last TTT;
		}
		else {
			$fin=$n;
		}
	}

	for my $n (0..($maxlen-1)){
		my ($cov,$max) = &total (\%{$hash{$n}});
		if ($max/$cov < 0.5 ){print "$n position may be heterozygous site at $file\n"}
		$eve+=$max;
		TNT:for my $key (keys %{$hash{$n}}){
			if ($hash{$n}{$key}==$max){
				$seq.="$key";
				last TNT;
			}
		}
	}
	my $ecov;
	if ($fin!=0){
		$ecov=int($eve/$fin);
	}
	else {
		$ecov=0;
		$seq="O";
	}

	return ($ecov,$seq);
}
sub total {
	my $su = $_[0];
	my %sh=%$su;
	my $tn=0;
	my $mn=0;
	for my $sk (keys %sh){
		$tn+=$sh{$sk};
		$mn = $sh{$sk} if ($mn < $sh{$sk})
	}
	return ($tn,$mn)
}
