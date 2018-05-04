#! usr/bin/perl -w
use strict;
use File::Basename qw(dirname basename);
use FindBin qw($Bin);

die "perl <4outrawends> <pre> <odir> <cpu>  " unless (@ARGV==4);
my %cluf;
my %clur;
my $dir=$ARGV[2];
my $idir=dirname $ARGV[0];
my $cpu=$ARGV[3];
system "mkdir $dir" unless (-d $dir);
system "mkdir $dir/ends" unless (-d "$dir/ends");
system "mkdir $dir/asm" unless (-d "$dir/asm");
system "mkdir $dir/mid" unless (-d "$dir/mid");

open FA, "$ARGV[0]" || die $!;
my (%depth,%for,%rev,%sum,%val);

while (<FA>){
	chomp;
	my $name=$_;
	$name=~s/\>//g;
	my ($l1,$l2)=(split /;/,$name)[0,1];
	my ($id,$sort)=(split /_/,$l1)[0,1];
		
	chomp ($for{$id}{$sort}=<FA>);
	chomp ($rev{$id}{$sort}=<FA>);

	my $lfor=length $for{$id}{$sort};
	my $lrev=length $rev{$id}{$sort};
	if (($for{$id}{$sort}=~/O|N/) || ($lfor <130)){delete $for{$id}{$sort};}
	if (($rev{$id}{$sort}=~/O|N/) || ($lrev <130)){delete $rev{$id}{$sort};}  #140 to 130

	my ($c1,$c2)=(split /_/,$l2)[0,1];
	$depth{$id}{$sort}{1}=$c1;
	$depth{$id}{$sort}{2}=$c2;
	
}
close FA;

#get value for each pair
open (TA,">","$dir/ends/$ARGV[1].depth.xls");
my $idn=0;
TTT:for my $id (sort {$a<=>$b} keys %for){
		for (my $sor1=1;$sor1<=2;$sor1++){
			 for (my $sor2=1;$sor2<=2;$sor2++){
				my $sor=join ("_",$sor1,$sor2);
				$val{$id}{$sor}=&val($sor1,$sor2);
				if ((! exists $for{$id}{$sor1}) || (! exists $rev{$id}{$sor2})) {  
					delete $val{$id}{$sor} ;
					next;
				}
				$idn++;
				my $dep=join ("_",$depth{$id}{$sor1}{1},$depth{$id}{$sor2}{2});
				print TA "$idn\t$id\t$sor\t$dep\t$val{$id}{$sor}\n";
		}
	}
}
close TA;

#get topn pairs
for my $id (sort {$a<=>$b} keys %val){
	for my $sor (sort {$val{$id}{$b}<=>$val{$id}{$a}} keys %{$val{$id}}){
		my ($sor1,$sor2)=(split /\_/,$sor)[0,1];
		my $len1=length $for{$id}{$sor1};
		my $len2=length $rev{$id}{$sor2};
		my $cov=int(($depth{$id}{$sor1}{1}*$len1+$depth{$id}{$sor2}{2}*$len2)/($len1+$len2));
		my $val=$val{$id}{$sor};
		open OUT, ">>$dir/ends/$ARGV[1].$val.fa" || die $!;
		print OUT ">$id\_$sor1\-$sor2\;$val;$depth{$id}{$sor1}{1}\_$depth{$id}{$sor2}{2};size=$cov\n$for{$id}{$sor1}\n$rev{$id}{$sor2}\n";
		close OUT;
	}
}

open (MIDLIS,">","$dir/mid/mid.lis") ||die $!;
print MIDLIS ">\nf=$dir/mid/mid.fa\n";
close MIDLIS;

open (BA,">>","$dir/../final.sh") || die $!;
#obtain mid.fa
print BA "$Bin/cmr -a $dir/../step1/$ARGV[1]\_mid_1.fasta  -b $dir/../step1/$ARGV[1]\_mid_2.fasta  -o $dir/mid/mid.fa -2 $dir/mid/$ARGV[1].midfail1 -3 $dir/mid/$ARGV[1].midfail2 -l 30 -u 120 -c 0.95 -m 0\n";
#gap close
for (my $i=1;$i<=4;$i++){
	print BA "$Bin/barcode  -e $dir/ends/$ARGV[1].$i.fa -r $dir/mid/mid.lis  -o $dir/asm/bar$i  -x 470 -n 220 -l 100 -v 8 -k 127 -t $cpu \n" ;
	print BA "perl $Bin/6_rename_kmer.pl $dir/asm/bar$i.contig $dir/ends/$ARGV[1].$i.fa $dir/asm\n";
}

print BA "perl $Bin/7_final.pl $dir/asm  bar $dir/../Result\n" ;
close BA;

sub val{
	my ($sor1,$sor2)=@_;
	my $val;
	if ($sor1==1 && $sor2==2){$val=3;}
	if ($sor1==2 && $sor2==1){$val=4;}
	if ($sor1==2 && $sor2==2){$val=2;}
	if ($sor1==1 && $sor2==1){$val=1;}
	return $val;
}

