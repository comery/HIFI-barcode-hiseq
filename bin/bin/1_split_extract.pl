#! usr/bin/perl -w
use strict;
=head1 Description

        Usage: perl extract.pl <parameter>

        --fq1           "xx.1.fq"
        --fq2           "xx.2.fq"
        --pri           "primer set with index ahead, format in "ID \t seq";
        --out           "prefix of the output file"
        --len           "the length of index | by default =0"
        --help          "print out this information"

=cut

use Getopt::Long;
use FindBin qw($Bin $Script);
use strict;


my ($Fq1,$Fq2,$Pri,$Out,$Len,$Help);

my%primer=(
        "V" => "(A|C|G)",
        "D" => "(A|T|G)",
        "B" => "(T|G|C)",
        "H" => "(A|T|C)",
        "W" => "(A|T)",
        "S" => "(C|G)",
        "K" => "(T|G)",
        "M" => "(A|C)",
        "Y" => "(C|T)",
        "R" => "(A|G)",
        "N" => "(A|T|C|G)",
);
GetOptions(
        "fq1:s"=>\$Fq1,
        "fq2:s"=>\$Fq2,
        "out:s"=>\$Out,
        "pri:s"=>\$Pri,
        "len:i"=>\$Len,
         "help"=>\$Help
);
die `pod2text $0` if ($Help || !defined ($Fq1) || !defined ($Fq2)||!defined ($Pri));
$Len ||= 0;
if ($Fq1=~/gz$/){
                open (FFQ, "<:gzip",$Fq1) || die $!;
        }else{
                open (FFQ, $Fq1) || die $!;
        }
if ($Fq2=~/gz$/){
                open (RFQ,"<:gzip",$Fq2) || die $!;
        }else{
                open (RFQ, $Fq2) || die $!;
        }

my $priwc;
$priwc = `less -S $Pri|wc -l`;
die "the primer file need to have each forward and reverse primer"
unless ($priwc%4 == 0);

die "the primer file contain more than 2 primer sets, however, without index length hava been seted" if ($priwc > 4 && $Len==0);
print "the len is $Len\n";

my %indp;
open PRI, "$Pri" || die $!;
my ($plenf, $plenr,%prif,%prir);
while (<PRI>){
	chomp;
	my @a=split /\t/;
	die "primer set is not in correct format\n" unless (@a==2);
	my $ipr = convert ($a[-1]);
	if (exists $indp{$a[0]}){
		warn "double exits $a[0] in primer set\n";
	}else{
		$indp{$a[0]}=$ipr;
	}
}
close PRI;
my %FH;

open (LS,">","$Out.splitlist");
for my $key (keys %indp){
	print LS "$Out\_$key\.fasta\n";
}
close LS;

for my $key (keys %indp){
	open ($FH{"$key"}, ">", "$Out\_$key\.fasta") || die $!;
}
sub convert {
        my $old = $_[0];
        my @a = split /\s*/,$old;
        my $primer;
        for (my $i=0; $i<=$#a; $i++){
                if (exists $primer{$a[$i]}){
                        $primer.=$primer{$a[$i]};
                }else{
                        $primer.=$a[$i];
                }
        }
        return $primer;
}
open ERR, ">$Out\_err.fasta" || die $!;
open MIF, ">$Out\_mid_1.fasta" || die $!;
open MIR, ">$Out\_mid_2.fasta" || die $!;
my $seqnum=0;
TNT:while(my $t1=<FFQ>){
        $seqnum++;
        chomp($t1);
        chomp(my $seqf=<FFQ>);
        <FFQ>;
        <FFQ>;
        chomp (my $t2=<RFQ>);
        chomp (my $seqr=<RFQ>);
        <RFQ>;
        <RFQ>;
	my $fmatch=0;
	my $rmatch=0;
	my $fpri;
	my $rpri;
	TTT:for my $key (keys %indp){
		if ($seqf=~/^$indp{$key}/){
			$fmatch++;
			$fpri=$key;
			my $ff=$seqf;
			$ff=~s/^$indp{$fpri}//;
			if ($ff=~/$indp{$key}/){
				$fmatch=2;
			}
		}
		elsif ($seqf =~/$indp{$key}/){
			$fmatch=2;
		}
		if ($seqr=~/^$indp{$key}/){
			$rmatch++;
			 $rpri=$key;
			 my $rr=$seqr;
			 $rr=~s/^$indp{$key}//;
			 if ($rr=~/$indp{$key}/){
				 $rmatch=2;
			 }
		}
		elsif ($seqr =~/$indp{$key}/){
			$rmatch=2;
		}
	}
	if ($fmatch ==1 and $rmatch ==0){
		$seqf=~s/^$indp{$fpri}//;
		print {$FH{$fpri}} ">$t1\n$seqf\n";
		next TNT;
	}elsif($fmatch ==0 and $rmatch ==1){
		$seqr=~s/^$indp{$rpri}//;
		print {$FH{$rpri}} ">$t2\n$seqr\n";
		next TNT;
	}elsif ($fmatch >=1 or $rmatch >=1){
		print ERR ">$seqnum\_f\n$seqf\n>$seqnum\_r\n$seqr\n";
		next TNT;
	}else{
		print MIF ">$seqnum\_f\n$seqf\n";
		print MIR ">$seqnum\_r\n$seqr\n";
	}
}	
close ERR;
close MIR;
close MIF;
#for my $key (keys %FH){
#	close {$FH{$key}}
#}
print "done\n";
