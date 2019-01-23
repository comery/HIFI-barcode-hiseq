#! usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename qw(dirname basename);
die "perl <unique_list generated in step2> <fq1> <fq2> <outdir>" unless (@ARGV==4);

open IN, "$ARGV[0]" || die $!;

if ($ARGV[1]=~/gz$/){
                open (FQ1, "<:gzip",$ARGV[1]) || die $!;
        }else{
                open (FQ1, $ARGV[1]) || die $!;
        }
if ($ARGV[2]=~/gz$/){
                open (FQ2,"<:gzip",$ARGV[2]) || die $!;
        }else{
                open (FQ2, $ARGV[2]) || die $!;
        }

my $outdir=$ARGV[3];
system "mkdir $outdir " unless(-d "$outdir");
system "mkdir $outdir/03.overlap" unless (-d "$outdir/03.overlap");
system "mkdir $outdir/03.overlap/reads" unless (-d "$outdir/03.overlap/reads");
system "mkdir $outdir/lis " unless(-d "$outdir/lis");


my %hash;
my %mid;
my %cmr;
my %ha;
my $n=1;

while (<IN>){
	chomp;
	my $name=(split /\t/)[0];
	$name=~s/.*\///g;
	my $midt = (split /\.|\_/,$name)[1];
	if (exists $ha{$midt}){
		$n++;}
	else {$n=1;$ha{$midt}=1;}
	$midt="$midt\_$n";
	my $mid1 = "$midt\_1";
	my $mid2 = "$midt\_2";
	$cmr{$midt}[0]=$mid1;
	$cmr{$midt}[1]=$mid2;
	$mid{$mid1}=1;
	$mid{$mid2}=2;
	chomp (my $read=<IN>);
	my @a=split /\@/,$read;
	TN:foreach my $ele (@a){
		$ele=~s/\t$//;
		$ele=~s/\s+$//;
		$ele=~s/^\s+//;
		if ($ele eq "" or $ele eq " "){next TN};
		if (exists $hash{$ele}){
			warn "$ele exists more than once\n"
		}else{
			$hash{$ele}="$midt";
		}
	}
	
}
close IN;
my $num=keys %hash;
print "there are total $num reads\n";

my $seqnum=0;
my $exinum=0;
TNT:while (my $t1=<FQ1>){
	$seqnum++;
	chomp($t1);
	chomp (my $f=<FQ1>);
	<FQ1>;
	<FQ1>;
	chomp (my $t2=<FQ2>);
	chomp (my $r=<FQ2>);
	<FQ2>;
	<FQ2>;
	$t1 =~s/^@//;
	$t2 =~s/^@//;
	if (exists $hash{$t1}){
		$exinum++;
		open (O1,">>","$outdir/03.overlap/reads/$hash{$t1}\_1.fasta");
		open (O2,">>","$outdir/03.overlap/reads/$hash{$t1}\_2.fasta");
		
		print O1 ">$t1\n$f\n";
		print O2 ">$t2\n$r\n";
		close O1;
		close O2;
		next TNT;
	}
	if (exists $hash{$t2}){
		$exinum++;
		open (O1,">>","$outdir/03.overlap/reads/$hash{$t2}\_1.fasta");
		open (O2,">>","$outdir/03.overlap/reads/$hash{$t2}\_2.fasta");
		print O1 ">$t2\n$r\n";
		print O2 ">$t1\n$f\n";
		close O1;
		close O2;
	}
}
print "there find $exinum reads\n";
close FQ1;
close FQ2;
open (OT,">","$outdir/lis/03overfas.lis");
for my $key (keys %cmr){
	my $file1="$outdir/03.overlap/reads/$cmr{$key}[0]\.fasta";
	my $file2="$outdir/03.overlap/reads/$cmr{$key}[1]\.fasta";
	my $out="$outdir/03.overlap/$key\.fasta";
	`$Bin/cmr -a $file1 -b $file2 -o $out -2 $outdir/03.overlap/fail1 -3 $outdir/03.overlap/fail2 -l 30 -u 120 -c 0.95 -m 0`;#add by yangchentao
	print "$key has been overlapted\n";
	print OT "$out\n";
}
close OT;
print "done\n"
