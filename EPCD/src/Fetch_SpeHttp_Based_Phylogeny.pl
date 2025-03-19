###   此程序输出某一个分类下所有已经测序的物种

use strict;
use warnings;

my $phylogeny = $ARGV[0];

my $Http_File = $ARGV[1];

my $Taxion_File = "Database/Taxion/TaxPhylogeny_Name.txt";
if ( $ARGV[2] ) {
	$Taxion_File = $ARGV[2];
}

my %PList = (); 
open INPUT,$Taxion_File;
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_);
	if ( $_ =~ /\t$phylogeny\t/ or $_ =~ /\t$phylogeny$/ ) {
		$PList{$a[$#a]} = $_;
	}
}
close(INPUT);

my %HList = ();

open OUTPUT,">$phylogeny\_SeqSpe.txt";

open INPUT,$Http_File;
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_);
	my $species = $a[7];
	if ( $species and $PList{$species} and !$HList{$a[19]} ) {
		$HList{$a[19]} = 1;
		my $GCF_Http = $a[19];
		if ($a[19]=~ /(GCA\_.*?)\_/) {
			my $GCA=$1;
			if ($a[17]=~ /GCF/) {
				$GCF_Http=~ s/\/GCA\//\/GCF\//;
				$GCF_Http=~ s/$GCA/$a[17]/;
			}
		}
		print OUTPUT $species,"\t",$GCF_Http,"\t",$PList{$species},"\n";
	}
}
close(INPUT);

close(OUTPUT);
