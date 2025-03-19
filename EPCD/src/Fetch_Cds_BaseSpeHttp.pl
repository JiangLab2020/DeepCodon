###   此程序用于基于物种名下载对应的cds文件
###   物种输入格式为Fetch_SpeHttp_Based_Phylogeny.pl输出的格式, 第一列为物种名, 第二列为http下载链接


use strict;
use warnings;


my $SpeHttp = $ARGV[0];

my $Path = `pwd`;  $Path=~s/\n$//;

my @b=split(/[\/\.]/,$SpeHttp);  my $Target = $b[$#b-1];

my $TargetPath=$Path;
if ( $SpeHttp =~ /\// ) {
	$TargetPath=join("/",@b[0..$#b-2]);
}

my $Type    = 3;   ###  1表示下载gff文件, 2表示下载genome文件, 3表示下载cds文件, 4表示下载gbff文件, 5表示下载pep文件, 6表示下载所有蛋白功能注释文件
if ( $ARGV[1] ) {
	$Type   = $ARGV[1];
}

my %TList = ();

$TList{1}="gff";  $TList{2}="genome";  $TList{3}="cds";  $TList{4}="gbff";  $TList{5}="pep";  $TList{6}="pep";

my $work_path = $Target."_".$TList{$Type};
if ( !(-e $work_path) ) {
	`mkdir $work_path`;
}

my $k = 0;

print "Download $TList{$Type} files #############################################################################","\n\n";

my %AList = ();

if ( -e "$TargetPath/$Target\.cds" ) {
	open IN,"$TargetPath/$Target\.cds";
	while (<IN>) {
		chomp($_);
		if ( $_ =~ /http\=\"(.*)\"/ ) {
			$AList{$1}+=1;
		}
	}
	close(IN);
}

open OUTPUT,">>$TargetPath/$Target\.cds";

open INPUT,$SpeHttp;
while (<INPUT>) {
	chomp($_);
	my @a = split(/\t/,$_);
	if ( !$AList{$a[1]} ) {

		$AList{$a[1]} = 2;

		my $species = $a[0];
		my $ftp = $a[1];
		my @f = split(/\//,$ftp);

		my $Http = "";  my $HttpFile = "";

		if ( $Type eq 1 ) {
			$Http = "$ftp/$f[$#f]\_genomic.gff.gz";
			$HttpFile  = "$f[$#f]\_genomic.gff.gz";
		}
		if ( $Type eq 2 ) {
			$Http = "$ftp/$f[$#f]\_genomic.fna.gz";
			$HttpFile  = "$f[$#f]\_genomic.fna.gz";
		}
		if ( $Type eq 3 ) {
			$Http = "$ftp/$f[$#f]\_cds_from_genomic.fna.gz";
			$HttpFile  = "$f[$#f]\_cds_from_genomic.fna.gz";
		}
		if ( $Type eq 4 ) {
			$Http = "$ftp/$f[$#f]\_genomic.gbff.gz";
			$HttpFile  = "$f[$#f]\_genomic.gbff.gz";
		}
		if ( $Type eq 5 or $Type eq 6 ) {
			$Http = "$ftp/$f[$#f]\_protein.faa.gz";
			$HttpFile  = "$f[$#f]\_protein.faa.gz";
		}

		print $Http,"\n";

		`wget -q $Http`;
		
		if ( -e "$HttpFile" ) {
			`gunzip $HttpFile`;  $HttpFile =~ s/\.gz$//;
			`mv $HttpFile $work_path/$Target\_$k\.$TList{$Type}`;

			chdir($work_path);

			if ( $Type eq 1 ) {
				`NCBI_GffMaxGene.pl $Target\_$k\.gff`;
				`Gff2Loc.pl $Target\_$k\_Max.gff $Target\_$k\.txt`;

				my %TList = ();
				open IN,"$Target\_$k\_Max.txt";
				while (<IN>) {
					chomp($_);
					my @a = split(/\t/,$_);
					$TList{$a[0]}=$a[3];
				}
				close(IN);
				
				open IN,"$Target\_$k\.txt";
				while (<IN>) {
					chomp($_);
					my @a = split(/\t/,$_);
					if ( $TList{$a[0]} ) {
						$a[0]=$TList{$a[0]};
						print OUTPUT join("\t",@a),"\n";
					}
				}
				close(IN);

				unlink("$Target\_$k\.gff");
				unlink("$Target\_$k\_Max.gff");
				unlink("$Target\_$k\.txt");
				unlink("$Target\_$k\_Max.txt");
			}

			if ( $Type eq 3 ) {
				open IN,"$Target\_$k\.cds";
				while (<IN>) {
					chomp($_);
					if ( $_ =~ /^\>(.*?) / ) {
						print OUTPUT ">",$1," ","[$species]"," ","http=\"$ftp\"","\n";
					}
					else {
						print OUTPUT $_,"\n";
					}
				}
				close(IN);
			}

#			if ( $Type eq 5 ) {
#				my $seqio = new Bio::SeqIO(-file => "$Target\_$k\.pep", -format => 'fasta');
#				while ( my $seq = $seqio->next_seq ) {
#					print OUTPUT ">",$seq->display_id,"\t",$seq->description,"\n",$seq->seq(),"\n";
#				}
#			}
#			if ( $Type eq 6 ) {
#				my $seqio = new Bio::SeqIO(-file => "$Target\_$k\.pep", -format => 'fasta');
#				while ( my $seq = $seqio->next_seq ) {
#					print OUTPUT $seq->display_id,"\t",$seq->description,"\n";
#				}
#			}
			
			if ( $k > 100 ) {
				`rm -f $Target\_$k\.$TList{$Type}`;
			}
		}
		else {
			print $a[1],"\n";
		}
		chdir($Path);
		$k++;
	}
}
close(INPUT);


print $k,"\n";

close(OUTPUT);

