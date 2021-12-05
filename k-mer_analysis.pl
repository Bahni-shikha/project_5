#!perl/bin/usr
#Bahnishikha
#24th Nov 2021
#this this the analysis for tetramer
#one can analyse for other k-mers also just by changing the number in place of '4' to the number of their interest in line no. 27, 29, 41 and 43. 
open FILE, "intron_IME_data.fasta" or die "error $!";
open OUT, "> pb_p5_out.txt" or die "error $!";
chomp(my @intron=<FILE>);
my $seq=join("", @intron);
my @seq = split(">", $seq);
shift(@seq);

for($i=0;$i<=scalar(@seq);$i++)
{
if (@seq[$i]=~m/i1_/){push(@early, @seq[$i]);}  #push function will store each of i1 seq in @ealy array even @ealy is printed outside the loop
elsif (@seq[$i]=~m/i2_/){push(@late, @seq[$i]);}
}

#for early intron
for($i=0;$i<=scalar(@early);$i++)
{
  @sequence = split('CDS|5UTR|3UTR',$early[$i]);  #it will split the sequence information and the sequence of early intron
 push(@tog,$sequence[1]);
}

$nuc=join("",@tog);
for ($n=0;$n<=length($nuc)-4;$n++)
{
$str=substr($nuc,$n,4);
push (@str,$str);
}

#for late intron
for($i=0;$i<=scalar(@late);$i++)
{
  @Lt = split('CDS|5UTR|3UTR',$late[$i]);  #it will split the sequence information and the sequence of early intron
 push(@toglt,$Lt[1]);
}

$nuclt=join("",@toglt);
for ($n=0;$n<=length($nuclt)-4;$n++)
{
$string=substr($nuclt,$n,4);
push (@strings,$string);
}


foreach $str(@str){$count1{$str}++;}
foreach $str(@strings){$count2{$str}++;}
print OUT "k-mer"."  \t"."i1"." \t "."Not_i1"." \t "."Odds"."  \t"."Log odds"."\n";
foreach $str(sort keys %count1,sort keys %count2){print OUT $str ." \t ". $count1{$str}."  \t". $count2{$str}."  \t",sprintf ("%.3f",($count1{$str}/$count2{$str})). "  \t",
 sprintf ("%.3f", log($count1{$str}/$count2{$str})/log(2))."\n";}


close (FILE);
close (OUT);
