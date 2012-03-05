#!/usr/bin/perl -w
 
use strict;
use List::MoreUtils qw/all any each_array/;
use Scalar::Util qw/looks_like_number/;
use Statistics::Basic qw/mean stddev/;
 
######################
#Microarray Filter and Fold Change Finder
######################
 
#Open data file and read into array:
 
print "\nMicroarray Filter and Analysis Tool:\n";
 
 
if (@ARGV != 1) {
        die ("\nUse: perl microarray.pl <Input datafile.txt>");
}
 
my $file = $ARGV[0];
 
open(my $in, "<", $file) or die "\nCouldn't open file $file: $!";
 
my $topLine;
my @nameArray;
my @dataArray =();
my $indexLine = 0;
my $sampleNum = 0;
 
 
while (my $line = <$in>) {
        if ($line =~ /^probes/) {
                chomp $line;
                $topLine = $line;
                $indexLine = -1;
        }
        elsif ($line =~ /probes/) {
            die "File '$file' is malformed at line $.";
        }
        else {
                chomp $line;
                my @tempArray = split(" ",$line);
                $nameArray[$indexLine] = shift @tempArray;
                $sampleNum = @tempArray;
                die "File '$file' contains non-numeric data at line $."
                        unless all { looks_like_number($_) } @tempArray;
                @{$dataArray[$indexLine]} = @tempArray;
        }
        $indexLine++;
}
 
close $in;
 
my $geneNumber = $indexLine;
 
my @filterNames;
my @filterData;
 
my $ea = each_array(@nameArray, @dataArray);
while (my ($name, $gene) = $ea->())  {
        if (any { $_ > 300 } @$gene) {
                push (@filterNames,$name);
                push (@filterData,$gene);
        }
}
 
my $filterNumber = scalar @filterNames;
 
print "\nThere are $filterNumber genes that meet filter criteria.\n";
 
my $fldScore = 0;
my %scoreHash = ();
my $reporter = 0;
my $incrementor = 0;
 
 
for (my $i = 0; $i < $filterNumber; $i++) {
        my @controlArray;
        my @sampleArray;
        for (my $j = 0; $j < 20; $j++) {
                push (@controlArray,$filterData[$i][$j]);
        }
        for (my $k = 20; $k < 41; $k++) {
                push(@sampleArray, $filterData[$i][$k]);
        }
        my $controlMean = mean(\@controlArray);
        my $controlSD = stddev(\@controlArray);
        my $sampleMean = mean(\@sampleArray);
        my $sampleSD = stddev(\@sampleArray);
        my $fldNum = $controlMean - $sampleMean;
        my $fldDenom = $controlSD + $sampleSD;
        $fldScore = $fldNum / $fldDenom;
        $scoreHash{$fldScore} = $i;
        $reporter++;
        $incrementor++;
        print "\nFLD score: $fldScore";
        if ($incrementor == 100) {
                print "\nCurrent cycle: $reporter";
                $incrementor = 0;
        }
}
 
my $scoreCounter = 1;
 
 
foreach my $key (sort keys %scoreHash) {
        print "\nTop Ranking Differentially Expressed Genes:";
        print "\n$scoreCounter. $filterNames[$scoreHash{$key}]";
        $scoreCounter++;
}
 
