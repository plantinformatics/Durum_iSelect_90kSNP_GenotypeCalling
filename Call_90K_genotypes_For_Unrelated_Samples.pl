#!/usr/bin/perl -w

use IO::File;
use File::Spec::Functions;
use Math::Trig;
use Getopt::Std;

#######################################################################################
# Get input filenames from the DOS command line.
#######################################################################################

## Get options from ARGV
my %Options = ();
getopts("c:f:d:s:r:o:h", \%Options);

my $CTMLfilename = $Options{c};
my $FDTfilename = $Options{f};
my $DistanceLimit = $Options{d};
my $ConfidenceCutoff = $Options{s};
my $MAF = $Options{r};
my $OutputFilename = $Options{o};

if ($Options{h} or !defined ($Options{c} and $Options{f} and $Options{r} and $Options{o})) {
	print <<EOF;

##################################################################################################################################

	Call_90K_genotypes_For_Unrelated_Samples.pl

	This scripts generates genotype calls for unrelated samples assayed using the iSelect 90K wheat SNP bead chip.
	
	Genotype calls (A and B, where A denotes adenosine or thymine; and B denotes cytosine and guanine) are generated for 
	clusters that unambigiously associate with genetically map loci. 

	Options for running this perl script
	
	  -c   Folder path and name of file containing cluster information; *_ClusterToMappedLoci.txt
	  
	  -f   Folder path and name of FullDataTable exported from GenomeStudio
	       Format: Text-delimited <SNPindex><SNPname><Sample1.theta><Sample1.normR><Sample2.theta><Sample2.normR> etc
		   
	  -d   Distance limit for cluster assignment (in standard deviations). Default 2
	  
	  -s   Confidence score for genotype call(value between 0 and 1, where 1 is highest) Default 0.7
	  
	  -r   Only report polymorphic SNPs in *_GENO_InformativeSNPtable and *_GENO_InformativeSNPreport
	       when they have MAF > user-specified value. Use 0 to report all SNPs
		   
	  -o   Folder path and name of output file
	  
##################################################################################################################################

EOF
	exit(0);
}	  

	  
#######################################################################################
# Fixed parameters
#######################################################################################

if (!defined $DistanceLimit) { $DistanceLimit = 2; } 
if (!defined $ConfidenceCutoff) { $ConfidenceCutoff = 0.7; } 


#######################################################################################
# Main Program
#######################################################################################

# Create output files
open (GENOreport, '>', $OutputFilename."_GENO_CompleteReport.txt") or die;
print GENOreport $OutputFilename."_GENO_CompleteReport.txt\n"."Distance limit for cluster assignment = ".$DistanceLimit.";   Confidence Score for genotype call = ".$ConfidenceCutoff."\n";
print GENOreport "Index\tName\tLocus\tComment\t%Samples assigned to cluster\t%Samples with genotypes called\t#AA\t#BB\t#NC\n";
open (GENOtable, '>', $OutputFilename."_GENO_CompleteTable.txt") or die;
open (GENOinfoReport, '>', $OutputFilename."_GENO_InformativeSNPreport.txt") or die;
print GENOinfoReport $OutputFilename."GENO_InformativeSNPreport.txt\n"."Distance limit for cluster assignment = ".$DistanceLimit.";   Confidence Score for genotype call = ".$ConfidenceCutoff.";   Minor Allele Frequency cut-off = ".$MAF."\n";
print GENOinfoReport "Index\tName\tLocus\tComment\t%Samples assigned to cluster\t%Samples with genotypes called\t#AA\t#BB\t#NC\n";
open (GENOinfoTable, '>', $OutputFilename."_GENO_InformativeSNPtable.txt") or die;
open (ClusterPostions, '>', $OutputFilename."_ClusterPositions.txt") or die;
print ClusterPostions "SNPindex\tSNPname\tCluster#\tChromosome\tC T Mean\t C T Dev\t C R Mean\tC R Dev\tAllelicState\n";

# Process ClusterToMappedLoci file
if ($CTMLfilename =~ /\.gz$/) { open (CTML, "<:gzip", $CTMLfilename) or die; } else { open (CTML, '<', $CTMLfilename) or die; }
my $CTMLheader = <CTML>;  
while (<CTML>) {
	chomp $_;
	@_ = split ('\t', $_);
	push @ClusterData, [@_]; 
}

my $ClusterTrackingNumber = 1;
@ClusterData = sort { $a->[0] <=> $b->[0] } @ClusterData;

until (@ClusterData == 0) {
	push @Temp, shift @ClusterData;
	while (@ClusterData >0 and $Temp[0]->[0] == $ClusterData[0]->[0]) { push @Temp, shift @ClusterData; }
	@Temp = sort { $a->[3] <=> $b->[3] || $a->[5] <=> $b->[5]} @Temp;
	my $CurrentSNP = $Temp[0]->[0]; 
	push @{ $MappedClusterDetails{$ClusterTrackingNumber} }, [$Temp[0]->[2], $Temp[0]->[8] ]; 
	unshift @{ $CurrentSNP }, shift @Temp;
	${ $CurrentSNP }[0]->[9] = $ClusterTrackingNumber; 
	until (@Temp == 0) {
		if (${ $CurrentSNP }[0]->[3] != $Temp[0]->[3] and ${ $CurrentSNP }[0]->[5] != $Temp[0]->[5]) {
			$ClusterTrackingNumber++;
			push @{ $MappedClusterDetails{$ClusterTrackingNumber} }, [$Temp[0]->[2], $Temp[0]->[8] ]; 
			unshift @{ $CurrentSNP }, shift @Temp;
			${ $CurrentSNP }[0]->[9] = $ClusterTrackingNumber; 
		} else {
			push @{ $MappedClusterDetails{$ClusterTrackingNumber} }, [ $Temp[0]->[2], $Temp[0]->[8] ]; 
		}
	}
	$ClusterTrackingNumber++;
}	

# Process FullDataTable
if ($FDTfilename =~ /\.gz$/) { open (FDT, "<:gzip", $FDTfilename) or die; } else { open (FDT, '<', $FDTfilename) or die; }
my $FTDtableHeader = <FDT>; 
@FDTheader = split ('\t', $FTDtableHeader);  
print GENOtable $FDTheader[0]."\t".$FDTheader[1]."\tLocus"; 
print GENOinfoTable $FDTheader[0]."\t".$FDTheader[1]."\tLocus"; 
for (my $i=2; $i < @FDTheader; $i+=2) {   
	$SampleName = $FDTheader[$i];
	$SampleName =~ s/.Theta//;
	print GENOtable "\t".$SampleName.".Cluster\t".$SampleName.".Score\t".$SampleName.".Geno"; 
	print GENOinfoTable "\t".$SampleName.".Geno"; 
}
print GENOtable "\n";
print GENOinfoTable "\n";

my $SampleSize = (@FDTheader)/2 -1;
while (<FDT>) {
	chomp $_;
	@_ = split ('\t', $_);
	&AssignSamplesToCluster; 
}


#######################################################################################
# AssignSamplesToCluster subroutine
#######################################################################################

sub AssignSamplesToCluster {
	my %ClusterTypes = ();
	my %MappedChr = ();
	my %ClustersUsed = ();
	@ClusterDetails = sort { $a->[3] <=> $b->[3] || $a->[5] <=> $b->[5] } @{$_[0]};	
	
	foreach $c (@ClusterDetails) {
		for (my $i=0; $i < @{ $MappedClusterDetails{$c->[9]} }; $i++) {
			$ClusterTypes{uc($MappedClusterDetails{$c->[9]}[$i]->[1])}++; 
			if (exists $ClusterTypes{"A"} or exists $ClusterTypes{"B"}) { $Comment = "Mapped"; } else { $Comment = "Uninformative"; }		
			if ($MappedClusterDetails{$c->[9]}[$i]->[1] ne "C" and $MappedClusterDetails{$c->[9]}[$i]->[1] ne "U") { $MappedChr{substr($MappedClusterDetails{$c->[9]}[$i]->[0],0,2)}++; } else { $MappedChr{"n/a"}++; }
		}
	}

	foreach $Locus (sort keys %MappedChr) {
		next if ($Comment eq "Mapped" and $Locus eq "n/a");
		$MinorAlleleFreq = 0;
		$GenoCall{"A"} = 0;
		$GenoCall{"B"} = 0;
		$GenoCall{"NC"} = 0;
		$NumClustered = 0;
		@Genotypes = ();
		print GENOtable $_[0]."\t".$_[1]."\t".$Locus;  
		print GENOreport $_[0]."\t".$_[1]."\t".$Locus."\t".$Comment;  
		for ($Sample=2; $Sample	< @_; $Sample += 2) {
			my $ClosestCluster = 0;
			my $ClosestCluster_nED = 1000;
			my $NextClosestCluster = 0;
			my $NextClosestCluster_nED = 1000;
			my $ClusterNumber = 0;
			my $NCC_nED = 1000;
			my $NNCC_nED = 1000;
			my $ClosestCluster_Chr = "";
			my $ClosestCluster_TgtSNPnuc = "";
			if ($_[$Sample] ne "NaN" and $_[$Sample+1] ne "NaN") {
				for ($c=0; $c < @ClusterDetails; $c++) {
					$ClusterNumber++;
					if ($ClusterDetails[$c]->[4] != 0 and $ClusterDetails[$c]->[6] != 0) {
						$nED = sqrt(((($_[$Sample] - $ClusterDetails[$c]->[3])/$ClusterDetails[$c]->[4])**2) + ((($_[$Sample+1] - $ClusterDetails[$c]->[5])/$ClusterDetails[$c]->[6])**2));
					} else {
						$nED = sqrt(($_[$Sample] - $ClusterDetails[$c]->[3])**2 + ($_[$Sample+1] - $ClusterDetails[$c]->[5])**2);
					}
					$NextClosestCluster = $ClosestCluster, $NextClosestCluster_nED = $ClosestCluster_nED, $ClosestCluster = $ClusterNumber, $ClosestCluster_nED = $nED if ($nED < $ClosestCluster_nED);
					$NextClosestCluster = $ClusterNumber, $NextClosestCluster_nED = $nED if ($nED < $NextClosestCluster_nED and $nED > $ClosestCluster_nED);
				}
				foreach $ClusterAssocLoci (@{ $MappedClusterDetails{$ClusterDetails[$ClosestCluster-1]->[9]} }) {
					$AssocLocus = substr($ClusterAssocLoci->[0],0,2);
					$AssocTgtSNPnuc = $ClusterAssocLoci->[1];
					push @{ $ClustersUsed{$ClosestCluster."_".$AssocLocus} }, ( $AssocLocus, $ClusterDetails[$ClosestCluster-1]->[3], $ClusterDetails[$ClosestCluster-1]->[4], $ClusterDetails[$ClosestCluster-1]->[5], $ClusterDetails[$ClosestCluster-1]->[6], $AssocTgtSNPnuc ) unless exists $ClustersUsed{$ClosestCluster."_".$AssocLocus}; 
				}
				if ($ClosestCluster != 0) {
					for (my $i=0; $i < @{ $MappedClusterDetails{$ClusterDetails[$ClosestCluster-1]->[9]} }; $i++) {
						$ClosestCluster_Chr = substr($MappedClusterDetails{$ClusterDetails[$ClosestCluster-1]->[9]}[$i]->[0],0,2), $ClosestCluster_TgtSNPnuc = $MappedClusterDetails{$ClusterDetails[$ClosestCluster-1]->[9]}[$i]->[1] if (substr($MappedClusterDetails{$ClusterDetails[$ClosestCluster-1]->[9]}[$i]->[0],0,2) eq $Locus);
					}
					for ($c=0; $c < @ClusterDetails; $c++) {
						my $SameMappedLocus = 0;
						my $Cluster_TgtSNPnuc = "";
						for (my $i=0; $i < @{ $MappedClusterDetails{$ClusterDetails[$c]->[9]} }; $i++) {
							$SameMappedLocus = 1, $Cluster_TgtSNPnuc = $MappedClusterDetails{$ClusterDetails[$c]->[9]}[$i]->[1] if (substr($MappedClusterDetails{$ClusterDetails[$c]->[9]}[$i]->[0],0,2) eq $Locus);
						}	
						next if ((uc($Cluster_TgtSNPnuc) eq "A" or uc($Cluster_TgtSNPnuc) eq "B") 
								 and $SameMappedLocus == 1 
								 and uc($Cluster_TgtSNPnuc) eq uc($ClosestCluster_TgtSNPnuc)); 
						if ($ClusterDetails[$c]->[4] != 0 and $ClusterDetails[$c]->[6] != 0) {
							$nED = sqrt(((($ClusterDetails[$ClosestCluster-1]->[3] - $ClusterDetails[$c]->[3])/$ClusterDetails[$c]->[4])**2) + ((($ClusterDetails[$ClosestCluster-1]->[5] - $ClusterDetails[$c]->[5])/$ClusterDetails[$c]->[6])**2));
						} else {
							$nED = sqrt(($ClusterDetails[$ClosestCluster-1]->[3] - $ClusterDetails[$c]->[3])**2 + ($ClusterDetails[$ClosestCluster-1]->[5] - $ClusterDetails[$c]->[5])**2);
						}
						$NNCC_nED = $NCC_nED, $NCC_nED = $nED if ($nED < $NCC_nED);
						$NNCC_nED = $nED if ($nED < $NNCC_nED and $nED > $NCC_nED);
					}
				}		
				$ConfidenceScore = (pi/2 - atan($ClosestCluster_nED/$NNCC_nED))/(pi/2);			
			}
			if ($ConfidenceScore >= $ConfidenceCutoff and $ClosestCluster_nED <= sqrt($DistanceLimit**2 + $DistanceLimit**2)) {
				$NumClustered++; 
				print GENOtable "\tC".$ClosestCluster."\t".$ConfidenceScore;  
				if ($ClosestCluster_TgtSNPnuc ne "C" and $ClosestCluster_TgtSNPnuc ne "U" and $ClosestCluster_Chr eq $Locus) {
					print GENOtable "\t".$ClosestCluster_TgtSNPnuc.$ClosestCluster_TgtSNPnuc;
					push @Genotypes, $ClosestCluster_TgtSNPnuc.$ClosestCluster_TgtSNPnuc;
					$GenoCall{uc($ClosestCluster_TgtSNPnuc)}++;
				} else {
					print GENOtable "\tNC";	
					push @Genotypes, "NC";
					$GenoCall{"NC"}++;  
				}
			} else { $GenoCall{"NC"}++; push @Genotypes, "NC"; print GENOtable "\tNC\t0\tNC"; }
		}
		print GENOtable "\n"; 
		print GENOreport "\t".(sprintf '%.1f', ($NumClustered/$SampleSize*100))."\t".(sprintf '%.1f', ($GenoCall{"A"}+$GenoCall{"B"})/$SampleSize*100)."\t".$GenoCall{"A"}."\t".$GenoCall{"B"}."\t".$GenoCall{"NC"}."\n";
		if ($GenoCall{"A"} > 0 and $GenoCall{"B"} > 0) {
			if ($GenoCall{"A"} < $GenoCall{"B"}) { $MinorAlleleFreq = $GenoCall{"A"}/$SampleSize; } else { $MinorAlleleFreq = $GenoCall{"B"}/$SampleSize; }
		}
		if ($Comment eq "Mapped" and $MinorAlleleFreq >= $MAF) {
			if ($GenoCall{"A"}>0 and $GenoCall{"B"}>0) {
				print GENOinfoReport $_[0]."\t".$_[1]."\t".$Locus."\t".$Comment;
				print GENOinfoTable $_[0]."\t".$_[1]."\t".$Locus;
				foreach $g (@Genotypes) { print GENOinfoTable "\t".$g; } print GENOinfoTable "\n";
				print GENOinfoReport "\t".(sprintf '%.1f', ($NumClustered/$SampleSize*100))."\t".(sprintf '%.1f', ($GenoCall{"A"}+$GenoCall{"B"})/$SampleSize*100)."\t".$GenoCall{"A"}."\t".$GenoCall{"B"}."\t".$GenoCall{"NC"}."\n";
			}
		} elsif ($Comment eq "Mapped") {
			print GENOinfoReport $_[0]."\t".$_[1]."\t".$Locus."\t".$Comment;
			print GENOinfoTable $_[0]."\t".$_[1]."\t".$Locus;
			foreach $g (@Genotypes) { print GENOinfoTable "\t".$g; } print GENOinfoTable "\n";
			print GENOinfoReport "\t".(sprintf '%.1f', ($NumClustered/$SampleSize*100))."\t".(sprintf '%.1f', ($GenoCall{"A"}+$GenoCall{"B"})/$SampleSize*100).	"\t".$GenoCall{"A"}."\t".$GenoCall{"B"}."\t".$GenoCall{"NC"}."\n";
		}
	}
	foreach $c (sort keys %ClustersUsed) {
		@HashKeyDetails = split ('_', $c);
		print ClusterPostions $_[0]."\t".$_[1]."\t"."C".$HashKeyDetails[0]."\t"; 
		my @Summary = @{ $ClustersUsed{$c} };
		print ClusterPostions join ("\t", @Summary)."\n";
	}	
}
