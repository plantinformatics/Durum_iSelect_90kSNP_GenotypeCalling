# Durum_iSelect_90kSNP_GenotypeCalling
Perl script for generating genotype calls for unrelated durum wheat samples assayed using the Illumina iSelect 90K SNP bead chip array

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
