# Description of the VSC

[VSC](https://melibrun.shinyapps.io/combined_the_segments/)  enables the automatic concatenation of sequences written in *GenBank* or *FASTA* format. Please note that for correct analysis of FASTA sequences, the header should be written as follows: “>((strain name(subtype))_segment_segment-number” (e.g. ">(A/waterfowl/Korea/S005/2014(H5N5))_segment_3" ). Sample files in Genbank and FASTA format can be found at https://github.com/melibrun/Concatenation-of-segmented-viral-genomes-for-reassortment-analysis/tree/main/examples.   

To visualise the separation between the segments, an option has also been added to add a certain number of “N” symbols.  

The default settings of the application are set for the influenza A viruses, whose genome consists of eight segments of the specified length. If the virus of interest has a different number of segments with a different length, this should be specified in the corresponding fields.  

The “Permissible difference” fields indicate the range of segment length that will be treated as a complete segment sequence. For example, if “Permissible difference” is 200 and the length of segment 1 is 2300, sequences with a length between 2100 and 2500 will be treated as a complete segment length.  

Using the annotated genebank file allows the use of additional options. For example, only protein-coding regions can be concatenated if the corresponding button has been selected.  


