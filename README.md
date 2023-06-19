# Concatenation of segmented viral genomes for reassortment analysis
This [app](https://melibrun.shinyapps.io/combined_the_segments/) concatenate segmented viral genomes for reassortment analysis.  
- input format file -> genbank  
- output format file -> fasta  

In order to concatenate segments of your virus, you should upload a **gb** file format to the application. You should also set the required number of segments in the field (Number of segments). Next, you should set an acceptable spread of nucleotides in the field (Permissible difference). After that, you should specify the number of nucleotides in each of the segments. Be careful.  

If you want to concatenate segments of a virus that does not consist of 8 segments, you need to enter the number of your segments in the field (Number of segments) and change the length of the segments according to the subtype of the virus. The remaining segments can be left untouched.  

If you want to have only protein-coding sequences in the final file, switch "Turn on if you want only protein-coding sequence".  

If your virus has a suspiciously small number of sequences, try switching the "Turn on if the name of your virus is contained in the descriptor 'organization'" field. Perhaps the name of your virus is hidden in the "organism" field.



