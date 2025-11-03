Scripts and modules for creating a DNA binding site prediction model. The model
will take as input the sequence and output a predicted DNA binding sequence.

#############################################################################

#############################################################################


bsp_utils.py:

    contains useful functions/classes 

###########################################################################
analyzeSiteData.py:

analyze a data file containing cognate site data.
format should be two columns separated by a delimiter
    RE Name (1)     site sequence (2)
    
output:
    number of entries
    site length stats
    histogram of site lengths
    histogram of character use in site strings
    
###########################################################################
analyzeSequenceData.py:

analyze a data file containing protein sequences and possible sites
format should be FASTA
    
output for (1) all sequences, and (2) all sequences with site data:
    statistics on lengths
    histogram of lengths (only for (2))
    
option for saving data in csv format as either 

    (1) site data file:   2 columns:  RE    SITE
    (2) combined data:    3 columns:  RE    SITE    SEQUENCE
    
###########################################################################
combineSiteFiles.py:

read in a list of site data files and create a list of all unique RE names
that occur at least once in all the files. Data files should all be in same
directory. site data file format should be two columns separated by a 
delimiter:

    RE Name (1)     site sequence (2)

will create a data frame and save if requested that contains :

1) the unique RE name 
2) the site sequences 
3) list of files contained in
 
for now, will remove cut site symbol from strings, so character list is 
assumed to be 

    'ACGTNBDHKMRSUVWY'

Will need to remove non-recognized characters and convert to upper case

note: should check for consistency. If not, what? perhaps save in a separate 
dateframe that records all sequences and file names for inspection

finally, will print out:
    number of total entries
    number of unique entries
    number of unique entries with inconsistent site sequences
    site length stats
    histogram of site lengths
    histogram of character use in site strings
    
notes: 
    
    1) format 15 has multiple entries for promiscuous REs, eschewing use of
    ambiguity code

    2) format 2 has flanking Ns and lists reverse complement for some
    
    3) format 13 & 15 has reverse complements and leading/trailing Ns
    
    4) format 33 has leading spaces

    5) formats 17 & 23 appear to have some entries reverse complement of 
    each other
    
    6) true inconsistencies seem limited to lpnpi,hpyum037x
    
    7) formats 17 and 23 appear to have reverse complements of each
    other for some entries
    
###########################################################################
createData.py:

combine a site data file and a sequence data file. site data file must be in csv format with one column 'RE' and the second one 'site'
sequence data must be in fasta format with name of enzyme leading in title 
line. 

for example:

>AatII   GACGTC  345 aa

Code will find intersection of RE names and create combined file with only
those also prints out stats on resulting dataframe
 
########################################################################### 
ORF_P_analysis.py:



###########################################################################

