Running Concoct for the BLUEPRINT project
=========================================


Version 1
---------
Using the reads from the LMO 2012 project. 
The fasta files did not work straight away, it made picard markduplicates fail and complain about duplicated names. 
Tried upgrading bowtie2 and changing parameters but no luck. Sed did the trick, adding a /1 or /2 tag to each read id in the fasta file to indicate which mate it was.


