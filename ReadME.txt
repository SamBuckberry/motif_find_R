### How to get the data in the format you need

1. Go to UCSC table browser.

2. Select organism

3. Paste identifiers into the "identifiers(names/accessions)"

4. In output format section, select "sequence" option

5. Select genomic option

6. Check box for Promoter/Upstream and add number of bases upstream to search for motif

7. Now you have the .fasta file we need

8. Repeat steps above, but select BED file output

9. Create one bed record per gene

10. Open the .fasta file and remove the "mm10_knownGene_" from header
You can do this using the replace option in gedit.
The header line must read >geneidentifier 

