# gene_track_plots


## Get protein domain info:
- Download data from [UCSC data integrator](https://genome.ucsc.edu/cgi-bin/hgIntegrator)
- genome: Human
- assembly: hg38
- region to annotate: genome
- track group: Genes and Gene Predictions
- track: UniProt
- subtracks: added each individually, downloaded each individually
- (downloaded as .gz into ./data/raw/ucsc/)
- (post download, manually uncomment the column headers for each file)
- NOTE: signal peptides have to be downloaded with the TABLE browser not integrator, for some reason...

You can click “view schema” for each subtrack, e.g., transmembrane domains:
https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=uniprot&hgta_table=unipLocTransMemb&hgta_doSchema=1

Now that I am looking at this data again, I see one of the columns is UniProt ID, so, one could map UniProt ID -> Gene ID...
