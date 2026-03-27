**heatmaps_from_m8.py**: creates heatmaps based on mmseqs2 similarity values between two group of proteins homologues. It dereplicates first. then compares group 1 vs group 1, then  group 2 vs group 2,
and finally  group 1 vs group 2. Runs as follow
python heatmaps_from_m8.py \
  --m8 all_vs_all.m8 \
  --prefix-a ProteinA_ \
  --prefix-b ProteinB_ \
  --label-a ProteinA \
  --label-b ProteinB \
  --outdir heatmaps_proteinAvsProtein \
  --derep-id 0.95 \
  --min-id 0.30 \
  --min-aln 50 \
  --max-evalue 1e-3

  Prerequirements:
  1.Get a group of homologues for your protein A and B using blastp or similar. Add a prefix to each homologue in the fasta header (i.e. >ProteinA_WPxxx)
  2. Run Mmseqs2 as follows:
  mmseqs easy-search all.faa all.faa all_vs_all.m8 tmp --format-output "query,target,fident,alnlen,evalue,bits"
**heatmap-neighbourhoods.py**: runs as :-i all.xlsx --mode zscore --n-left 120 --n-right 80
*Check input files folder for formatting of all.xlsx input file*

**AthenianV2.1.sh**:Better version of https://github.com/GM110Z/Phage-defence-scripts/blob/main/Defence-islands-finder/ATHENIAN.sh
Retrieves the neighbourhood of a specified protein Refseq ID list. Runs as:
 python athenianv2.py \                                                                
    -i input.athenian.txt \
    -k 25 \ #this is the Kb you want to allow each side of seed protein
    -t 8 #threads
