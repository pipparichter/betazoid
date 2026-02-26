I have a situation where I think there are enough reads present in a sample to assemble a genome for an organism, but there is so much strain variation that the scaffolds are getting broken up into pieces with regions of very high identity. I would essentially want a read assembler like a de Brujin graph that allows "jumps" between edges, i.e. infers scaffold order using reads which are fairly similar, but undoubtedly belong to different strains. However, I do not want to create a chimeric genome, i.e. I want to try to capture strain-level variation. I think the thing to do might be to manually look at the assembly graph and see where strain "bubbles" might be forming.

Tool or pipeline for graphical curation
1. Iterative read recruitment. Starting with a seed contig, recruit all reads and partners, and then iteratively extract all reads which map to the initial reads to construct the fully-connected subgraph of the assembly graph. By limiting to a subset (probably would require longer k-mer chunks), this *should* be computationally-tractable.
2. Manually curate resulting graph by constructing a handful of possible contigs and eliminating paths which obviously produce chimeras. (I wonder if some kind of machine learning approach could be applied here, e.g. using a gLM to choose most likely paths through the graph).
3. Visual inspection of the graph, seeing if any bubbles collapse back into single nodes. This could perhaps also be used to see if lots of recombination occurred, in which case there would be lots of branch points that don't obviously recombine.

This approach should help reduce dependence on coverage for grouping contigs that belong to the Betazoids, which seems to be problematic considering some genes are nearly identical and therefore have inflated coverage. 

# Research plan 

1. Create a set of complete or mostly-complete Betazoid genomes from different environmental samples. 
2. Create a set of *Methanoperedens* genomes which co-occur with the Betazoids. 
    *How complete do the genomes need to be in order to reliably place them in a tree? Are one or two marker genes enough, and can I just BLAST against 
    sample contigs instead of binning?* 
3. Characterize similaritites between different types of Betazoids. Is there a common gene set, or common characteristics among the genes which are present? Can we then use this information to broaden the search?* 
4. Are there any common *Methanoperedens* clades which seem to occur with Betazoids?
5. Revisit analysis of SR-VP bioreactor data to see if I can say anything about what they might be doing. 

## Building *Methanoperedens* tree

Before looking through the projects, use the GTDB aligned marker proteins to build a reference tree of *Methanoperedens*. I should download the entire genome for each reference, as it seems like GTDB-Tk typically expects MAGs. Use these genomes to construct a reference tree using GTDB-Tk. 

1. Go through each project in which Betazoids are found and download the tables of Archaeal single-copy genes and ribosomal proteins. These will also include taxonomy, as well as the GC content and coverage for the contig on which the gene is found. 
2. Filter for genes which are assigned Eukaryarchaeota taxonomy.
3. Use the coverage and GC content information to nominally group the contigs into bins. Filter the pseudo-bins according to GC content, with the expectation that *Methanoperedens* GC content is around 40-45%. This should also help reduce contamination from Borgs. 
4. Download all contigs on which SCGs and ribosomal proteins are identified and group them into FASTA files. 
5. Run GTDB-Tk to build a tree from scratch using the reference *Methanoperedens* genomes and the partial bins from the samples. 

```
gtdbtk identify --genome_dir ref/ --out_dir markers/ --extension fasta # Extract GTDB markers from the complete genomes. 
gtdbtk align --identify_dir markers/ --out_dir alignments/ # Align the extracted markers. 
iqtree2 -s alignments/alignments/gtdbtk.ar53.msa.fasta -m LG+G -bb 1000 -nt AUTO # Contstruct a reference tree. 
```
Use the same process as above to construct alignments for the partial bins constructed using contigs from the SCG contigs. 

(https://academic.oup.com/sysbio/article/68/2/365/5079844?login=false )
```
epa-ng --ref-msa ref_align/alignments/gtdbtk.ar53.msa.fasta --tree ref.treefile --query alignments/alignments/gtdbtk.ar53.msa.fasta --model LG+G --outdir out/
```