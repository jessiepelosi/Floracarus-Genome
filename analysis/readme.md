# Analysis of Mite Genomes

## Synteny 

We used SynMap2 within CoGe to identify syntenic gene pairs between Eriophyoidae mite genomes (Aculops, Floracarus, and RARGMITE). Using the output from SynMap2, we calculated the number of syntenic blocks, average number of genes per syntenic block, and the number of syntenic gene pairs for each comparison. 

## Phylogenomics of Genome Downsizing 

We used OrthoFinder v3.1.5 to identify low-copy nuclear orthologs among 24 taxa. 
We used only the longest isoform per gene identified with AGAT v1.4.2 as input.

```
```

Single-copy orthogroups present in at least 75% of taxa were aligned with MACSE v2.07 and sites with less than 50% taxon occupancy were removed with trimAl v1.3.
```
```

We then estimated the best substitution model with ModelFinderPlus and built maximum-likelihood phylogenies with IQTREE v3.0.1 with 1000 ultrafast bootstraps. These gene trees were input to ASTRAL IV to generate a species tree under the multi-species coalescent.
```
```

To first determine whether there was phylogenetic signal in genome size across these mites, we used the ‘phylosig’ function from phytools v2.5-2 to estimate and assess the significance of Pagel’s λ with 5000 simulations for the randomization test and 100 iteration for likelihood optimization of log(genome size) with the species tree. See R script. 

We then used a custom python script (orthofinder_to_malin.py) to parse the output from Orthofinder. Briefly, this script edits the multiple sequence alignment headers from OrthoFinder and the input GFFs to the format required for Malin. Malin (Csurös 2008) was used to analyze intron evolution across 21,491 orthogroups. We used the default parameters to generate an initial model of rate evolution (e.g., gains and losses of introns), allowing each branch to have its own rate, and from these initial estimates we optimized the likelihood-based model without rate variation. Then we examined the shared presence of introns and examined gain-loss across the phylogeny with Dollo parsimony and posterior probability. Finally, we generated uncertainty estimates around each of the inferred rate parameters using 500 bootstrap replicates. 
