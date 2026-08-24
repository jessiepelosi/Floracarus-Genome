## phylo_signal.R 

library(phytools)
library(readxl)
library(ggtree)
library(ggplot2)
library(dplyr)
library(picante)
library(geiger)
library(tibble)
library(viridisLite)
library(caper)

### Import tree

spTree <- read.tree("astral_sp.nwk")
plotTree(spTree)
# For now, let's prune Oppia b/c the repeat modeling is being problematic 

spTree.prune <- drop.tip(spTree, tip = "Oppia_nitens")

### Import repeat data 

repeats <- read_excel("Repeats.xlsx")

# Clean up repeat data 
repeats$Length <- gsub(" bp","", repeats$Length)
repeats$Length <- as.numeric(repeats$Length)
repeats$`Repeat Subfamily` <- gsub(":","", repeats$`Repeat Subfamily`)

# Replace space in taxon name with underscore to make them match tree tip names 
repeats$Tasxon <- repeats$Taxon <- gsub(" ", "_", repeats$Taxon)

# Filter out Oppia 
repeats <- filter(repeats, Taxon != "Oppia_nitens")

# Write a function that will run the phylogenetic signal tests: 

filter_run_phylosig <- function(dataframe, repeatfamily, repeatsubfamily, phylogeny) {
  df <- filter(dataframe, `Repeat Family` == repeatfamily & `Repeat Subfamily` == repeatsubfamily) 
  rownames(df) <- df$Taxon
  df_data <- setNames(df$PercGenome, rownames(df))
  contMap(phylogeny, df_data, mode=c("edges","tips","nodes"), palette="rainbow")
  lambda<- phylosig(phylogeny, df_data, method="lambda", test=TRUE, nsim=5000,
           se=NULL, start=NULL, control=list(), niter=10)
  return(lambda)
}

## 1. Total percent repetitive
all_sum_df <- repeats %>% 
  group_by(Taxon) %>% 
  summarize(total = sum(PercGenome, na.rm = T))

rownames(all_sum_df) <- all_sum_df$Taxon

all_sum <- setNames(all_sum_df$total, rownames(all_sum_df))

contMap(spTree.prune, all_sum, mode=c("edges","tips","nodes"), palette="rainbow")

phylosig(spTree.prune, all_sum, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)

# Retrotransposons:

retros_all <- filter_run_phylosig(repeats, "Retroelements", "All", spTree.prune)
retros_sine <- filter_run_phylosig(repeats, "Retroelements", "SINEs", spTree.prune)
retros_lines <- filter_run_phylosig(repeats, "Retroelements", "LINEs", spTree.prune)
retros_ltrs <- filter_run_phylosig(repeats, "Retroelements", "LTR elements", spTree.prune)
retros_copia <- filter_run_phylosig(repeats, "Retroelements", "Ty1/Copia", spTree.prune)
retros_gypsy <- filter_run_phylosig(repeats, "Retroelements", "Gypsy/DIRS1", spTree.prune)

# DNA Transposons: 
dna_all <- filter_run_phylosig(repeats, "DNA transposons", "All", spTree.prune)
dna_Tc1 <- filter_run_phylosig(repeats, "DNA transposons", "Tc1-IS630-Pogo", spTree.prune) # Sig
dna_hobo <- filter_run_phylosig(repeats, "DNA transposons", "hobo-Activator", spTree.prune)
dna_mule <- filter_run_phylosig(repeats, "DNA transposons", "MULE-MuDR", spTree.prune) # Sig
dna_tour <- filter_run_phylosig(repeats, "DNA transposons", "Tourist/Harbinger", spTree.prune)

# Simple Repeats 
df <- filter(repeats, `Repeat Family` == "Simple repeats:")  # sig
rownames(df) <- df$Taxon
df_data <- setNames(df$PercGenome, rownames(df))
n_col <- n_distinct(df_data)
contMap(spTree.prune, df_data, mode=c("edges","tips","nodes"), palette="rainbow")
lambda<- phylosig(spTree.prune, df_data, method="lambda", test=TRUE, nsim=5000,
                  se=NULL, start=NULL, control=list(), niter=10)

# Adjust for multiple comparisons 
p.adjust(p = c(retros_all$P,retros_sine$P, retros_lines$P, retros_ltrs$P, 
         retros_copia$P, retros_gypsy$P, dna_all$P, dna_Tc1$P,
         dna_hobo$P, dna_mule$P, dna_tour$P), method = "BH" )




#### Other genome metrics (structural gene annotations)

genome_metrics <- as.data.frame(read_excel("genome_metrics.xlsx"))
genome_metrics$Taxon <- genome_metrics$Taxon <- gsub(" ", "_", genome_metrics$Taxon)

# Remove Oppia -- RepeatModeler Failed
genome_metrics.prune <- filter(genome_metrics, Taxon != "Oppia_nitens")

rownames(genome_metrics.prune) <- genome_metrics.prune$Taxon

# Explore the data a bit: 

ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = TotalCDSLength, color = Eriophyoidae)) + 
  geom_point(size = 2) + theme_bw() + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size (bp))") + ylab("Total Length of Coding Sequences (bp)") +
  theme(legend.position = "none")


ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = TotalRepeatContentProp, color = Eriophyoidae)) + 
  geom_point(size = 2) + theme_bw() + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size (bp))") + ylab("Total Proportion of Genome Occupied by Repetitive Elements") +
  theme(legend.position = "none")


ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = log(TotalRetroelementProp), color = Eriophyoidae)) + 
  geom_point(size = 2) + theme_bw(base_size = 16) + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size)") + ylab("log(Prop. of Genome Occupied by Retroelements)") +
  theme(legend.position = "none")

ggsave("gs_retros.pdf", height = 6, width = 6)

ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = log(LINEProp), color = Eriophyoidae)) + 
  geom_point(size = 2) + theme_bw(base_size = 16) + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size))") + ylab("log(Prop. of Genome Occupied by LINEs)") +
  theme(legend.position = "none")

ggsave("gs_LINES.pdf", height = 6, width = 6)


ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = log(LTRProp), color = Eriophyoidae)) + 
  geom_point(size =2) + theme_bw() + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size))") + ylab("log(Total Proportion of Genome Occupied by LTRs)") +
  theme(legend.position = "none")


ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = log(TotalDNAelementProp), color = Eriophyoidae)) + 
  geom_point(size =2) + theme_bw() + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size))") + ylab("log(Total Proportion of Genome Occupied by DNA Transposable Elements)") +
  theme(legend.position = "none")


ggplot(data = genome_metrics.prune, mapping = aes(x = log(GenomeSize), y = log(PropSingleExon), color = Eriophyoidae)) + 
  geom_point(size = 2) + theme_bw() + scale_color_manual(values = c("Y" = "limegreen", "N" = "black")) +
  xlab("log(Genome Size))") + ylab("log(Proportion of Single-Exon Genes in Genome)")



# Reorder dataframe to match species tree order:

genome_metrics.prune <- genome_metrics.prune[match(spTree.prune$tip.label, genome_metrics.prune$Taxon),]

# Calculate phylogenetic signal for log(genome size) across the phylogeny

pdf(file = "gs_phylosig.pdf", width = 8, height = 8)

gs <- dplyr::select(genome_metrics.prune, Taxon, GenomeSize)
rownames(gs) <- gs$Taxon
gs <- setNames(gs$GenomeSize, rownames(gs))

cm <- contMap(spTree.prune, log(gs), mode=c("edges","tips","nodes") , plot = F)

cm$tree <- ladderize.simmap(cm$tree, right = F)

cm <- setMap(cm, hcl.colors(n = 50, palette = "Spectral", rev = T))
plot(cm, leg.txt = "log(Genome Size)", offset = 1)

phylosig(spTree.prune, log(gs), method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)


dev.off()

# Proportion of single-exon genes: 

pdf(file = "intron_phylosig.pdf", width = 8, height = 8)

prop_intron <- dplyr::select(genome_metrics.prune, Taxon, PropSingleExon)
rownames(prop_intron) <- prop_intron$Taxon
prop_intron <- setNames(prop_intron$PropSingleExon, rownames(prop_intron))
cm_intron <- contMap(spTree.prune, prop_intron, mode=c("edges","tips","nodes"), plot = F)
cm_intron <- setMap(cm_intron, hcl.colors(n = 50, palette = "Spectral", rev = T))
plot(cm_intron, leg.txt = "Prop. Intronless Genes", offset = 1)

phylosig(spTree.prune, prop_intron, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)

dev.off()

# Total repeat content proportion: 

all_reps <- dplyr::select(genome_metrics.prune, Taxon, TotalRepeatContentProp)
rownames(all_reps) <- all_reps$Taxon
all_reps <- setNames(all_reps$TotalRepeatContentProp, rownames(all_reps))
contMap(spTree.prune, all_reps, mode=c("edges","tips","nodes"), palette="rainbow")
phylosig(spTree.prune, all_reps, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)

retros <- dplyr::select(genome_metrics.prune, Taxon, TotalRetroelementProp)
rownames(retros) <- retros$Taxon
retros <- setNames(retros$TotalRetroelementProp, rownames(all_reps))
contMap(spTree.prune, retros, mode=c("edges","tips","nodes"), palette="rainbow")
phylosig(spTree.prune, retros, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)

DNA <- dplyr::select(genome_metrics.prune, Taxon, TotalDNAelementProp)
rownames(DNA) <- DNA$Taxon
DNA <- setNames(DNA$TotalDNAelementProp, rownames(DNA))
contMap(spTree.prune, DNA, mode=c("edges","tips","nodes"), palette="rainbow")
phylosig(spTree.prune, DNA, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)


# CDS 
prop_cds <- genome_metrics.prune %>% 
  mutate(prop_cds = TotalCDSLength/GenomeSize)
prop_cds <- dplyr::select(prop_cds, Taxon, prop_cds)
rownames(prop_cds) <- prop_cds$Taxon
prop_cds <- setNames(prop_cds$prop_cds, rownames(prop_cds))
cm_cds <- contMap(spTree.prune, prop_cds, mode=c("edges","tips","nodes"), plot = F)
cm_cds <- setMap(cm_cds, hcl.colors(n = 50, palette = "Spectral", rev = T))
plot(cm_cds, leg.txt = "Total Coding Sequence Legnth", offset = 1)

phylosig(spTree.prune, prop_cds, method="lambda", test=TRUE, nsim=5000,
         se=NULL, start=NULL, control=list(), niter=100)



## PGLS



# caper: 
comp.data <- comparative.data(spTree, genome_metrics, names.col = "Taxon")

model <- pgls(TotalCDSLength ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model)

model <- pgls(TotalGeneLength ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model)

model_caper <- pgls(PropSingleExon ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model_caper)

model_caper <- pgls(TotalRetroelementProp ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model_caper)

model_caper <- pgls(LINEProp ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model_caper)

model_caper <- pgls(LTRProp ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model_caper)

model_caper <- pgls(TotalDNAelementProp ~log(GenomeSize), data = comp.data, lambda = "ML")
summary(model_caper)
