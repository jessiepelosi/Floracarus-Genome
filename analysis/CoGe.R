# plotCOGE.R

library(dplyr)
library(ggplot2)
library(ggh4x)

# Get contig lengths: 

Floracarus_fai <- read.delim("Floracarus_perrepae_v1.fa.fai", header = F)
Aculops_fai <- read.delim("comparitive/Aculy_genome.tfa.fai", header = F)
RARGMITE_fai <- read.delim("comparitive/RARGMITEv1.1.fasta.fai", header = F)

## Floracarus - RARGMITE

cds_synmap <- read.delim("comparitive/CoGe_SynMap_CDS_ks.txt", header = F, sep = "")

cds_synmap_filter <- cds_synmap[!grepl('^#',cds_synmap$V1),]

cds_synmap_blocks <- cds_synmap[grepl("^#[0-9]+", cds_synmap$V1),]

mean(as.numeric(cds_synmap_blocks$V6), na.rm = T)

cds_synmap_filter_main <- cds_synmap_filter %>% 
  filter(V3 == "a69908_Floracarus_perrepae_u1402870ctg" | V3 == "a69908_Floracarus_perrepae_u2121237ctg")

cds_synmap_filter_main$V3 <- gsub("a69908_", "", cds_synmap_filter_main$V3)
cds_synmap_filter_main$V7 <- gsub("b69909_", "", cds_synmap_filter_main$V7)

ggplot(data = cds_synmap_filter_main, mapping = aes(x = as.numeric(V5), y = as.numeric(V9))) + geom_point(size = 0.5) +
  facet_grid(V7~V3, scales = "free", space = "free") + theme_bw() + 
  theme(panel.spacing = unit(0, "lines"), axis.text = element_blank(), axis.ticks = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(face = "bold")) +
  facetted_pos_scales(
    x = list(
      V3 == "Floracarus_perrepae_ptg000008l" ~ scale_x_continuous(limits = c(0, 785328)),
      V3 == "Floracarus_perrepae_u13737484ctg" ~ scale_x_continuous(limits = c(0, 1762482)),
      V3 == "Floracarus_perrepae_u1402870ctg" ~ scale_x_continuous(limits = c(0, 10295163)), 
      V3 == "Floracarus_perrepae_u2121237ctg" ~ scale_x_continuous(limits = c(0,10182834))
    ),
    y = list(
      V7 == "RARGMITECHR1" ~ scale_y_continuous(limits = c(0, 17366404)),
      V7 == "RARGMITECHR2" ~ scale_y_continuous(limits = c(0,16977790))
    )
  ) +
  xlab(expression(paste(italic("Floracarus perrepae"), " Genomic Position"))) + 
  ylab("RARGMITE Genomic Position")

ggsave("F-R_dotplot.pdf", height = 6, width = 6)


## Floracarus - Aculops
FA_cds_synmap <- read.delim("comparitive/Aculops_Floracarus_CoGe_SynMap_CDS_ks.txt", header = F, sep = "")

FA_cds_synmap_filter <- FA_cds_synmap[!grepl('^#',FA_cds_synmap$V1),]

FA_cds_synmap_blocks <- FA_cds_synmap[grepl("^#[0-9]+", FA_cds_synmap$V1),]

mean(as.numeric(FA_cds_synmap_blocks$V6), na.rm = T)

FA_cds_synmap_filter_main <- FA_cds_synmap_filter %>% 
  filter(V1 == "a69908_Floracarus_perrepae_u1402870ctg" | V1 == "a69908_Floracarus_perrepae_u2121237ctg")

FA_cds_synmap_filter_main$V1 <- gsub("a69908_", "", FA_cds_synmap_filter_main$V1)
FA_cds_synmap_filter_main$V5 <- gsub("b69911_", "", FA_cds_synmap_filter_main$V5)

ggplot(data = FA_cds_synmap_filter_main, mapping = aes(x = as.numeric(V3), y = as.numeric(V7))) + geom_point(size = 0.5) +
  facet_grid(V5~V1, scales = "free", space = "free") + theme_bw() + 
  theme(panel.spacing = unit(0, "lines"), axis.text = element_blank(), axis.ticks = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(face = "bold")) +
  facetted_pos_scales(
    x = list(
      V1 == "Floracarus_perrepae_ptg000008l" ~ scale_x_continuous(limits = c(0, 785328)),
      V1 == "Floracarus_perrepae_u13737484ctg" ~ scale_x_continuous(limits = c(0, 1762482)),
      V1 == "Floracarus_perrepae_u1402870ctg" ~ scale_x_continuous(limits = c(0, 10295163)), 
      V1 == "Floracarus_perrepae_u2121237ctg" ~ scale_x_continuous(limits = c(0,10182834))
    ),
    y = list(
      V5 == "scaffold00001" ~ scale_y_continuous(limits = c(0, 12438527)),
      V5 == "scaffold00002" ~ scale_y_continuous(limits = c(0,10499462)),
      V5 == "scaffold00003" ~ scale_y_continuous(limits = c(0,3663898)),
      V5 == "scaffold00004" ~ scale_y_continuous(limits = c(0,3565256)),
      V5 == "scaffold00005" ~ scale_y_continuous(limits = c(0,2355566))
    )
  ) +
  xlab(expression(paste(italic("Floracarus perrepae"), " Genomic Position"))) + 
  ylab(expression(paste(italic("Aculops lycopersici"), " Genomic Position"))) 

ggsave("F-A_dotplot.pdf", height = 12, width = 12)
ggsave("F-A_dotplot.png", height = 12, width = 12)



## Aculops - RARGMITE

AR_cds_synmap <- read.delim("comparitive/Aculops_RARGMITE_CoGe_SynMap_CDS_ks.txt", header = F, sep = "")

AR_cds_synmap_filter <- AR_cds_synmap[!grepl('^#',AR_cds_synmap$V1),]

AR_cds_synmap_blocks <- AR_cds_synmap[grepl("^#[0-9]+", AR_cds_synmap$V1),]

mean(as.numeric(AR_cds_synmap_blocks$V6), na.rm = T)

#AR_cds_synmap_filter_main <- AR_cds_synmap_filter %>% 
#  filter(V3 == "a69908_Floracarus_perrepae_u1402870ctg" | V3 == "a69908_Floracarus_perrepae_u2121237ctg")

ggplot(data = AR_cds_synmap_filter, mapping = aes(x = as.numeric(V3), y = as.numeric(V7))) + geom_point(size = 0.5) +
  facet_grid(V5~V1, scales = "free") + theme_bw() + xlab("RARGMITE Genomic Position (bp)") + ylab("Aculops lycopersici Genomic Position (bp)")

AR_cds_synmap_filter$V1 <- gsub("a69909_", "", AR_cds_synmap_filter$V1)
AR_cds_synmap_filter$V5 <- gsub("b69911_", "", AR_cds_synmap_filter$V5)

ggplot(data = AR_cds_synmap_filter, mapping = aes(x = as.numeric(V3), y = as.numeric(V7))) + geom_point(size = 0.5) +
  facet_grid(V5~V1, scales = "free", space = "free") + theme_bw() + 
  theme(panel.spacing = unit(0, "lines"), axis.text = element_blank(), axis.ticks = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(face = "bold")) +
  facetted_pos_scales(
    x = list(
      V1 == "RARGMITECHR1" ~ scale_x_continuous(limits = c(0, 17366404)),
      V1 == "RARGMITECHR2" ~ scale_x_continuous(limits = c(0,16977790))
    ),
    y = list(
      V5 == "scaffold00001" ~ scale_y_continuous(limits = c(0, 12438527)),
      V5 == "scaffold00002" ~ scale_y_continuous(limits = c(0,10499462)),
      V5 == "scaffold00003" ~ scale_y_continuous(limits = c(0,3663898)),
      V5 == "scaffold00004" ~ scale_y_continuous(limits = c(0,3565256)),
      V5 == "scaffold00005" ~ scale_y_continuous(limits = c(0,2355566))
    )
  ) +
  xlab("RARGMITE Genomic Position") + 
  ylab(expression(paste(italic("Aculops lycopersici"), " Genomic Position"))) 

ggsave("A-R_dotplot.pdf", height = 12, width = 12)
ggsave("A-R_dotplot.png", height = 12, width = 12, dpi = 300)


# Are any homeobox genes in these syntenic regions? 

FP_hox <- read.delim("comparitive/Fp_hox.genes", header =F)

grepl(FP_hox, cds_synmap_filter)
