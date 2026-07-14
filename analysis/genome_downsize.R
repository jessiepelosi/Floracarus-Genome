## genome_downsize.R
##
##



library(readxl)
library(ggplot2)
library(dplyr)

genome_metrics <- read_excel("genome_metrics.xlsx")

# Remove Mesostigmata 

genome_metrics_filt <- genome_metrics %>% 
  filter(Taxon != "Varroa destructor" & Taxon != "Dermanyssus gallinae" & 
           Taxon != "Stratiolaelaps scimitus" & Taxon != "Galendromus occidentalis") %>% 
  filter(N50 > 10000)

ggplot(data = genome_metrics_filt, mapping = aes(x = GenomeSize, y = PropSingleExon )) + 
  geom_point() + theme_bw() + geom_smooth(method = "lm")

ggplot(data = genome_metrics_filt, mapping = aes(y = PropSingleExon, x = N50 )) + 
  geom_point() + theme_bw() 

linear_model <- lm(GenomeSize ~ PropSingleExon, data = genome_metrics_filt)
summary(linear_model)



