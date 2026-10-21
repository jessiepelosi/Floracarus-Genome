# Intergenic distances 

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(wesanderson)

file_paths <- list.files(path = "comparitive/intergenic_distances/", pattern = "\\.txt$", full.names = TRUE)

# Read files into a list
data_list <- lapply(file_paths, read.delim, header = F)

# Optional: Name each item in the list by its file name for easy access
names(data_list) <- basename(file_paths)


median(data_list[["bycatch_intergenic_dist.txt"]]$V4)

ggplot() + 
  geom_boxplot(data=data_list[["Floracarus_intergenic_dist.txt"]], mapping = aes(x= V4, y = -1), fill = "lightblue") +
  geom_boxplot(data=data_list[["RARGMITE_intergenic_dist.txt"]], mapping = aes(x= V4, y = 1), fill = "red") +
  geom_boxplot(data=data_list[["Aculus_intergenic_dist.txt"]], mapping = aes(x= V4, y = 0), fill = "forestgreen") +
  geom_boxplot(data=data_list[["Aculops_intergenic_dist.txt"]], mapping = aes(x = V4, y = 2), fill = "orange") +
  xlim(0, 10000) + xlab("Integenic Distance (bp)") + ylab("Species") + theme_bw()

ggsave("intergenic_distances.pdf", width = 7, height = 5)


# Introns

file_paths <- list.files(path = "comparitive/intron_len/", pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(file_paths, read.delim, header = F)

names(data_list) <- basename(file_paths)

ggplot() + 
  geom_boxplot(data=data_list[["Floracarus_intron_len.txt"]], mapping = aes(x= V4, y = -1), fill = "lightblue") +
  geom_boxplot(data=data_list[["RARGMITE_intron_len.txt"]], mapping = aes(x= V4, y = 1), fill = "red") +
  geom_boxplot(data=data_list[["Aculus_intron_len.txt"]], mapping = aes(x= V4, y = 0), fill = "forestgreen") +
  geom_boxplot(data=data_list[["Aculops_intron_len.txt"]], mapping = aes(x = V4, y = 2), fill = "orange") +
  geom_boxplot(data=data_list[["Fragariocoptes_intron_len.txt"]], mapping = aes(x = V4, y = -2), fill = "pink") +
  geom_boxplot(data=data_list[["Setoptus_intron_len.txt"]], mapping = aes(x = V4, y = -3), fill = "purple") +
  xlim(0, 1000) + xlab("Intron Length (bp)") + ylab("Species") + theme_bw()

ggsave("intergenic_distances.pdf", width = 7, height = 5)

# I want to plot densities, so need to merge lists into a single big dataframe

intron_all <- bind_rows(data_list, .id = "sample") |>
  filter(V4 > 1 & sample != "Setoptus_intron_len.txt" & sample != "Fragariocoptes_intron_len.txt") |>
  mutate(sample = factor(sample, levels = rev(c("Aculops_intron_len.txt", "RARGMITE_intron_len.txt", 
                                            "Aculus_intron_len.txt", "Floracarus_intron_len.txt"))))

meds <- intron_all |>
  group_by(sample) |>
  summarise(med = median(V4))

means <- intron_all |>
  group_by(sample) |>
  summarise(mean = mean(V4))

ggplot(intron_all, aes(x = V4, y = sample, fill = sample)) +
  geom_density_ridges(
    scale = 1,              
    rel_min_height = 0.001, 
    color = "black", linewidth = 0.3, alpha = 0.9
  ) +
  geom_segment(data = meds,
               aes(x = med, xend = med, y = as.numeric(sample),
                   yend = as.numeric(sample) + 0.6),
               inherit.aes = FALSE, linewidth = 0.8) +
  scale_fill_manual(values = wes_palette("Darjeeling1", type = "discrete")) +
  labs(x = "Intron length (bp)", y = NULL) + xlim(0, 1000) +
  theme_bw(base_size = 13)

ggsave("intron_lengths_1000.pdf", width = 9, height = 7)



ggplot(intron_all, aes(x = V4, y = sample, fill = sample)) +
  geom_density_ridges(
    scale = 1,              
    rel_min_height = 0.001, 
    color = "black", linewidth = 0.3, alpha = 0.9
  ) +
  geom_segment(data = meds,
               aes(x = med, xend = med, y = as.numeric(sample),
                   yend = as.numeric(sample) + 0.6),
               inherit.aes = FALSE, linewidth = 0.8) +
  scale_fill_manual(values = wes_palette("Darjeeling1", type = "discrete")) +
  labs(x = NULL , y = NULL) + xlim(0, 10000) +
  theme_bw() + theme(legend.position = "none", axis.text.y = element_blank())

ggsave("intron_lengths_10000.pdf", width = 4, height = 4)


