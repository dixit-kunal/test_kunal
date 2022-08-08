library(microbiome)
library(dplyr)
setwd("C:\\Users\\ncmr\\Desktop\\Kunal_D\\NCGS\\Run_wise_analysis\\Phyloseq\\Analysis_Feb2022\\merged\\plots\\core_microbiome")


ps <- readRDS("merged.ps.rds")
samdf <- read.table(file = "merged.metadata.txt", sep = "\t", header = T, row.names = 1, check.names = F)
row.names(samdf)
sample_data(ps) <- samdf
ps

sample_names(ps)
ps.baseline <- subset_samples(physeq = ps, Time_point == 'Baseline')
ps.pg <- subset_samples(physeq = ps, Time_point == 'Post_GFD')
ps.stool <- subset_samples(physeq = ps.baseline, Sample_Type == 'Whole_gut')
ps.duo <- subset_samples(physeq = ps.baseline, Sample_Type == 'Duodenal')
ps.sig <- subset_samples(physeq = ps.baseline, Sample_Type == 'Sigmoid')

ps.ibs.stool.all <- subset_samples(physeq = ps.stool, Sample_Category == 'IBS(AGA_Negative)')
ps.ibs.duo <- subset_samples(physeq = ps.duo, Sample_Category == 'IBS(AGA_Negative)')
ps.ibs.sig <- subset_samples(physeq = ps.sig, Samples == 'IBS(AGA_Negative)')

ps.healthy <- subset_samples(physeq = ps.stool, Sample_Category == 'Healthy_Control')

ps.ncgs <- subset_samples(physeq = ps.baseline, Sample_Category == 'NCGS')
ps.ncgs.pg <- subset_samples(physeq = ps.pg, Sample_Category == 'NCGS')
ps.ncgs.stool.pg <- subset_samples(physeq = ps.ncgs.pg, Sample_Type == 'Whole_gut')
ps.ncgs.duo.pg <- subset_samples(physeq = ps.ncgs.pg, Sample_Type == 'Duodenal')
ps.ncgs.sig.pg <- subset_samples(physeq = ps.ncgs.pg, Sample_Type == 'Sigmoid')


#################### aggregate_rare - function ##############
aggregate_rare <- function (x, level, detection, prevalence, include.lowest=FALSE, ...) {
  
  x <- aggregate_taxa(x, level)
  
  rare <- rare_members(x, detection, prevalence, include.lowest)
  
  tax <- tax_table(x)
  
  inds <- which(rownames(tax) %in% rare)
  
  tax[inds, level] <- "Other"
  
  tax_table(x) <- tax
  
  tt <- tax_table(x)[, level]
  
  tax_table(x) <- tax_table(tt)
  
  aggregate_taxa(x, level)
  
}

############## core microbiota analysis ##############
ps.rel <- microbiome::transform(ps.healthy, "compositional")
taxa_names(ps.rel)[1:3]
library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(ps.rel))
names(dna) <- taxa_names(ps.rel)
ps.rel <- merge_phyloseq(ps.rel, dna)
taxa_names(ps.rel) <- paste0("ASV", seq(ntaxa(ps.rel)))
# now check again
taxa_names(ps.rel)[1:3]
core.taxa.standard <- core_members(ps.rel, detection = 0.0001, prevalence = 50/100)
ps.core <- core(ps.rel, detection = 0.0001, prevalence = .5)
taxa_names(ps.core)
core.taxa <- taxa(ps.core)
core.taxa
class(core.taxa)
tax.mat <- tax_table(ps.core)
tax.df <- as.matrix(tax.mat)
tax.df <- as.data.frame(tax.df)
tax.df$OTU <- rownames(tax.df)
tax.data <- as.data.frame(tax.df)
core.taxa.class <- dplyr::filter(tax.data, rownames(tax.data) %in% core.taxa)
knitr::kable(head(core.taxa.class))

####################### core line plot ##################
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(ps.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

#################### core heatmaps ######################
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
library(viridis)
library(RColorBrewer)
p1 <- plot_core(ps.rel, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance)") + theme_bw()+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.title = element_text(size = 12, face = "bold"),plot.title = element_text(color="black", size=14, face="bold"),strip.text.x = element_text(color = "black", size = 12, face="bold"), legend.title = element_text(color = "black", size = 12, face="bold")) + ggtitle("Core Microbiome - ASVs relative abundance")
p1 <- p1 + ylab("ASVs")+ scale_fill_viridis()
p1


##################### genus - core heatmap ##################
prevalences <- seq(.05, 1, .05)

detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
ps.genus <- aggregate_taxa(ps.rel, "Genus")
taxa_names(ps.genus)
ps.pruned.genus <- prune_taxa(taxa = taxa_names(ps.genus) != 'Unknown', ps.genus)
p1 <- plot_core(ps.pruned.genus, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = 0.67) +
  xlab("Relative Abundance") + ggtitle("Core Microbiome - Genus relative abundance")+ theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), 
                                                                                            axis.title.y = element_text(size = 16, face = "bold"), 
                                                                                            axis.title.x = element_text(size = 16, face = "bold"),
                                                                                            legend.title = element_text(size = 16,  face = "bold",colour = "black",), 
                                                                                            legend.text = element_text(size = 16, face = "bold", colour = "black"),plot.title = element_text(color="black", size=14, face="bold"),
                                                                                            axis.text.y = element_text(size = 12, face = "bold",colour = "black"))
p1 <- p1 + ylab("Genus")+ scale_fill_viridis(alpha = 1,
                                             begin = 0,
                                             end = 1,
                                             direction = 1,
                                             discrete = FALSE,
                                             option = "H")
p1

write.table(x = p1$data,file = 'ibs_stool_core.txt',quote = F,sep = '\t',col.names = NA)

##################### family - core heatmap ################
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)
ps.phylum <- aggregate_taxa(ps.rel, "Phylum")
ps.pruned.phylum <- prune_taxa(taxa = taxa_names(ps.phylum) != 'Unknown', ps.phylum)

p1 <- plot_core(ps.pruned.phylum, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = 0.67) +
  xlab("Relative Abundance") + ggtitle("Core Microbiome - Phylum relative abundance")+ theme(axis.text.x = element_text(angle = 90, size = 12, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), 
                                                                                            axis.title.y = element_text(size = 16, face = "bold"), 
                                                                                            axis.title.x = element_text(size = 16, face = "bold"),
                                                                                            legend.title = element_text(size = 16,  face = "bold",colour = "black",), 
                                                                                            legend.text = element_text(size = 16, face = "bold", colour = "black"),plot.title = element_text(color="black", size=14, face="bold"),
                                                                                            axis.text.y = element_text(size = 12, face = "bold",colour = "black"))

p1 <- p1 + ylab("Genus")+ scale_fill_viridis(alpha = 1,
                                             begin = 0,
                                             end = 1,
                                             direction = 1,
                                             discrete = FALSE,
                                             option = "H")
p1

write.table(x = p1$data,file = 'healthy_stool_phylum_core.txt',quote = F,sep = '\t',col.names = NA)
