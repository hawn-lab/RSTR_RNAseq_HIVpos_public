library(tidyverse)
library(ggvenn)
library(ggrepel)
library(ggpubr)
library(patchwork)

#### Main model data ####
attach("results/HIVpos_lme.RData")

HIV_pos_results <- model_intr$lme %>% 
  #Select main model variables
  filter(variable != "(1 | FULLIDNO)") %>% 
  #remove hgnc list formatting
  unnest(symbol) %>% 
  select(model, variable, gene, symbol, estimate, FDR) %>% 
  #rename variables
  mutate(var_name = case_when(variable == "condition" ~ "+Mtb vs media        ",
                              variable == "Sample_Group" ~ "        RSTR vs LTBI",
                              variable == "condition:Sample_Group" ~ "Mtb:RSTR interaction"),
         var_name = factor(var_name, 
                           levels = c("+Mtb vs media        ", "        RSTR vs LTBI",
                                      "Mtb:RSTR interaction"))) %>% 
  arrange(var_name)

#### Venn FDR < 0.2 ####
fdr.cutoff <- 0.2
venn_dat <- list()

for (v in unique(HIV_pos_results$var_name)) {
  venn_dat[[v]] <- HIV_pos_results %>% 
    filter(FDR < fdr.cutoff) %>% 
    filter(var_name == v) %>% 
    pull(gene) %>% unique()
}

p1 <- ggvenn(venn_dat, 
             show_percentage = FALSE, text_size = 4, set_name_size = 4, 
             stroke_size = 0.5, fill_color = rep("white", length(venn_dat)))
# p1

#### Volcano ####
dat_volc <- HIV_pos_results %>% 
  #select variables of interest
  filter(var_name %in% c("   RSTR vs LTBI",
                         "Mtb:RSTR interaction")) %>% 
  #Remove extra spaces in variable levels
  mutate(var_name = recode(var_name, "   RSTR vs LTBI"="RSTR vs LTBI")) %>% 
  #Set significance color groups
  mutate(label = ifelse(FDR < 0.2, symbol, NA),
         col.group = dplyr::case_when(
           FDR < fdr.cutoff & estimate < 0 ~ "down", 
           FDR < fdr.cutoff & estimate > 0 ~ "up", 
           TRUE ~ "NS")) 

#plot
p2 <- dat_volc %>% 
  ggplot(aes(x = estimate, y = -log10(FDR))) + 
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) + 
  labs(y = "-log10( FDR )", 
       x = "Log2 fold change") + 
  # geom_hline(yintercept = -log10(fdr.cutoff), lty = "dashed", color = "grey90") +
  facet_wrap(~var_name) +
  #NS data points 
  ##Need to be under DEG labels
  geom_point(data=filter(dat_volc, col.group == "NS"),
             aes(color = col.group, shape = col.group)) + 
  #label DEG
  ##Up
  geom_text_repel(data=filter(dat_volc, col.group == "up"), 
                  aes(label = label), size=3,
                  nudge_x = 6.2 - filter(dat_volc, col.group == "up")$estimate, 
                  segment.size = 0.2, angle = 0, hjust = 1, box.padding = 0.5,
                  segment.color = "black", direction = "y",
                  show.legend = FALSE, max.overlaps = Inf, seed=32) +
  ##down
  geom_text_repel(data=filter(dat_volc, col.group == "down"),
                  aes(label = label), size=3,
                  nudge_x = -3.8 - filter(dat_volc, col.group == "down")$estimate,
                  segment.size = 0.2, angle = 0, hjust = 0, box.padding = 0.5,
                  segment.color = "black", direction = "y",
                  show.legend = FALSE, max.overlaps = Inf, seed=32) +
  geom_vline(xintercept = 0, color = "black") +
  #Signif data points 
  ##Need to be on top of DEG labels
  geom_point(data=filter(dat_volc, col.group != "NS"),
             aes(color = col.group, shape = col.group)) + 
  #color and shape
  scale_color_manual(name = "FDR < 0.2",
                     values = c(down = "#005AB5", 
                                NS = "grey90", 
                                up = "#DC3220"), 
                     na.value = "grey90") + 
  scale_shape_manual(name = "FDR < 0.2",
                     values = c(down = "triangle", 
                                NS = "circle", 
                                up = "square")) +
  #x-axis
  scale_x_continuous(limits = c(-4,6.3),
                     breaks = seq(-4,4,2))

# p2

#### Expression data ####
genes.OI <- c("ITGB8", "CDC42BPG")

#Check all genes of interest are DEG
genes.OI %in% unlist(model_con$lme.contrast$symbol)

attach("data_clean/HIV_pos_voom.RData")
dat <- as.data.frame(dat.voom$E) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  pivot_longer(-ensembl_gene_id, names_to = "libID") %>% 
  left_join(dat.voom$targets) %>% 
  left_join(dat.voom$genes) %>% 
  unnest(symbol) %>% 
  filter(symbol %in% genes.OI) %>% 
  mutate(condition = recode(condition,
                            "MEDIA"="Media",
                            "TB"="Mtb"),
         xlab = paste(condition,Sample_Group, sep="\n")) %>% 
  mutate(symbol = factor(symbol, levels=genes.OI))

#### Contrast model data ####
#pull FDR from lm model if not available in lme
miss.fdr <- model_con$lme.contrast %>% 
  filter(FDR == "NaN") %>% 
  select(gene, variable, contrast_ref, contrast_lvl)

HIV_pos_con_lm <- model_con$lm.contrast %>% 
  inner_join(miss.fdr)

HIV_pos_contrast <- model_con$lme.contrast %>% 
  anti_join(miss.fdr) %>% 
  bind_rows(HIV_pos_con_lm) %>% 
  unnest(symbol) %>% 
  filter(symbol %in% genes.OI) %>% 
  select(variable, contrast_ref, contrast_lvl, gene, symbol, FDR) 

#### FDR annotation ####
dat.anno <- HIV_pos_contrast %>% 
  #Make symbols for plots
  mutate(sig_symbol = case_when(FDR <= 0.05 ~"**",
                            FDR <= 0.2~"*")) %>% 
  filter(!is.na(sig_symbol)) %>% 
  #x labels
  mutate(across(contrast_ref:contrast_lvl,
                ~recode(.,
                        "MEDIA LTBI"="Media\nLTBI",
                        "MEDIA RSTR"="Media\nRSTR",
                        "TB LTBI"="Mtb\nLTBI",
                        "TB RSTR"="Mtb\nRSTR"))) %>%
  rename(group1=contrast_ref, group2=contrast_lvl,
         ensembl_gene_id=gene) %>% 
  dplyr::select(ensembl_gene_id, symbol, group1, group2, sig_symbol)

#Add y location for pval based on max expression in plot
dat.anno.max <- dat %>% 
  #Max expression per gene and experiment
  group_by(ensembl_gene_id, symbol) %>% 
  summarise(max.e = max(value, na.rm=TRUE)) %>% 
  ungroup() %>% 
  right_join(dat.anno) %>% 
  arrange(symbol, group1, group2)

#first entry per gene
#Set y position to 1
first <- dat.anno.max %>% 
  group_by(symbol) %>% 
  slice(1) %>% 
  mutate(y.position1 = 1)

#Add first position data back and fill in remaining
dat.anno.max.ord <- dat.anno.max %>% 
  full_join(first) %>% 
  
  #Fill in positions 2 - N
  mutate(symbol = factor(symbol, levels = genes.OI)) %>% 
  group_by(symbol) %>% 
  mutate(y.position2 = lag(y.position1)+1,
         y.position3 = lag(y.position2)+1,
         y.position4 = lag(y.position3)+1) %>% 
  #Collapse positions into 1 column
  mutate(y.position = ifelse(!is.na(y.position1),y.position1,
                             ifelse(!is.na(y.position2),y.position2,
                                    ifelse(!is.na(y.position3),y.position3,
                                           ifelse(!is.na(y.position4),y.position4,NA))))) %>% 
  #Scale to max expression
  mutate(y.position = max.e+y.position) %>% 
  ungroup()

#### Gene boxplots ####
p3 <- dat %>% 
  
  ggplot(aes(x=xlab, y=value, color=Sample_Group)) +
  geom_jitter(width=0.1, height=0) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.25) +
  stat_summary(fun=mean, geom="errorbar", 
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               color="black", width=0.5) +
  facet_wrap(~symbol, scales="free", nrow=1) +
  #Add FDR
  stat_pvalue_manual(data=dat.anno.max.ord, 
                     label="sig_symbol", 
                     xmin="group1", xmax="group2") +
  #Beautify
  theme_bw() +
  labs(x="", y="Normalized log2 expression") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background =element_rect(fill="white")) +
  scale_color_manual(values = c("#E1BE6A", "#40B0A6"))

# p3 

#### Check colors ####
# library(colorblindr)
# cvd_grid(p2)
# cvd_grid(p3)

#### Save ####

p <- p1 + p2 + p3 + 
  plot_annotation(tag_levels = "A",
                  tag_prefix = "(", tag_suffix = ")") +
  plot_layout(widths = c(1.2,2,1.5))
# p

ggsave("publication/Fig1.DEG.pdf", p, width = 12, height = 3.5)
