library(tidyverse)
library(BIGpicture)
library(patchwork)

#colors
col.vec <-c("#56B4E9","#009E73","#0072B2", "#D55E00","#CC79A7", "grey80")
names(col.vec) <- c("extracellular matrix organization",
                    "IFN signaling",
                    "systemic lupus erythematosus",
                    "cell adhesion molecules (CAMs)",
                    "metabolism of vitamins and cofactors",
                    "none")

#### Enrichment data ####
attach("results/enrich/RSTR_DEG_hypergeo.RData")

enrich_signif <- enrich_C2 %>% 
  #DEG results
  filter(group == "All_DEG") %>% 
  #Significant FDR < 0.3 and overlap > 1
  filter(FDR < 0.3 & group_in_pathway > 1) %>% 
  select(group, gs_cat, gs_subcat, pathway, `k/K`, 
         FDR, group_in_pathway, genes) %>%
  #Format facet labels
  mutate(gs_cat = recode(gs_cat, "C7"="MSigDB C7")) %>% 
  mutate(gs_subcat = gsub("^CP:","",gs_subcat)) %>% 
  #Clean pw names
  mutate(pathway = gsub("^REACTOME_|^KEGG_","", pathway),
         pathway = tolower(pathway),
         pathway = gsub("_"," ", pathway),
         pathway = gsub("interferon","IFN", pathway),
         pathway = gsub(" cams$"," (CAMs)", pathway)) %>% 
  #color groups
  mutate(Significance = case_when(
    FDR < 0.05 ~ "FDR < 0.05",
    FDR < 0.3 ~ "FDR < 0.3"),
    Significance = factor(Significance))

#### Enrichment plot ####
p1 <- enrich_signif %>% 
  ggplot(aes(stats::reorder(pathway, `k/K`), `k/K`)) +
  geom_segment(aes(stats::reorder(pathway, `k/K`), 
                   xend = pathway, 
                   y = 0, yend = `k/K`)) + 
  geom_point(size = 3, aes(fill = Significance), 
             shape = 21, stroke = 1) + 
  geom_hline(yintercept = 0) + 
  facet_grid(gs_subcat~., scales = "free", space = "free") +
  scale_fill_manual(values=c("FDR < 0.3" = "#D55E00",
                             "NS" = "grey90"))+
  coord_flip() + 
  labs(x="", y = "Proportion enriched ( k/K )") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(lineheight = 0.8)) +
  scale_y_continuous(breaks = c(0, 0.005, 0.010))
# p1

#### STRING data ####
#Format enrichment to symbol
attach("data_clean/HIV_pos_voom.RData")
gene_key <- dat.voom$genes %>% 
  select(ensembl_gene_id, symbol)

enrich_signif_symbol <- enrich_signif %>% 
  unnest(genes) %>% 
  left_join(gene_key, by = c("genes"="ensembl_gene_id")) %>% 
  select(-genes) %>% 
  group_by(across(-symbol)) %>%
  summarise(genes = paste(unique(symbol), 
                          collapse = ",", sep = ",")) %>% 
  mutate(genes = strsplit(genes, split = ",")) %>% 
  ungroup()

#list DEG
attach("results/HIVpos_lme.RData")

deg_rstr <- model_intr$lme %>% 
  filter(grepl("Sample_Group", variable) & FDR < 0.2) %>%
  unnest(symbol) %>% 
  pull(symbol)

#Map to string
map.rstr.400 <- map_string(genes = deg_rstr, score_threshold = 400)

#### STRING plot ####
p2 <- plot_string(map.rstr.400, enrichment = enrich_signif_symbol, 
                  layout = "grid",
                  fdr.cutoff = 0.3,
                  colors = col.vec, text_size = 3, node_size = 2)
# p2

#### Save ####
p.all <- p2+p1 + 
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(widths = c(2,1))

ggsave(filename = "publication/FigS1.RSTR.enrich.string.pdf",
       p.all, width=13, height=2.5)
