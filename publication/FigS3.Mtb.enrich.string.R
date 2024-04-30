library(tidyverse)
library(BIGpicture)
library(patchwork)

#### Enrichment data ####
attach("results/enrich/HIVpos_specific_hypergeo.RData")

#C2 select top signif only. This is all HIV- and a subset of HIV+
#Ends up being FDR < 0.05 for HIV+
enrich_C2_top <- enrich_C2 %>% 
  select(group, gs_cat, gs_subcat, pathway, `k/K`, 
         FDR, group_in_pathway, genes) %>%
  #Significant FDR < 0.05 and overlap > 1
  filter(FDR < 0.05 & group_in_pathway > 1) %>% 
  #top hits to reduce HIV+ list
  group_by(group, gs_cat, gs_subcat) %>% 
  slice_min(order_by = FDR, n = 13, with_ties = TRUE)

enrich_signif <- bind_rows(enrich_C2_top, enrich_unip) %>% 
  #Significant FDR < 0.05 and overlap > 1
  filter(FDR < 0.05 & group_in_pathway > 1) %>% 
  
  #Format facet labels
  mutate(gs_cat = recode(gs_cat, "C2"="MSigDB C2", "UniProt_kw"="UniProt")) %>% 
  mutate(gs_subcat = gsub("^CP:","",gs_subcat),
         gs_subcat = ifelse(gs_subcat=="","keyword",gs_subcat)) %>% 
  #Clean pw names
  mutate(pathway = gsub("^REACTOME_|^KEGG_","", pathway),
         pathway = tolower(pathway),
         pathway = gsub("_"," ", pathway),
         pathway = gsub("dna-","DNA-", pathway),
         pathway = gsub("dna ","DNA ", pathway),
         pathway = gsub("rho gtpase","Rho GTPase", pathway),
         pathway = gsub(" hox "," HOX ", pathway),
         pathway = gsub("beta catenin tcf","Beta-catenin TCF", pathway),
         pathway = gsub("sirt1","SIRT1", pathway),
         pathway = gsub(" rrna "," rRNA ", pathway),
         pathway = gsub("pkn1","PKN1", pathway),
         pathway = gsub("klk2","KLK2", pathway),
         pathway = gsub("klk3","KLK3", pathway),
         pathway = gsub("ap site","AP site", pathway),
         pathway = gsub("runx1","RUNX1", pathway),
         pathway = gsub(" orc "," ORC ", pathway),
         pathway = gsub("ps1","PS1", pathway),
         pathway = gsub("sm path","SM path", pathway),
         pathway = gsub("hes path","HES path", pathway),
         pathway = gsub("hats ace","HATs ace", pathway),
         pathway = gsub("ar androgen receptor","androgen\nreceptor", pathway)) %>%
  mutate(pathway = gsub(" with ","\nwith ",pathway),
         pathway = gsub(" differentiation","\ndifferentiation",pathway),
         pathway = gsub(" development","\ndevelopment",pathway)) %>% 
  #color groups
  mutate(Significance = case_when(
    FDR < 0.05 ~ "FDR < 0.05"),
    Significance = factor(Significance))

#remove terms in both HIV groups
HIVpos_signif <- enrich_signif %>% 
  filter(group == "HIVpos_specific_Mtb") %>% 
  pull(pathway)

HIVneg_signif <- enrich_signif %>% 
  filter(group == "HIVneg_specific_Mtb") %>% 
  pull(pathway)

enrich_signif_spec <- enrich_signif %>% 
  filter(!pathway %in% intersect(HIVpos_signif, HIVneg_signif))

#### Enrichment plot: HIV+ ####
##C2
p1b <- enrich_signif_spec %>% 
  filter(group == "HIVpos_specific_Mtb") %>% 
  filter(gs_cat == "MSigDB C2") %>% 
  
  ggplot(aes(stats::reorder(pathway, `k/K`), `k/K`)) +
  geom_segment(aes(stats::reorder(pathway, `k/K`), 
                   xend = pathway, 
                   y = 0, yend = `k/K`)) + 
  geom_point(size = 3, aes(fill = Significance), 
             shape = 21, stroke = 1) + 
  geom_hline(yintercept = 0) + 
  # facet_grid(gs_subcat~., scales = "free", space = "free") +
  scale_fill_manual(values=c("FDR < 0.3" = "#D55E00",
                             "FDR < 0.05" = "#FBD769",
                             "NS" = "grey90"), drop=FALSE)+
  coord_flip() + 
  labs(x="", y = "Proportion enriched ( k/K )",
       subtitle = "MSigDB C2 REACTOME") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(lineheight = 0.8)) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0,0.24))

##UniProt
p1a <- enrich_signif_spec %>% 
  filter(group == "HIVpos_specific_Mtb") %>% 
  filter(gs_cat == "UniProt") %>% 
  
  ggplot(aes(stats::reorder(pathway, `k/K`), `k/K`)) +
  geom_segment(aes(stats::reorder(pathway, `k/K`), 
                   xend = pathway, 
                   y = 0, yend = `k/K`)) + 
  geom_point(size = 3, aes(fill = Significance), 
             shape = 21, stroke = 1) + 
  geom_hline(yintercept = 0) + 
  # facet_grid(gs_subcat~., scales = "free", space = "free") +
  scale_fill_manual(values=c("FDR < 0.3" = "#D55E00",
                             "FDR < 0.05" = "#FBD769",
                             "NS" = "grey90"))+
  coord_flip() + 
  labs(x="", y = "",
       title = "(A) PLWH, Mtb-HIV specific DEG",
       subtitle = "UniProt keyword") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(lineheight = 0.8),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0,0.24))

#### Enrichment plot: HIV- ####
##C2
p2b <- enrich_signif_spec %>% 
  filter(group == "HIVneg_specific_Mtb") %>% 
  filter(gs_cat == "MSigDB C2") %>% 
  
  ggplot(aes(stats::reorder(pathway, `k/K`), `k/K`)) +
  geom_segment(aes(stats::reorder(pathway, `k/K`), 
                   xend = pathway, 
                   y = 0, yend = `k/K`)) + 
  geom_point(size = 3, aes(fill = Significance), 
             shape = 21, stroke = 1) + 
  geom_hline(yintercept = 0) + 
  # facet_grid(gs_subcat~., scales = "free", space = "free") +
  scale_fill_manual(values=c("FDR < 0.3" = "#D55E00",
                             "FDR < 0.05" = "#FBD769",
                             "NS" = "grey90"), drop=FALSE)+
  coord_flip() + 
  labs(x="", y = "Proportion enriched ( k/K )",
       subtitle = "MSigDB C2 REACTOME") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(lineheight = 0.8)) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0,0.24))

##UniProt
p2a <- enrich_signif_spec %>% 
  filter(group == "HIVneg_specific_Mtb") %>% 
  filter(gs_cat == "UniProt") %>% 
  
  ggplot(aes(stats::reorder(pathway, `k/K`), `k/K`)) +
  geom_segment(aes(stats::reorder(pathway, `k/K`), 
                   xend = pathway, 
                   y = 0, yend = `k/K`)) + 
  geom_point(size = 3, aes(fill = Significance), 
             shape = 21, stroke = 1) + 
  geom_hline(yintercept = 0) + 
  # facet_grid(gs_subcat~., scales = "free", space = "free") +
  scale_fill_manual(values=c("FDR < 0.3" = "#D55E00",
                             "FDR < 0.05" = "#FBD769",
                             "NS" = "grey90"))+
  coord_flip() + 
  labs(x="", y = "",
       title = "(C) non-PLWH, Mtb-HIV specific DEG",
       subtitle = "UniProt keyword") + 
  theme_bw() +
  theme(strip.background =element_rect(fill="white"),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(lineheight = 0.8),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0,0.24))


#### STRING data ####
#Format enrichment to symbol
attach("data_clean/HIV_pos_voom.RData")
gene_key <- dat.voom$genes %>% 
  select(ensembl_gene_id, symbol) %>% 
  unnest(symbol)

enrich_signif_symbol <- enrich_signif_spec %>% 
  ungroup() %>% 
  unnest(genes) %>% 
  left_join(gene_key, by = c("genes"="ensembl_gene_id")) %>% 
  mutate(symbol = ifelse(symbol == "NULL" | symbol == "NA" |
                           is.na(symbol), genes, symbol)) %>%
  select(-genes) %>%
  group_by(across(-symbol)) %>%
  summarise(genes = paste(unique(symbol),
                          sep = ",", collapse = ",")) %>% 
  mutate(genes = strsplit(genes, split = ",")) %>% 
  ungroup()

#collapse like terms
txn_genes <- filter(enrich_signif_symbol, 
                    pathway %in% c("transcription",
                                   "transcription regulation")) %>% 
  unnest(genes) %>% 
  pull(genes) %>% unique()
zn_genes <- filter(enrich_signif_symbol, 
                   pathway %in% c("zinc", "zinc-finger")) %>% 
  unnest(genes) %>% 
  pull(genes) %>% unique()

enrich_signif_symbol_coll <- enrich_signif_symbol %>% 
  add_row(group="HIVpos_specific_Mtb", gs_cat="UniProt", gs_subcat="keyword",
          pathway="transcription / transcription regulation", 
          group_in_pathway = Inf, FDR = 0, 
          genes = list(txn_genes)) %>% 
  add_row(group="HIVpos_specific_Mtb", gs_cat="UniProt", gs_subcat="keyword",
          pathway="zinc / zinc-finger", 
          group_in_pathway = Inf, FDR = 0, 
          genes = list(zn_genes)) %>% 
  filter(!pathway %in% c("transcription",
                         "transcription regulation",
                         "zinc", "zinc-finger"))


#list DEG
deg <- read_csv("results/HIVspecific_DEG.csv")

#convert to symbol
deg_symbol <- deg %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "group", values_to = "ensembl_gene_id") %>% 
  drop_na(ensembl_gene_id) %>% 
  left_join(gene_key)

deg_sym_HIVpos <- deg_symbol %>% 
  filter(group == "HIVpos_only_Mtb_DEG") %>% 
  pull(symbol) %>% unique()

#Map to string
map.HIVpos.400 <- map_string(genes = deg_sym_HIVpos, score_threshold = 400)
map.HIVpos.700 <- map_string(genes = deg_sym_HIVpos, score_threshold = 700)

#Get main cluster
library(bc3net)
map.HIVpos.400.lrg <- map.HIVpos.400
map.HIVpos.400.lrg$subgraph <- getgcc(map.HIVpos.400.lrg$subgraph)

#Get small clusters
map.HIVpos.400.sm <- map.HIVpos.400
map.HIVpos.400.sm$subgraph <- delete_vertices(map.HIVpos.400.sm$subgraph,
                               vertex_attr(map.HIVpos.400.lrg$subgraph)$name)

#### STRING plot ####
col.vec <- c('#E69F00','#76CBFB','#009E73','#F0E442','#4498C7','grey90')
names(col.vec) <- c("DNA-binding","mitochondrion nucleoid","nucleus",
                    "transcription / transcription regulation",
                    "zinc / zinc-finger", "none")
                      
p3.orphan <- plot_string(map.HIVpos.400, 
                         enrichment = filter(enrich_signif_symbol_coll,
                                             group == "HIVpos_specific_Mtb" &
                                               gs_cat == "UniProt"), 
                         fdr.cutoff = 0.05, discard = "cluster",
                         enriched.only = TRUE, layout = "graphopt",
                         colors = col.vec[c(1,3,4,5)], 
                         text_size = 2, node_size = 0.5) +
  theme(legend.position = "none")

p3.lrg <- plot_string(map.HIVpos.400.lrg, 
                      enrichment = filter(enrich_signif_symbol_coll,
                                          group == "HIVpos_specific_Mtb" &
                                            gs_cat == "UniProt"), 
                      fdr.cutoff = 0.05, discard = "orphan",
                      layout = "lgl",
                      colors = col.vec, text_size = 2, node_size = 0.5) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p3.sm <- plot_string(map.HIVpos.400.sm, 
                     enrichment = filter(enrich_signif_symbol_coll,
                                         group == "HIVpos_specific_Mtb" &
                                           gs_cat == "UniProt"), 
                     fdr.cutoff = 0.05, discard = "orphan",
                     layout = "fr",
                     colors = col.vec, text_size = 2, node_size = 0.5) +
  theme(legend.position = "none")

#### Save ####
p.enrich <- p1a/p1b/p2a/p2b+
  plot_layout(heights = c(3,6,1,1))
p.enrich

ggsave(filename = "publication/FigS3.Mtb.enrich.pdf",
       p.enrich, width=7, height=10)

p.string <- p3.lrg+p3.sm+p3.orphan
p.string

ggsave(filename = "publication/FigS3.Mtb.string.pdf",
       p.string, width=22, height=10)
