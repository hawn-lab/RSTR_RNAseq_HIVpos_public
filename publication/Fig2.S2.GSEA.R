library(tidyverse)
library(patchwork)

#### Data ####
fdr.cutoff <- 0.1

#HIV+ results
attach("results/enrich/HIVpos_GSEA.RData")

GSEA_HIVpos <- GSEA_H %>% 
  #Make facet and offset groups
  mutate(group = case_when(group == "MEDIA RSTR-MEDIA LTBI" ~ 
                             "RSTRinMEDIA",
                           group == "TB RSTR-TB LTBI" ~ 
                             "RSTRinTB",
                           group == "TB LTBI-MEDIA LTBI" ~ 
                             "TBinLTBI",
                           group == "TB RSTR-MEDIA RSTR" ~ 
                             "TBinRSTR")) %>% 
  mutate(cohort = "Uganda", HIV = "PLWH") %>% 
  select(-leadingEdge)

#HIV- results
GSEA_HIVneg <- read_csv("results/HIVneg/h_GSEA_U_rstr_neg.result.csv") %>% 
  mutate(cohort = "Uganda") %>% 
  # bind_rows(
  #   read_csv("results/HIVneg/h_GSEA_SA_rstr_neg.result.csv") %>% 
  #     mutate(cohort = "South Africa")
  # ) %>% 
  mutate(HIV = "non-PLWH") %>% 
  rename_all(~gsub("fgsea.","",.)) %>% 
  select(-leadingEdge)

#### PW of interest ####
#PW signif in HIV+
rstr_signif <- GSEA_HIVpos %>% 
  filter(group %in% c("RSTRinMEDIA", "RSTRinTB")) %>% 
  filter(FDR < fdr.cutoff) %>% 
  pull(pathway) %>% unique()

#PW signif in U and SA HIV-
attach("results/HIVneg/GSEA.signif.RData")

neg_signif <- h.signif %>% 
  pull(pathway) %>% unique()

#All unique signif PW
all_signif <- unique(c(rstr_signif, neg_signif))
## Remove non-consistent direction from HIV- results (section iv) that are also not signif in HIV+
all_signif2 <- all_signif[!all_signif %in%
                            c("HALLMARK_KRAS_SIGNALING_DN",
                              "HALLMARK_HYPOXIA",
                              "HALLMARK_IL2_STAT5_SIGNALING")]

#New terms for HIV+ data
rstr_signif[!rstr_signif %in% neg_signif]

#Combine, filter to pw of interest and beautify names
GSEA_clean <- bind_rows(GSEA_HIVpos, GSEA_HIVneg) %>% 
  filter(pathway %in% all_signif2) %>% 
  ##lowercase and then correction
  mutate(pathway = gsub("HALLMARK_","",pathway),
         pathway = gsub("_"," ", tolower(pathway)),
         pathway = gsub("tnfa","TNFA",pathway),
         pathway = gsub("nfkb","NF-kB",pathway),
         pathway = gsub("interferon","IFN",pathway),
         pathway = gsub("il2","IL2",pathway),
         pathway = gsub("stat5","STAT5",pathway),
         pathway = gsub("kras","KRAS",pathway),
         pathway = gsub(" dn"," down",pathway),
         pathway = gsub("tgf beta","TGF-beta",pathway),
         pathway = gsub("notch","NOTCH",pathway)) %>%
  ## order
  mutate(pathway = factor(pathway,
                          levels=c(
                            #Section i
                            "apical junction",
                            "inflammatory response",
                            "TNFA signaling via NF-kB",
                            "TGF-beta signaling",
                            #Section ii
                            "oxidative phosphorylation",
                            "adipogenesis",
                            #Section iii
                            "IFN gamma response",
                            "IFN alpha response",
                            "allograft rejection",
                            #Section iv
                            "IL2 STAT5 signaling",
                            "myogenesis",
                            "hypoxia",
                            "KRAS signaling down",
                            "epithelial mesenchymal transition",
                            "coagulation",
                            "NOTCH signaling"
                          ))) %>%
  mutate(facet_y = case_when(
    #Differences in HIV+ vs -
    pathway %in% c("TNFA signaling via NF-kB",
                   "inflammatory response",
                   "apical junction",
                   "TGF-beta signaling") ~ "i",
    #same in HIV+ and -
    pathway %in% c("adipogenesis",
                   "oxidative phosphorylation") ~ "ii",
    pathway %in% c("IFN alpha response",
                   "allograft rejection",
                   "IFN gamma response") ~ "iii",
    #Inconsistent HIV-
    TRUE ~ "iv")) %>%
  #Make significance group (color)
  mutate(Significance = case_when(FDR < fdr.cutoff ~ 
                                    paste("FDR <", fdr.cutoff),
                                  TRUE ~ "NS")) %>% 
  #order cohorts
  mutate(cohort = factor(cohort,
                         levels = c("Uganda",
                                    "South Africa")))

#format facet names
facet_lab <- c('RSTRinMEDIA'="Down in RSTR <-  -> Up in RSTR",
               'RSTRinTB'="Down in RSTR <-  -> Up in RSTR",
               'TBinLTBI'="Down in +Mtb <-  -> Up in +Mtb",
               'TBinRSTR'="Down in +Mtb <-  -> Up in +Mtb",
               'i'='i','ii'='ii','iii'='iii','iv'='iv','v'='v')

#### RSTR plot ####
p1_dat <- GSEA_clean %>% 
  filter(group %in% c("RSTRinMEDIA","RSTRinTB")) %>% 
  droplevels()

#plot lims
plot.lim <- max(abs(p1_dat$NES)) + 0.15

## HIV+ 
p1 <- p1_dat %>% 
  filter(HIV == "PLWH") %>% 
  
  ggplot(aes(x = pathway, y = NES)) +
  geom_linerange(aes(ymin = 0, ymax = NES, group = cohort), size = 0.5,
                 position = position_dodge(width = 0.4)) +
  geom_point(aes(fill = Significance),
             stroke = 0.7, shape = 21,
             position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=c("FDR < 0.1" = "#D55E00",
                             "NS" = "grey90")) +
  lims(y=c(-plot.lim,plot.lim)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized enrichment score (NES)",
       fill = "Significance", shape = "Cohort",
       title = "(A) PLWH Media               (B) PLWH +Mtb") + 
  facet_grid(facet_y ~ group, 
             labeller = as_labeller(facet_lab),
             scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

##HIV-
p2 <- p1_dat %>% 
  filter(HIV == "non-PLWH") %>% 
  
  ggplot(aes(x = pathway, y = NES)) +
  geom_linerange(aes(ymin = 0, ymax = NES, group = cohort), size = 0.5,
                 position = position_dodge(width = 0.4)) +
  geom_point(aes(fill = Significance),
             stroke = 0.7, shape = 21,
             position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=c("FDR < 0.1" = "#D55E00",
                             "NS" = "grey90")) +
  lims(y=c(-plot.lim,plot.lim)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized enrichment score (NES)",
       fill = "Significance", 
       title = "(C) non-PLWH Media        (D) non-PLWH +Mtb") + 
  facet_grid(facet_y ~ group, 
             labeller = as_labeller(facet_lab),
             scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank(),
        legend.position = "none")

# p1/p2

#### Mtb plot ####
p2_dat <- GSEA_clean %>% 
  filter(group %in% c("TBinRSTR","TBinLTBI")) %>% 
  mutate(Sample_Group = case_when(group=="TBinRSTR"~"RSTR",
                                  group=="TBinLTBI"~"LTBI"),
         group2 = paste(HIV, Sample_Group),
         group3 = "Down in +Mtb <-  -> Up in +Mtb  "
         ) %>% 
  droplevels()

#plot lims
plot.lim2 <- max(abs(p2_dat$NES)) + 0.15

## HIV+ 
p3 <- p2_dat %>% 

  ggplot(aes(x = pathway, y = NES, group = group2)) +
  geom_linerange(aes(ymin = 0, ymax = NES), size = 0.5,
                 position = position_dodge(width = 0.8)) +
  geom_point(aes(fill = Significance, shape = group2),
             stroke = 0.7, 
             position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=c("FDR < 0.1" = "#D55E00",
                             "NS" = "grey90")) +
  scale_shape_manual(values = c(23,24,22,25)) +
  lims(y=c(-plot.lim2,plot.lim2)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized enrichment score (NES)",
       fill = "Significance", shape = "Cohort") + 
  facet_grid(facet_y~group3, 
             scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.direction = "vertical") +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(ncol=2, reverse = TRUE))
# p3

##HIV-
# p4 <- p2_dat %>% 
#   filter(HIV == "non-PLWH") %>% 
#   
#   ggplot(aes(x = pathway, y = NES)) +
#   geom_linerange(aes(ymin = 0, ymax = NES, group = cohort), size = 0.5,
#                  position = position_dodge(width = 0.4)) +
#   geom_point(aes(fill = Significance),
#              stroke = 0.7, shape = 21,
#              position = position_dodge(width = 0.4)) +
#   geom_hline(yintercept = 0) +
#   scale_fill_manual(values=c("FDR < 0.1" = "#D55E00",
#                              "NS" = "grey90")) +
#   lims(y=c(-plot.lim2,plot.lim2)) +
#   coord_flip() +
#   labs(x = "Pathway", y = "Normalized enrichment score (NES)",
#        fill = "Significance", shape = "Cohort",
#        title = "(C) non-PLWH LTBI           (D) non-PLWH RSTR") + 
#   facet_grid(facet_y ~ group, 
#              labeller = as_labeller(facet_lab),
#              scales = "free_y", space = "free_y") +
#   theme_bw() +
#   theme(strip.background = element_rect(fill="white"),
#         panel.grid.minor = element_blank(),
#         legend.position = "none")

# p3/p4

#### Check colors ####
# library(colorblindr)
# cvd_grid(p1)
# cvd_grid(p2)

#### Save ####
ggsave(filename="publication/Fig2.GSEA.RSTR.pdf",
       p1/p2, width = 6.8, height=7)
ggsave(filename="publication/Fig2.GSEA.RSTR.png",
       p1/p2, width = 6.8, height=7)

# ggsave(filename="publication/FigS2.GSEA.MTB.pdf",
#        p3/p4, width = 6.8, height=7)
# ggsave(filename="publication/FigS2.GSEA.MTB.png",
#        p3/p4, width = 6.8, height=7)

ggsave(filename="publication/FigS2.GSEA.MTB.pdf",
       p3, width = 5, height=5.5)
# ggsave(filename="publication/FigS2.GSEA.MTB.png",
#        p3, width = 5, height=5.5)

#### Save pw for STRING ####
gsea_pw <- all_signif2

save(gsea_pw, file = "publication/GSEA_pw.RData")
