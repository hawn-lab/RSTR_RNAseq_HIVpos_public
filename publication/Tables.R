library(tidyverse)
library(limma)
library(openxlsx)

#### Table 1. demographics ####
attach("data_clean/HIV_pos_voom.RData")
#Numeric summary
dat.voom$targets %>% 
  filter(condition == "MEDIA") %>% 
  select(Sample_Group, KCHCA_AGE_YR_CURRENT, RISK_SCORE, CD4) %>% 
  pivot_longer(-Sample_Group) %>% 
  group_by(Sample_Group, name) %>% 
  summarise(mean = mean(value),
            stdev = sd(value),
            med = median(value),
            min = min(value),
            max = max(value))
#Numeric stats
dat.voom$targets %>% 
  filter(condition == "MEDIA") %>% 
  select(Sample_Group, KCHCA_AGE_YR_CURRENT, RISK_SCORE, CD4) %>% 
  summarise(across(KCHCA_AGE_YR_CURRENT:CD4, 
                   funs(wilcox.test(.[Sample_Group == "LTBI"], 
                                    .[Sample_Group == "RSTR"])$p.value)))

#Categorical
dat.voom$targets %>% 
  filter(condition == "MEDIA") %>% 
  select(Sample_Group, M0_KCVSEX) %>% 
  count(Sample_Group, M0_KCVSEX) %>% 
  pivot_wider(names_from = M0_KCVSEX, values_from = n) %>% 
  mutate(percent_female = F/(F+M)*100)
#Categorical stats
dat.voom$targets %>% 
  filter(condition == "MEDIA") %>% 
  count(Sample_Group, M0_KCVSEX) %>% 
  pivot_wider(names_from = M0_KCVSEX, values_from = n) %>% 
  column_to_rownames("Sample_Group") %>% 
  chisq.test()

#### Table S1. model results ####
#model results
attach("results/HIVpos_lme.RData")

# Main model
fdr.main <- model_intr$lme %>% 
  rename(ensembl_gene_id=gene) %>% 
  select(symbol, ensembl_gene_id, 
         variable, estimate, pval, FDR) %>% 
  mutate(variable = gsub("Sample_Group","RSTR",variable),
         variable = gsub("condition","Mtb",variable))

#Contrast model
## Fill in failed lme fit with lm
miss.fdr <- model_con$lme.contrast %>% 
  filter(FDR == "NaN") %>% 
  select(gene, variable, contrast_ref, contrast_lvl)

cont_lm <- model_con$lm.contrast %>% 
  inner_join(miss.fdr)

fdr.cont <- model_con$lme.contrast %>% 
  anti_join(miss.fdr) %>% 
  bind_rows(cont_lm) %>% 
  rename(ensembl_gene_id=gene) %>% 
  select(symbol, ensembl_gene_id, 
         variable, contrast_ref, contrast_lvl, estimate, pval, FDR) %>% 
  mutate(variable = gsub("Sample_Group","RSTR",variable),
         variable = gsub("condition","Mtb",variable))

fdr.main.int <- fdr.main %>% 
  filter(variable == "Mtb:RSTR" & FDR < 0.2)

fdr.main.rstr <- fdr.main %>% 
  filter(variable == "RSTR" & FDR < 0.2)

fdr.main.mtb <- fdr.main %>% 
  filter(variable == "Mtb" & FDR < 0.2)

# Risk model
fdr.risk <- read_csv("results/model_fitting/HIVpos_risk_model.csv") %>% 
  left_join(dat.voom$genes, by = c("gene"="ensembl_gene_id")) %>% 
  select(condition, symbol, variable, estimate, pval, FDR)


sheets <- list("(A) RSTR_FDR<0.2"=fdr.main.rstr,
               "(B) interaction_FDR<0.2"=fdr.main.int,
               "(C) Mtb_FDR<0.2"=fdr.main.mtb,
               "(D) main_model"=fdr.main,
               "(E) contrast_model"=fdr.cont,
               "(F) risk_model"=fdr.risk)

write.xlsx(sheets, file = "publication/TableS1.linear_model_results.xlsx")

#### Table S2. DEG enrichment and GSEA ####
#ENRICHMENT
attach("results/enrich/RSTR_DEG_hypergeo.RData")

enrich_H %>% 
  filter(group == "All_DEG" & FDR < 0.3 & group_in_pathway > 1)
#None therefore, create place holder
H <- data.frame(group = "RSTR_and_interaction_DEG",
                gs_cat = "H",
                pathway = "No significant gene sets")

C2 <- enrich_C2 %>% 
  filter(group == "All_DEG" & FDR < 0.3 & group_in_pathway > 1) %>% 
  select(group:gs_subcat, pathway, group_in_pathway, `k/K`, pval, FDR, genes) %>% 
  mutate(genes = as.character(genes)) %>% 
  mutate(group = "RSTR_and_interaction_DEG")

#None for Uniprot as well
U <- data.frame(group = "RSTR_and_interaction_DEG",
                gs_cat = "Uniprot_keyword",
                pathway = "No significant gene sets")

#GSEA
attach("results/enrich/HIVpos_GSEA.RData")
H_gsea <- GSEA_H %>% 
  filter(FDR < 0.1) %>% 
  select(group, gs_cat, pathway, NES, pval, FDR, leadingEdge) %>% 
  mutate(leadingEdge = as.character(leadingEdge)) 

#Save
sheets_deg <- list("(A) RSTR_DEG_Hallmark" = H,
                   "(B) RSTR_DEG_C2" = C2,
                   "(C) RSTR_DEG_UniProt" = U,
                   "(D) RSTR_GSEA_Hallmark" = H_gsea)

write.xlsx(sheets_deg, file = "publication/TableS2.RSTR.DEG.enrich.gsea.xlsx")

#### Table S3. Mtb-HIV specific DEG and enrichment ####
#gene anno
attach("data_clean/HIV_pos_voom.RData")
gene.key <- dat.voom$genes %>% 
  select(ensembl_gene_id, symbol)

#specific DEG
mtb_deg <- read_csv("results/HIVspecific_DEG.csv") %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "group", values_to = "ensembl_gene_id") %>% 
  left_join(gene.key) %>% 
  mutate(group = recode(group,
                        "HIVpos_only_Mtb_DEG"="PLWH_specific_Mtb",
                        "HIVneg_only_Mtb_DEG"="woHIV_specific_Mtb")) %>% 
  select(-rowname) %>% 
  drop_na(ensembl_gene_id) %>% 
  arrange(group, ensembl_gene_id)

#ENRICHMENT
attach("results/enrich/HIVpos_specific_hypergeo.RData")

enrich_H %>% 
  filter(FDR < 0.3 & group_in_pathway > 1)
#None therefore, create place holder
H <- data.frame(group = c("PLWH_specific_Mtb","woHIV_specific_Mtb"),
                gs_cat = "H",
                pathway = "No significant gene sets")

C2 <- enrich_C2 %>% 
  filter(FDR < 0.3 & group_in_pathway > 1) %>% 
  select(group:gs_subcat, pathway, group_in_pathway, `k/K`, pval, FDR, genes) %>% 
  mutate(genes = as.character(genes)) %>% 
  mutate(group = recode(group,
                        "HIVpos_specific_Mtb"="PLWH_specific_Mtb",
                        "HIVneg_specific_Mtb"="woHIV_specific_Mtb"))

U <- enrich_unip %>% 
  filter(FDR < 0.3 & group_in_pathway > 1) %>% 
  select(group:gs_subcat, pathway, group_in_pathway, `k/K`, pval, FDR, genes) %>%
  mutate(genes = as.character(genes)) %>% 
  mutate(group = recode(group,
                        "HIVpos_specific_Mtb"="PLWH_specific_Mtb",
                        "HIVneg_specific_Mtb"="woHIV_specific_Mtb"))

#Save
sheets_deg2 <- list("(A) Mtb_DEG" = mtb_deg,
                    "(B) Mtb_DEG_Hallmark" = H,
                   "(C) Mtb_DEG_C2" = C2,
                   "(D) Mtb_DEG_UniProt" = U)

write.xlsx(sheets_deg2, file = "publication/TableS3.Mtb.DEG.enrich.xlsx")
