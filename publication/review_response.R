library(tidyverse)
library(kimma)
library(BIGpicture)

load("data_clean/HIV_pos_voom.RData")

#### TB risk scores ####
# RSTR vs LTBI risk scores
dat.voom$targets %>% 
  summarise(risk =t.test(RISK_SCORE[Sample_Group == "LTBI"], 
                         RISK_SCORE[Sample_Group == "RSTR"])$p.value)
dat.voom$targets %>% 
  group_by(Sample_Group) %>% 
  summarise(m=mean(RISK_SCORE),
            s=sd(RISK_SCORE))

dat.voom$targets %>%
  ggplot(aes(x=Sample_Group, y=RISK_SCORE)) +
  geom_violin() +
  geom_jitter(width=0.2, height = 0) +
  stat_summary(geom = "point", fun.y = "mean", col = "red",
               size = 3, shape = "square") +
  theme_classic()

# Risk score and gene expression

#overlap with previous DEGs?
attach("results/model_fitting/HIVpos_models.RData")
deg_intr <- model_intr$lme %>% 
  filter(FDR < 0.2 & variable %in% c("Sample_Group","condition:Sample_Group")) %>% 
  pull(gene) %>% 
  unique()



desk_risk <- risk_all %>% 
  filter(FDR<0.2) %>% 
  left_join(dat.voom$genes,by=c("gene"="ensembl_gene_id")) %>% 
  unnest(symbol) %>% 
  select(condition, variable, gene, symbol, estimate, pval, FDR)

deg_risk

# Plot overlap degs
p1 <- as.data.frame(dat.voom$E) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  pivot_longer(-ensembl_gene_id, names_to = "libID") %>% 
  left_join(dat.voom$genes) %>% 
  unnest(symbol) %>% 
  filter(symbol %in% deg_risk$symbol) %>% 
  left_join(dat.voom$targets) %>% 
  left_join(risk_all, by=c("condition"="condition","ensembl_gene_id"="gene")) %>% 
  mutate(facet_lab = paste(condition,"FDR =", signif(FDR, digits=2))) %>% 
  
  ggplot(aes(x=RISK_SCORE,y=value)) +
  geom_point(aes(color=Sample_Group)) +
  geom_smooth(method="lm", se=FALSE, formula = 'y ~ x', color="grey") + 
  # geom_smooth(aes(color=Sample_Group), method="lm", se=FALSE, formula = 'y ~ x') +
  theme_classic() +
  facet_wrap(symbol~facet_lab, scales="free") +
  labs(y="Log2 normalized expression", x="TB risk score", color="")
# p1

ggsave(p1, filename="publication/risk_score_deg.png", width=5.5, height=5)

# correlation
# library(ggpubr)
# 
# as.data.frame(dat.voom$E) %>% 
#   rownames_to_column("ensembl_gene_id") %>% 
#   pivot_longer(-ensembl_gene_id, names_to = "libID") %>% 
#   left_join(dat.voom$genes) %>% 
#   unnest(symbol) %>% 
#   filter(ensembl_gene_id %in% deg_intr$gene) %>% 
#   left_join(dat.voom$targets) %>% 
#   
#   ggplot(aes(x=RISK_SCORE,y=value,color=condition)) +
#   geom_point() +
#   geom_smooth(method="lm", se=FALSE, formula = 'y ~ x') + 
#   theme_classic() +
#   facet_wrap(~symbol, scales="free") +
#   labs(y="Log2 normalized expression", x="TB risk score", color="") +
#   stat_cor()
