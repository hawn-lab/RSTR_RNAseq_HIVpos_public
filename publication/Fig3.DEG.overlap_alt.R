library(tidyverse)
library(patchwork)

#### DATA ####
HIV_neg_results <- read_csv("results/HIVneg/RSTR.Mtb.model.results.anno.corrected.csv") %>% 
  mutate(model = "non-PLWH") %>% 
  select(model, variable, gene, symbol, FDR) %>% 
  mutate(variable = recode(variable, 
                           "conditionTB"="condition",
                           "conditionTB:Sample_GroupRSTR" = "condition:Sample_Group",
                           "Sample_GroupRSTR"="Sample_Group")) 

attach("results/HIVpos_lme.RData")
HIV_pos_results <- model_intr$lme %>% 
  mutate(model = "PLWH")

#### Format data ####

##Condition
dat1 <- full_join(
  select(HIV_neg_results, -symbol), 
  select(HIV_pos_results, -symbol),
  by=c("variable", "gene"))  %>% 
  rename(HIV=FDR.x, PLWH=FDR.y) %>% 
  filter(variable == "condition") %>% 
  mutate(across(c(HIV, PLWH), 
                ~case_when(. < 0.01 ~ "s",
                           . > 0.2 ~ "ns",
                           TRUE~"other"))) %>% 
  distinct() %>% 
  count(HIV, PLWH) %>% 
  #Remove NS in both
  filter(HIV=="s" | PLWH=="s") %>% 
  mutate(group="mtb")

##interaction
dat2 <- full_join(
  select(HIV_neg_results, -symbol), 
  select(HIV_pos_results, -symbol),
  by=c("variable", "gene")) %>% 
  rename(HIV=FDR.x, PLWH=FDR.y) %>% 
  filter(variable == "condition:Sample_Group") %>% 
  mutate(across(c(HIV, PLWH), 
                ~case_when(. < 0.2 ~ "s",
                           . > 0.4 ~ "ns",
                           TRUE~"other"))) %>% 
  distinct() %>% 
  count(HIV, PLWH) %>% 
  #Remove NS in both
  filter(HIV=="s" | PLWH=="s") %>% 
  mutate(group="interact")

##combine
dat_all <- bind_rows(dat1,dat2) %>% 
  #format labels for plot
  mutate(name = paste(HIV,PLWH,sep="_"),
         name = recode_factor(name,
                              "ns_s"="PLWH specific DEGs",
                              "other_s"="PLWH non-specific DEGs",
                              "s_s"="Shared DEGs",
                              "s_other"="non-PLWH non-specific DEGs",
                              "s_ns"="non-PLWH specific DEGs")) %>% 
  mutate(group = recode_factor(group,
                               "mtb"="Mtb-infected vs media",
                               "interact"="Mtb:RSTR interaction"))

#### Plot ####

p1 <- dat_all %>% 
  ggplot() +
  aes(x=group, y=n, fill=name, label=n) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  coord_flip() +
  facet_wrap(~group, scales = "free", ncol=1) +
  theme_classic()  +
  theme(strip.text.x = element_blank(),
        panel.border = element_rect(fill=NA))+
  labs(y="Total genes", x="", fill="") +
  scale_fill_manual(values = c('#DC267F','#FE6100','#FFB000','#648FFF','#785EF0'))
p1

#### Save ####
ggsave(p1, filename="publication/Fig3.DEG.overlap_alt.pdf",
       height=2, width=8.5)
ggsave(p1, filename="publication/Fig3.DEG.overlap_alt.png",
       height=2, width=8.5)
