# setup -------------------------------------------------------------------
library(tidyverse)
library(here)
library(forcats)
library(broom)
library(broom.mixed)
library(glmmTMB)
library(kableExtra)

options(scipen=999)
theme_set(theme_bw(base_size = 14) + theme(legend.position = "top"))
knitr::opts_chunk$set(echo = FALSE, message = FALSE, fig.width = 9,
                      fig.height = 4)

# load data ---------------------------------------------------------------
path_glia <- here("Data", "Cingulotomy_brains_analysis_NeuN-GFAP-Iba1.rds")
df_glia <- read_rds(path_glia)

path_GABA <- here("Data", "Cingulotomy_brains_analysis_NeuN-GABA-GAD67.rds")
df_GABA <- read_rds(path_GABA)

path_wm <- here("Data", "White-matter-measurements.rds")
df_wm <- read_rds(path_wm) %>% dplyr::mutate(Treatment =
                                               fct_relevel(Treatment, "control"))
df_histo <- dplyr::bind_rows(df_glia, df_GABA) %>%
  dplyr::filter(Antibody != "NeuN") %>%
  dplyr::mutate(ROI = fct_collapse(ROI, VPAG = c("LPAG", "VLPAG"),
                                   DPAG = c("DMPAG", "DLPAG")))

# Modeling ----------------------------------------------------------------
df_beta_models <- df_histo %>% dplyr::nest_by(Antibody, ROI) %>%
  dplyr::rowwise(ROI, Antibody) %>%
  dplyr::mutate(model = list(glmmTMB(Prop_nonNeuN_NeuN ~ Treatment + (1|Animal_ID),
                                     data = data, family = beta_family())))


beta_models_tidy <- df_beta_models %>% reframe(broom.mixed::tidy(model)) %>%
  dplyr::filter(term == "Treatmentcingulotomy") %>%
  dplyr::select(Antibody, ROI, estimate:p.value) %>%
  dplyr::group_by(Antibody) %>%
  dplyr::mutate(fdr.p.adjusted = p.adjust(p.value, method = "fdr"),
                Holm.p.adjusted = p.adjust(p.value, method = "holm"),
                Bonferroni.p.adjusted = p.adjust(p.value, method = "bonferroni"))

beta_models_tbl <- beta_models_tidy %>% kbl() %>%
  kable_styling(bootstrap_options = c("striped", "condensed")) %>%
  row_spec(which(beta_models_tidy$p.value < 0.05), color = "red")

beta_models_tbl

# NeuN cell counts --------------------------------------------------------

path_glia <- here("Data", "Cingulotomy_brains_analysis_NeuN-GFAP-Iba1.rds")
df_glia <- read_rds(path_glia)

path_GABA <- here("Data", "Cingulotomy_brains_analysis_NeuN-GABA-GAD67.rds")
df_GABA <- read_rds(path_GABA)

df_NeuN <- dplyr::bind_rows(df_glia, df_GABA) %>%
  dplyr::filter(Antibody == "NeuN")

df_beta_models <- df_NeuN %>% dplyr::nest_by(ROI) %>%
  dplyr::mutate(model = list(glmmTMB(NeuN_cell_counts ~ Treatment + (1|Animal_ID),
                                     data = data))) %>%
  reframe(broom.mixed::tidy(model)) %>%
  dplyr::filter(term == "Treatmentcingulotomy") %>%
  dplyr::select(ROI, estimate:p.value) %>%
  dplyr::mutate(fdr.p.adjusted = p.adjust(p.value, method = "fdr"),
                Holm.p.adjusted = p.adjust(p.value, method = "holm"),
                Bonferroni.p.adjusted = p.adjust(p.value, method = "bonferroni"))
df_beta_models