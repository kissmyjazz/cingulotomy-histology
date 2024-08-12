
# setup -------------------------------------------------------------------
library(tidyverse)
library(here)
library(forcats)
library(ggthemes)
library(glmmTMB)
library(kableExtra)
library(ggsignif)

options(scipen=999)
set_here(path='..')
theme_set(theme_grey(base_size = 14) + theme(legend.position = "top",
                                                panel.grid = element_blank()))


cust_palette <- c("#47ABD8", "#ED1C24")


# load data ---------------------------------------------------------------
path_glia <- here("Data", "Cingulotomy_brains_analysis_NeuN-GFAP-Iba1.rds")
df_glia <- read_rds(path_glia) %>%
  dplyr::mutate(Region = ROI)

path_GABA <- here("Data", "Cingulotomy_brains_analysis_NeuN-GABA-GAD67.rds")
df_GABA <- read_rds(path_GABA) %>%
  dplyr::mutate(Region = ROI)

path_wm <- here("Data", "White-matter-measurements.rds")
df_wm <- read_rds(path_wm) %>% dplyr::mutate(Treatment =
                                               fct_recode(Treatment,
                                                          ACC = "cingulotomy",
                                                          CON = "control"),
                                             Treatment =
                                               fct_relevel(Treatment, "CON"))

# aggregate GABA data
df_GABA_r <- df_GABA %>% dplyr::select(Animal_ID, Treatment, Region, Antibody,
                                       Slide, Section, Num_detections,
                                       Area_sq_mm) %>%
  dplyr::arrange(Animal_ID, Slide, Section, Region, Antibody) %>%
  dplyr::mutate(NeuN_detections = ifelse(Antibody == "NeuN",
                                         Num_detections, NA_real_)) %>%
  tidyr::fill(NeuN_detections) %>%
  dplyr::mutate(Region = dplyr::recode_factor(Region,
                                              Area_24 = "ACC",
                                              Area_23 = "ACC", DLPAG = "DPAG",
                                              DMPAG = "DPAG", LPAG = "VPAG", VLPAG = "VPAG",
                                              .default = levels(Region)),
                Treatment = fct_recode(Treatment,
                                       ACC = "cingulotomy",
                                       CON = "control")) %>%
  dplyr::group_by(Animal_ID, Treatment, Region, Antibody) %>%
  dplyr::summarise(Cells_per_sq_mm = sum(Num_detections) / sum(Area_sq_mm),
                   NeuN_detections = first(ifelse(Antibody == "NeuN",
                                                  sum(NeuN_detections) / sum(Area_sq_mm),
                                                  NA_real_))) %>%
  dplyr::arrange(Animal_ID, Region, Antibody) %>%
  tidyr::fill(NeuN_detections) %>% ungroup() %>%
  dplyr::mutate(Prop_nonNeuN_NeuN = Cells_per_sq_mm / (Cells_per_sq_mm + NeuN_detections))

# preparation for the graph
df_GABA_r2 <- df_GABA_r %>%
  dplyr::filter(Region != "ACC") %>%
  dplyr::mutate(Region = fct_relevel(Region, "Amy_Ce", "Amy_AB", "Amy_B", "Amy_L", "DPAG", "VPAG"))

# df for significance mapping
df_GABA_graph <- df_GABA %>%
  dplyr::mutate(Region = dplyr::recode_factor(Region,
                                              Area_24 = "ACC",
                                              Area_23 = "ACC", DLPAG = "DPAG",
                                              DMPAG = "DPAG", LPAG = "VPAG", VLPAG = "VPAG",
                                              .default = levels(Region)),
                Treatment = fct_recode(Treatment,
                                       ACC = "cingulotomy",
                                       CON = "control")) %>%
  dplyr::filter(Antibody == "GABA", Region != "ACC") %>%
  dplyr::mutate(Region = fct_relevel(Region, "Amy_Ce", "Amy_AB", "Amy_B", "Amy_L", "DPAG", "VPAG"))


# Modeling ----------------------------------------------------------------
df_beta_models <- df_GABA_r2 %>% dplyr::nest_by(Region) %>%
  dplyr::mutate(model = list(glmmTMB(Prop_nonNeuN_NeuN ~ Treatment,
                                     data = data, family = beta_family())))


beta_models_tidy <- df_beta_models %>% reframe(broom.mixed::tidy(model)) %>%
  dplyr::filter(term == "TreatmentACC") %>%
  dplyr::select(Region, estimate:p.value) %>%
  dplyr::mutate(fdr.p.adjusted = p.adjust(p.value, method = "fdr"),
                Holm.p.adjusted = p.adjust(p.value, method = "holm"),
                Bonferroni.p.adjusted = p.adjust(p.value, method = "bonferroni"))

beta_models_tbl <- beta_models_tidy %>% kbl() %>%
  kable_styling(bootstrap_options = c("striped", "condensed")) %>%
  row_spec(which(beta_models_tidy$p.value < 0.05), color = "red")

beta_models_tbl


# p-values ----------------------------------------------------------------

p_values <- df_GABA_graph %>%
  dplyr::filter(Antibody == "GABA") %>%
  dplyr::nest_by(Region) %>%
  dplyr::mutate(model = list(glmmTMB(Prop_nonNeuN_NeuN ~ Treatment + (1|Animal_ID),
                                     data = data, family = beta_family()))) %>%
  reframe(broom.mixed::tidy(model)) %>%
  dplyr::filter(term == "TreatmentACC") %>%
  dplyr::select(p.value) %>%
  dplyr::pull() %>% as.list()

p_values

# graphs ------------------------------------------------------------------

# list of facet labels
nuclei <- c(
  'DPAG' = "Dorsal\nperiaqueductal gray",
  'VPAG' = "Ventral\nperiaqueductal gray",
  'Amy_AB' = "Basomedial amygdala",
  'Amy_B' = "Basolateral amygdala",
  'Amy_Ce' = "Central amygdala",
  'Amy_L' = "Lateral amygdala"
)

g_GABA <- ggplot(df_GABA_r2 %>% dplyr::filter(Antibody == "GABA"),
                 aes(x = Treatment,
                     y = Prop_nonNeuN_NeuN)) +
  stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. - 0.3,
                                                             yend = ..y..),
               color = "gray10", linewidth = 1.6) +
  stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = ..x.. + 0.3,
                                                             yend = ..y..),
               color = "gray10", linewidth = 1.6) +
  geom_point(aes (colour = Treatment, shape = Animal_ID), size = 4, position = position_jitter(width = 0.0)) +
  scale_colour_manual(values = cust_palette) +
  scale_shape_manual(values = c(0:3, 15:18), guide = "none") +
  facet_wrap(~Region, nrow = 1, labeller = as_labeller(nuclei)) +
  labs(y = "Proportion of GABA-positive detections", x = NULL) +
  theme(strip.text = element_text(size = 12))
g_GABA

# the number of significance stars over dorsal PAG graph is needed to be adjusted manually
g <- g_GABA + geom_signif(data = df_GABA_graph,
                          comparisons = list(c("CON", "ACC")),
                          test = "wilcox.test",
                          map_signif_level = TRUE,
                          y_position = 0.54,
                          tip_length = 0.02,
                          color = "black")
g

ggsave(file = here("Figures", "gaba-graph.svg"), plot = g, width = 12,
       height = 6)

