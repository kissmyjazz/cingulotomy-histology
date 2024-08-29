library(here)
library(tidyverse)
library(ggdark)
library(lme4)
library(glmmTMB)
library(broom)
library(broom.mixed)
library(ggsignif)


theme_set(theme_grey(base_size = 14) + theme(legend.position = "top",
                                             panel.grid = element_blank()))


cust_palette <- c("#47ABD8", "#ED1C24")

# Define the folder path
folder_path <- here("Data", "cellpose_results")  # Replace with your actual folder path

# Get all excel files in the folder
files <- list.files(path = folder_path, pattern = "^[^~]*\\.tsv$",
                          full.names = TRUE)

named_files <- purrr::set_names(files, nm = basename(files))

df_ <- named_files |>
  map_dfr(read_tsv, .id = NULL, col_types = "c_________iiid_", name_repair = "universal")

# data frame for graphs -----------------------------------------------------
df <- df_ |> dplyr::mutate(Treatment = dplyr::case_when(str_detect(Image, "Chioma") ~ "CON",
                                                        str_detect(Image, "Keala") ~ "CON",
                                                        str_detect(Image, "Lehua") ~ "ACC",
                                                        str_detect(Image, "Lovita") ~ "ACC",
                                                        str_detect(Image, "Rashmi") ~ "ACC",
                                                        str_detect(Image, "Shashi") ~ "CON",
                                                        str_detect(Image, "Taiwo") ~ "CON",
                                                        str_detect(Image, "Yazhu") ~ "ACC",
                                                        .default = NA), Treatment = factor(Treatment, levels = c("CON", "ACC")),
                           ID = dplyr::case_when(str_detect(Image, "Chioma") ~ "Chioma",
                                                 str_detect(Image, "Keala") ~ "Keala",
                                                 str_detect(Image, "Lehua") ~ "Lehua",
                                                 str_detect(Image, "Lovita") ~ "Lovita",
                                                 str_detect(Image, "Rashmi") ~ "Rashmi",
                                                 str_detect(Image, "Shashi") ~ "Shashi",
                                                 str_detect(Image, "Taiwo") ~ "Taiwo",
                                                 str_detect(Image, "Yazhu") ~ "Yazhu",
                                                 .default = NA), Image = NULL) |>
  dplyr::summarise(across(everything(), ~mean(.x)), .by = c(Treatment, ID)) |>
  pivot_longer(cols = -c(ID, Treatment, Area.µm.2), names_to = "Cell_type", values_to = "Count") |>
  dplyr::mutate(Cell_type = str_remove(Cell_type, "Num."), NeuN = if_else(str_detect(Cell_type, "NeuN"),
                                                                          Count, NA_integer_)) |>
               tidyr::fill(NeuN, .direction = "up") |>
  dplyr::mutate(Ratio = Count / (Count + NeuN)) |>
  dplyr::mutate(Cell_type = fct_relevel(Cell_type, "NeuN", "GFAP", "Iba1"),
                Group = case_when(ID %in% c("Shashi", "Keala") ~ "Sham",
                                  ID %in% c("Taiwo", "Chioma") ~ "Control",
                                  ID %in% c("Lehua", "Lovita", "Rashmi", "Yazhu") ~ "Cingulotomy",
                                  TRUE ~ NA_character_, .ptype = "factor"))
# data frame for p-values -------------------------------------------------
df_p <- df_ |> dplyr::mutate(Treatment = dplyr::case_when(str_detect(Image, "Chioma") ~ "CON",
                                                        str_detect(Image, "Keala") ~ "CON",
                                                        str_detect(Image, "Lehua") ~ "ACC",
                                                        str_detect(Image, "Lovita") ~ "ACC",
                                                        str_detect(Image, "Rashmi") ~ "ACC",
                                                        str_detect(Image, "Shashi") ~ "CON",
                                                        str_detect(Image, "Taiwo") ~ "CON",
                                                        str_detect(Image, "Yazhu") ~ "ACC",
                                                        .default = NA), Treatment = factor(Treatment, levels = c("CON", "ACC")),
                           ID = dplyr::case_when(str_detect(Image, "Chioma") ~ "Chioma",
                                                 str_detect(Image, "Keala") ~ "Keala",
                                                 str_detect(Image, "Lehua") ~ "Lehua",
                                                 str_detect(Image, "Lovita") ~ "Lovita",
                                                 str_detect(Image, "Rashmi") ~ "Rashmi",
                                                 str_detect(Image, "Shashi") ~ "Shashi",
                                                 str_detect(Image, "Taiwo") ~ "Taiwo",
                                                 str_detect(Image, "Yazhu") ~ "Yazhu",
                                                 .default = NA), Image = NULL) |>
  pivot_longer(cols = -c(ID, Treatment, Area.µm.2), names_to = "Cell_type", values_to = "Count") |>
  dplyr::mutate(Cell_type = str_remove(Cell_type, "Num.")) |>
  dplyr::mutate(Cell_type = fct_relevel(Cell_type, "NeuN", "GFAP", "Iba1"))


# p-values ----------------------------------------------------------------
p_values <- df_p %>%
  dplyr::nest_by(Cell_type) %>%
  dplyr::mutate(model = list(glmmTMB(Count ~ Treatment + (1|ID),
                                     data = data))) %>%
  reframe(broom.mixed::tidy(model)) %>%
  dplyr::filter(term == "TreatmentCON") %>%
  dplyr::select(p.value) %>%
  dplyr::pull() %>% as.list()

p_values

gg_counts <- ggplot(df,
                   aes(x = Treatment,
                       y = Count)) +
  stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = after_stat(x) - 0.3,
                                                             yend = after_stat(y)),
               color = "gray10", linewidth = 1.6) +
  stat_summary(fun = "mean", geom = "segment", mapping = aes(xend = after_stat(x) + 0.3,
                                                             yend = after_stat(y)),
               color = "gray10", linewidth = 1.6) +
  geom_point(aes(colour = Treatment, shape = ID), size = 4, position = position_jitter(width = 0.0)) +
  scale_colour_manual(values = cust_palette) +
  scale_shape_manual(values = c(0:3, 15:18), guide = "none") +
  facet_wrap(~Cell_type, nrow = 1, scales = "free_y") +
  labs(y = "Cell counts", x = NULL)

gg_counts

gg <- gg_counts + geom_signif(data = df_p,
                          comparisons = list(c("CON", "ACC")),
                          test = "t.test",
                          map_signif_level = TRUE,
                          tip_length = 0.02,
                          color = "black")
gg

ggsave(file = here("Figures", "lesion-graph.svg"), plot = gg, width = 6,
       height = 5)
