# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("ggplot2", "hrbrthemes", "dplyr") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------
haema55 <- readr::read_csv("output/simulations/haema_55.csv")
haema51 <- readr::read_csv("output/simulations/haema_51.csv")

mansoni55 <- readr::read_csv("output/simulations/mansoni_55.csv")
mansoni51 <- readr::read_csv("output/simulations/mansoni_51.csv")

# HAEMATOBIUM ------------------------------------------------------------------
haema55$design <- "5 - 5"
haema51$design <- "5 - 1"
haema <- dplyr::bind_rows(haema55, haema51)
haema[, 5:22] <- haema[, 5:22] * 100 

# District vs Sub-district 

# Compare the proportion of districts correctly treated
prevalence <- plot_results(haema, x = "correct_school_d", y = "correct_school_subd", 
                           col = "prev_sim", 
                           lx = 7, ly = 7, lcol = 5,
                           xlab = "District",
                           ylab = "Sub-district") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  facet_wrap(~ design) + 
  labs(title = "Proportion of schools correctly classified (%)") +
  coord_equal() + 
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

prevalence
ggsave(prevalence, filename = "figs/combined/haema/schools.png",
       width = 9, height = 4.5, dpi = 600)

overtreated <- plot_results(haema, x = "overtreated_d", y = "overtreated_subd", 
                           col = "prev_sim", 
                           lx = 5, ly = 5, lcol = 5,
                           xlab = "District",
                           ylab = "Sub-district") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  facet_wrap(~ design) + 
  labs(title = "Proportion of schools overtreated (%)") +
  coord_equal() + 
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

overtreated
ggsave(overtreated, filename = "figs/combined/haema/overtreated.png",
       width = 9, height = 4.5, dpi = 600)

undertreated <- plot_results(haema, x = "undertreated_d", y = "undertreated_subd", 
                             col = "prev_sim", 
                             lx = 5, ly = 5, lcol = 5,
                             xlab = "District",
                             ylab = "Sub-district") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  facet_wrap(~ design) + 
  labs(title = "Proportion of schools undertreated (%)") +
  coord_equal() + 
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

undertreated
ggsave(undertreated, filename = "figs/combined/haema/undertreated.png",
       width = 9, height = 5, dpi = 600)

# ROC curves

# MANSONI ----------------------------------------------------------------------
mansoni55$design <- "5 - 5"
mansoni51$design <- "5 - 1"
mansoni <- dplyr::bind_rows(mansoni55, mansoni51)
mansoni[, 5:22] <- mansoni[, 5:22] * 100 

# District vs Sub-district 

# Compare the proportion of districts correctly treated
prevalence <- plot_results(mansoni, x = "correct_school_d", y = "correct_school_subd", 
                           col = "prev_sim", 
                           lx = 7, ly = 7, lcol = 5,
                           xlab = "District",
                           ylab = "Sub-district") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  facet_wrap(~ design) + 
  labs(title = "Proportion of schools correctly classified (%)") +
  coord_equal() + 
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

prevalence

ggsave(prevalence, filename = "figs/combined/mansoni/schools.png",
       width = 9, height = 4.5, dpi = 600)

overtreated <- plot_results(mansoni, x = "overtreated_d", y = "overtreated_subd", 
                            col = "prev_sim", 
                            lx = 5, ly = 5, lcol = 5,
                            xlab = "District",
                            ylab = "Sub-district") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  facet_wrap(~ design) + 
  labs(title = "Proportion of schools overtreated (%)") +
  coord_equal() + 
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

overtreated
ggsave(overtreated, filename = "figs/combined/mansoni/overtreated.png",
       width = 9, height = 4.5, dpi = 600)

undertreated <- plot_results(mansoni, x = "undertreated_d", y = "undertreated_subd", 
                             col = "prev_sim", 
                             lx = 5, ly = 5, lcol = 5,
                             xlab = "District",
                             ylab = "Sub-district") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 1) +
  facet_wrap(~ design) + 
  labs(title = "Proportion of schools undertreated (%)") +
  coord_equal() + 
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

undertreated
ggsave(undertreated, filename = "figs/combined/mansoni/undertreated.png",
       width = 9, height = 5, dpi = 600)

