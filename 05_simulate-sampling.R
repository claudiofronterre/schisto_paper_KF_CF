# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "dplyr") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Fake schools 
schools <- readr::read_csv("data/processed/fake_schools.csv")

# Model parameters
params_haema <- readr::read_csv("output/models_fitted/haema/fitted_parameters.csv")
params_mansoni <- readr::read_csv("output/models_fitted/mansoni/fitted_parameters.csv")


# HAEMATOBIUM ------------------------------------------------------------------

# Generate simulations scenarios

mu <- qlogis(params_haema$mu)

scenarios <- generate_scenarios(params = params_haema,
                                mu = quantile(mu, prob = seq(0, 1, l = 11)),
                                quantiles = c(0.25, 0.5, 0.75))

# Run simulations 5 x 1
results_haema <- run_scenarios(scenarios = scenarios,
                               nsims = 100,
                               schools = schools,
                               n_district = 5,
                               n_subd = 1)
# Save simulations
readr::write_csv(results_haema, "output/simulations/haema_51.csv")

# Run simulations 5 x 5
results_haema <- run_scenarios(scenarios = scenarios,
                               nsims = 100,
                               schools = schools,
                               n_district = 5,
                               n_subd = 5)
# Save simulations
readr::write_csv(results_haema, "output/simulations/haema_55.csv")

# MANSONI ----------------------------------------------------------------------

# Generate simulations scenarios

mu <- qlogis(params_mansoni$mu)

scenarios <- generate_scenarios(params = params_mansoni,
                                mu = quantile(mu, prob = seq(0, 1, l = 11)),
                                quantiles = c(0.25, 0.5, 0.75))

# Run simulations 5 x 1
results_mansoni <- run_scenarios(scenarios = scenarios,
                                 nsims = 100,
                                 schools = schools,
                                 n_district = 5,
                                 n_subd = 1)

# Save simulations
readr::write_csv(results_mansoni, "output/simulations/mansoni_51.csv")

# Run simulations  5 x 5
results_mansoni <- run_scenarios(scenarios = scenarios,
                                 nsims = 100,
                                 schools = schools,
                                 n_district = 5,
                                 n_subd = 5)

# Save simulations
readr::write_csv(results_mansoni, "output/simulations/mansoni_55.csv")
