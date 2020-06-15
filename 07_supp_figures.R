# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!re1quire("pacman")) install.packages("pacman")
pkgs = c("tmap", "sf") # package names
pacman::p_load(pkgs, character.only = T)

source("R/functions.R")

# LOAD DATA --------------------------------------------------------------------

# Prevalence surveys
haema <- readr::read_csv("data/processed/haema_for_models.csv")
mansoni <- readr::read_csv("data/processed/mansoni_for_models.csv")

# Shapefiles
africa <- st_read("data/processed/geodata/africa.gpkg")

# Fake districts 
districts <- st_read("data/processed/geodata/fake_country.gpkg", layer = "districts")
subd <- st_read("data/processed/geodata/fake_country.gpkg", layer = "sub_districts")
schools <- readr::read_csv("data/processed/fake_schools.csv")


# DATA PROCESSING --------------------------------------------------------------
haema$specie <- "Haematobium"
mansoni$specie <- "Mansoni"

schisto <- dplyr::bind_rows(haema, mansoni)
schisto_sp <- schisto %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

schools_sp <- schools %>% 
  st_as_sf(coords = c("x", "y"), crs = st_crs(districts))

schisto_sp$prev100 <- (schisto_sp$positive / schisto_sp$examined) * 100

params_haema <- readr::read_csv("output/models_fitted/haema/fitted_parameters.csv")


# MAPS -------------------------------------------------------------------------

# Survey locations by specie
t1 <- tm_shape(africa) +
  tm_fill("white") + 
tm_shape(schisto_sp) +
  tm_symbols(col = "specie", legend.col.show = F,
             size = .007, palette = c("#CC0C00FF", "#5C88DAFF"), 
             shape = 19, alpha = .5) +
tm_shape(africa) +
  tm_borders("black", lwd = .5) +
tm_add_legend(type = "symbol", labels = c("Haematobium", "Mansoni"), 
              col = c("#CC0C00FF", "#5C88DAFF"), shape = 19,
              size = .2, title = "Specie", 
              alpha = .6) +
tm_compass(position = "right") +
tm_layout(bg.color = "lightblue", design.mode = F, legend.bg.color = "white", 
          asp = 0,
          legend.frame = "black", scale = 1, legend.position = c("right", "top"),
          frame.double.line = F, legend.title.size = .95, outer.margins = 0) +
tm_scale_bar(position = "left") 

tmap_save(t1, filename = "figs/species_locations.png",
          width = 6.3, height = 7)


# Spatial distribution of prevalence 
n <- 8
prev100 <- schisto_sp$prev100[schisto_sp$prev100 > 0]
breaks <- as.numeric(quantile(prev100, probs = seq(0, 1, l = n + 1)))
digits <- 0
breaks[2:n] <- round(breaks[2:n], digits)
breaks[1] <- floor_dec(breaks[1], digits)
breaks[n + 1] <- ceiling_dec(breaks[n + 1], digits) 


# Generate labels for color legend
labs <- create_labels(breaks, greater = T, text_sep = " - ")


# Generate color palette
pal <- tmaptools::get_brewer_pal("-RdYlBu", n = length(labs), contrast = c(0, 1), plot = F)

t2 <- tm_shape(africa) +
  tm_fill("white") +
tm_shape(schisto_sp) +
  tm_symbols(col = "prev100", legend.col.show = F,  
             size = .005, shape = 20, alpha = .5, palette = pal,
             border.col = "black", border.lwd = 0.5,
             style = "fixed", breaks = breaks) +
  tm_add_legend(type = "symbol", labels = labs,
                shape = 19, col = pal, size = .3, alpha = .5, 
                title = "Prevalence (%)",
                is.portrait = T) +
  tm_facets(by = "specie", free.coords = F, ncol = 2, nrow = 1) +
tm_shape(africa) +
  tm_borders("black", lwd = .25) +
tm_compass(position = "right", size = 1.5, text.size = .7) +
tm_layout(bg.color = "lightblue", design.mode = F, legend.bg.color = "white", 
          asp = 0,
          legend.frame = "black", scale = 1,
          frame.double.line = F, legend.title.size = .95, outer.margins = 0) +
tm_scale_bar(position = "left", text.size = .3) 

t2 

tmap_save(t2, filename = "figs/species_prev.png",
          width = 10, height = 4)
        
# Fake country 
fc <- tm_shape(subd) +
  tm_polygons(border.col = "gray20", col = "gray85", lwd = 1.3) +
tm_shape(districts) +
  tm_borders(col = "black", lwd = 2) +
tm_add_legend(type = "fill",
              labels = "Districts",
              col = "white",
              border.col = "black", 
              lwd = 3) +
  tm_add_legend(type = "fill",
                labels = "Sub-districts",
                col = "gray85",
                border.col = "gray20") +
  tm_add_legend(type = "symbol",
                labels = "Schools",
                col = "black",
                shape = 19) +
tm_shape(schools_sp) +
  tm_symbols(col = "black", shape = 19, size = .005) +
tm_scale_bar() +
tm_layout(asp = 0, frame = T, legend.position = c("left", "top"),
          legend.stack = "horizontal")

fc

tmap_save(fc, filename = "figs/simulated_country.png",
          width = 7, height = 4.7)

params_haema <- readr::read_csv("output/models_fitted/haema/fitted_parameters.csv")
