# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sf", "dplyr") # package names
pacman::p_load(pkgs, character.only = T)

# LOAD DATA --------------------------------------------------------------------

# Fake districts 
districts <- st_read("data/processed/geodata/fake_country.gpkg", layer = "districts")

subd <- st_read("data/processed/geodata/fake_country.gpkg", layer = "sub_districts")

# Ius 
ius <- st_read("data/original/geodata/AFRO_IUs_201802.shp") %>% 
  st_transform(crs = "+init=epsg:3857 +units=km")

# SIMULATE SCHOOLS LOCATIONS ---------------------------------------------------

# Simulate schools according to the density of schools per area in Uganda
nschools_uganda <- 14000
area_uganda <- as.numeric(sum(st_area(ius[ius$ADMIN0 == "Uganda", ])))
nschools <- round(sum((nschools_uganda / area_uganda) * as.numeric(st_area(districts))))
nschools <- 15000
schools <- st_as_sf(st_sample(districts, size = nschools))

plot(districts)
plot(schools, add = T, cex = .1)

# Check the average distance between schools
coords <- st_coordinates(schools)
summary(spatstat::nndist(coords)) 

# Add ID to districts and subd
districts$ID <- 1:nrow(districts)
subd$ID <- 1:nrow(subd)

# Associate each school with a district and a sub district
schools$ID_district <- st_join(schools, districts)$ID 
schools$ID_subd <- st_join(schools, subd)$ID 

# Check that all districts have at least 5 schools
sum(table(schools$ID_district) < 5) 

# Check that all sub-districts have at least 1 school
sum(table(schools$ID_subd) < 5) 

# Create data-set with schools coordinates and IDs and save it
schools <- st_drop_geometry(schools)
schools$x <- coords[, 1]
schools$y <- coords[, 2]

readr::write_csv(schools, "data/processed/fake_schools.csv")
