# LOAD REQUIRED PACKAGES AND FUNCTIONS -----------------------------------------
if (!require("pacman")) install.packages("pacman")
pkgs = c("sp", "rgdal", "sf", "dplyr", "raster") # package names
pacman::p_load(pkgs, character.only = T)

# LOAD DATA --------------------------------------------------------------------

# Load IU data to extract boundaries and other values 
ius <- st_read("data/original/geodata/AFRO_IUs_201802.shp") %>% 
  st_transform(crs = "+init=epsg:3857 +units=km")


# GENERATE FAKE COUNTRY WITH DISTRICT AND SUB - DISTRICTS ----------------------

# For each country calculate:
#   - the average size (area in km2) of the IUs
#   - the total number of IUs
#   - the total area of the country
countries <- as.character(unique(ius$ADMIN0))

ius_summaries <- data.frame(country = rep(0, length(countries)), 
                            avg_size_iu = rep(0, length(countries)),
                            total_iu = rep(0, length(countries)),
                            total_area = rep(0, length(countries)))

for (i in 1:length(countries)) {
  ius %>% 
    filter(ADMIN0 == countries[i]) %>% 
    st_area() -> areaIU
  ius_summaries[i, 2] <- mean(areaIU)
  ius_summaries[i, 4] <- sum(areaIU)
  ius_summaries[i, 3] <- length(areaIU)
  ius_summaries[i, 1] <- countries[i]
}


# Generate fake country inside Africa with area = median_country_size
generate_country <- function(xmin, xmax, ymin, ymax, crs) {
  box <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons")
  proj4string(box) <- crs(crs)
  return(box)
}

median_country_size <- median(ius_summaries$total_area)

st_bbox(ius)
sqrt(median_country_size)
mycountry <- generate_country(xmin = 2100, xmax = 2793, ymin = 0, ymax = 400, 
                              crs = "+init=epsg:3857 +units=km")

area(mycountry) / median_country_size # Should be close to 1

plot(ius$geometry)
plot(mycountry, add = T, col = "red")
plot(mycountry)

# Generate number of districts based on iu to country ratio 
iu_country_ratio <- median(ius_summaries$total_iu / ius_summaries$total_area)
nius <- round(iu_country_ratio * area(mycountry))
nius <- 70

# Instead of random sampling use an inhibitory spatial point pattern fixing the minimum distance
# if you want to have sub-districts of approximately the same size
# pp <- spatstat::rMaternI(kappa = (nius * 3) / area(mycountry),
#                          r = 20, win = maptools::as.owin.SpatialPolygons(mycountry))

# First simulate the subdistrics then simulate the districts
# from the union of subdistrics

factor <- 5 # Sub-district are factor times more than districts
pp <- spsample(mycountry, n = nius * factor, type = "random")
crs(pp) <- crs(mycountry)

plot(mycountry)
points(pp)
subd <- dismo::voronoi(pp)
subd <- rgeos::intersect(subd, mycountry)
plot(subd)

# Apply geometry union to every factor sub-districs to generate the disctricts (IUs)
# coords <- coordinates(subd)
# rownames(coords) <- 1:length(subd)
# id <- NULL
# for (i in 1:nius) {
#   idTemp <- order(pointDistance(coords[1, ], coords, lonlat = F))[2:factor]
#   id <- as.numeric(c(id, rownames(coords)[1], rownames(coords)[idTemp]))
#   coords <- coords[-c(1, idTemp), ]
# }
# 
coords <- coordinates(subd)
# kclust <- kmeans(coords, centers = nius)
id <- cutree(hclust(dist(coords), method = "ward.D2"), k = nius)
t <- as.numeric(table(id))
mean((5 - t)^2)
districts <- rgeos::gUnaryUnion(spgeom = subd, id = id)



subd <- subd[id, ]
districts <- rgeos::gUnaryUnion(spgeom = subd, id = rep(1:nius, each = factor))
plot(subd)
plot(districts)
subd$id <- NULL

# Save districts and sub-districts in a geopackage
st_write(st_as_sf(districts),
         delete_dsn = T,
         dsn = "data/processed/geodata/fake_country.gpkg",
         layer = "districts")


st_write(st_as_sf(subd),
         dsn = "data/processed/geodata/fake_country.gpkg",
         layer = "sub_districts")
