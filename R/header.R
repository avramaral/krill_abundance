add_args <- commandArgs()
res_pred <- add_args[1]
if ((is.null(res_pred)) | (res_pred == "RStudio")) { res_pred <- 0.01 }
print(paste("The prediction resolution is: ", res_pred, ".", sep = ""))

suppressMessages(library("tidyverse"))
suppressMessages(library("patchwork"))
suppressMessages(library("cowplot"))
suppressMessages(library("viridis"))
suppressMessages(library("raster"))
suppressMessages(library("sf"))
suppressMessages(library("ggthemes"))
suppressMessages(library("gstat"))
suppressMessages(library("dismo"))
suppressMessages(library("ncdf4"))
suppressMessages(library("marmap"))
suppressMessages(library("stars"))
suppressMessages(library("terra"))
suppressMessages(library("spatialEco"))
suppressMessages(library("waveslim"))

suppressMessages(library("INLA"))

pal <- c("#00008FFF", "#0000F2FF", "#0063FFFF", "#00D4FFFF", "#46FFB8FF", "#B8FF46FF", "#FFD400FF", "#FF6300FF", "#F00000FF", "#800000FF")

# Read 48.3 shape
subarea_48_3_file <- "DATA/SHAPEFILES/subarea_48_3.RDS"
if (!file.exists(subarea_48_3_file)) {
  shape <- st_read(dsn = "DATA/SHAPEFILES/ssmus/", layer = "ssmusPolygon", quiet = TRUE)
  subarea_48_3 <- shape %>% dplyr::filter(grepl("^SSMU 48.3", GAR_Name))
  saveRDS(object = subarea_48_3, file = subarea_48_3_file)
} else {
  subarea_48_3 <- readRDS(file = subarea_48_3_file)
}

coord_box <- st_bbox(subarea_48_3)
LON <- unname(c(coord_box["xmin"], coord_box["xmax"]))
LAT <- unname(c(coord_box["ymin"], coord_box["ymax"]))

# Coordinate limits based on the area plotted in the initial proposal
LON_sub <- c(-42, -36)
LAT_sub <- c(-55, -52)

# Coordinate limits based on the area with acoustic data
LON_sub_sub <- c(-40, -37.5)
LAT_sub_sub <- c(-54.1, -53)

# Coordinates limits based on the area with intersection with the KRILLBASE
LON_krillbase <- c(-40, -35) # c(-40, -35)
LAT_krillbase <- c(-56, -53)

###############################
# Create boxes for study area #
###############################

buff <- 0.5
expanded_study_area <- matrix(c(LON[1] - buff, LON[1] - buff, LON[2] + buff, LON[2] + buff, LON[1] - buff, LAT[1] - buff, LAT[2] + buff, LAT[2] + buff, LAT[1] - buff, LAT[1] - buff), ncol = 2, byrow = FALSE) %>% st_linestring()
expanded_study_area <- expanded_study_area %>% st_line_sample(n = 10e3)
expanded_study_area <- expanded_study_area %>% st_cast("POLYGON") %>% st_sfc(crs = 4326)

study_area <- st_multipoint(c(st_point(c(LON[1], LAT[1])), 
                              st_point(c(LON[1], LAT[2])),
                              st_point(c(LON[2], LAT[2])),
                              st_point(c(LON[2], LAT[1])))) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326)

study_area_sub <- st_multipoint(c(st_point(c(LON_sub[1], LAT_sub[1])), 
                                  st_point(c(LON_sub[1], LAT_sub[2])),
                                  st_point(c(LON_sub[2], LAT_sub[2])),
                                  st_point(c(LON_sub[2], LAT_sub[1])))) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326)

study_area_sub_sub <- st_multipoint(c(st_point(c(LON_sub_sub[1], LAT_sub_sub[1])), 
                                      st_point(c(LON_sub_sub[1], LAT_sub_sub[2])),
                                      st_point(c(LON_sub_sub[2], LAT_sub_sub[2])),
                                      st_point(c(LON_sub_sub[2], LAT_sub_sub[1])))) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326)

study_area_krillbase <- st_multipoint(c(st_point(c(LON_krillbase[1], LAT_krillbase[1])), 
                                        st_point(c(LON_krillbase[1], LAT_krillbase[2])),
                                        st_point(c(LON_krillbase[2], LAT_krillbase[2])),
                                        st_point(c(LON_krillbase[2], LAT_krillbase[1])))) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326)

# Retrieve islands shapes
land <- st_read(dsn = "DATA/SHAPEFILES/coastline/", layer = "coastline", quiet = TRUE)
land <- st_transform(land, crs = 4326) 
land <- land$geometry
land <- st_intersection(land, st_buffer(study_area, 1))
land <- land %>% st_as_sf() %>% st_cast("POLYGON") %>% mutate(area = st_area(.)) %>% arrange(desc(area)) %>% slice(1) %>% dplyr::select(x)

tmp_subarea_48_3 <- subarea_48_3[1, ] %>% st_make_valid() %>% st_union() %>% st_sf()
tmp_subarea_48_3 <- st_cast(tmp_subarea_48_3, "POLYGON") %>% st_geometry() %>% st_cast("LINESTRING") %>% st_cast("POLYGON") %>% st_sf()
tmp_subarea_48_3 <- tmp_subarea_48_3[1, ]

water     <- st_difference(expanded_study_area, land)
water_alt <- st_difference(expanded_study_area, land)

water_sub       <- (st_intersection(water, study_area_sub)       %>% st_cast("POLYGON") %>%  {.[which.max(st_area(.))]} %>% st_sfc() %>% st_cast("LINESTRING"))[[1]] %>% st_line_sample(n = 10e3) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326)
water_sub_sub   <- (st_intersection(water, study_area_sub_sub)   %>% st_cast("POLYGON") %>%  {.[which.max(st_area(.))]} %>% st_sfc() %>% st_cast("LINESTRING"))[[1]] %>% st_line_sample(n = 10e3) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326) # Select the largest spatial object and convert it to `LINESTRING` so the boundaries are not a box
water_krillbase <- (st_intersection(water, study_area_krillbase) %>% st_cast("POLYGON") %>%  {.[which.max(st_area(.))]} %>% st_sfc() %>% st_cast("LINESTRING"))[[1]] %>% st_line_sample(n = 10e3) %>% st_cast("POLYGON") %>% st_sfc(crs = 4326) # Select the largest spatial object and convert it to `LINESTRING` so the boundaries are not a box

# Coordinates for prediction
res_pred <- res_pred # c(0.01, 0.5)
coord_pred_file <- paste("DATA/OTHERS/coord_pred_res_",  gsub("\\.", "_", as.character(res_pred)), ".RDS", sep = "")
if (!file.exists(coord_pred_file)) {
  pts_bdy_x <- LON + c(-0.5, 0.5)
  pts_bdy_y <- LAT + c(-0.5, 0.5)
  coord_pred <- expand.grid(x = seq(pts_bdy_x[1], pts_bdy_x[2], by = res_pred), y = seq(pts_bdy_y[1], pts_bdy_y[2], by = res_pred))
  colnames(coord_pred) <- c("longitude", "latitude")
  
  coord_pred_sf <- st_as_sf(coord_pred, coords = c("longitude", "latitude"), crs = st_crs(land))
  coord_pred_sf <- coord_pred_sf %>% mutate(intersection = as.integer(st_intersects(geometry, land))) 
  coord_pred_sf <- coord_pred_sf %>% dplyr::filter(is.na(intersection))
  
  coord_pred_filtered <- st_coordinates(coord_pred_sf)
  colnames(coord_pred_filtered) <- c("longitude", "latitude")
  
  coord_pred <- as.matrix(coord_pred)
  coord_pred_filtered <- as.matrix(coord_pred_filtered)
  
  coord_pred_all <- list(coord_pred = coord_pred, coord_pred_filtered = coord_pred_filtered)
  saveRDS(object = coord_pred_all, file = coord_pred_file)
} else {
  coord_pred_all <- readRDS(file = coord_pred_file)
  coord_pred <- coord_pred_all$coord_pred
  coord_pred_filtered <- coord_pred_all$coord_pred_filtered
}


