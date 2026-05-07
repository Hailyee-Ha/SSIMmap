## code to prepare `DATASET` dataset goes here

library(sf)
library(terra)

Toronto_SSIM <- sf::st_read("/Users/hailyeeha/Dropbox/hailee/Western/2022FallProject/Paper_2_SSIM/New\ Writing/GIS/Toronto_cran.shp", quiet = TRUE)

fwi_0816_bc<-terra::rast("/Users/hailyeeha/Dropbox/hailee/Western/2022FallProject/Paper_2_SSIM/RCode/Revision/example\ data/fwi_0816_bc.tif")
fwi_0818_bc<-terra::rast("/Users/hailyeeha/Dropbox/hailee/Western/2022FallProject/Paper_2_SSIM/RCode/Revision/example\ data/fwi_0818_bc.tif")
fwi_1101_bc<-terra::rast("/Users/hailyeeha/Dropbox/hailee/Western/2022FallProject/Paper_2_SSIM/RCode/Revision/example\ data/fwi_1101_bc.tif")

fwi_0816_bc <- terra::wrap(fwi_0816_bc)
fwi_0818_bc <- terra::wrap(fwi_0818_bc)
fwi_1101_bc <- terra::wrap(fwi_1101_bc)

usethis::use_data(Toronto_SSIM, overwrite = TRUE)
usethis::use_data(fwi_0816_bc, overwrite = TRUE)
usethis::use_data(fwi_0818_bc, overwrite = TRUE)
usethis::use_data(fwi_1101_bc, overwrite = TRUE)


