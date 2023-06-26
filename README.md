# SSIMmap

<!-- badges: start -->

<!-- badges: end -->

The goal of the *SSIMmap* package is to offer a comprehensive set of tools for working with Structural Similarity Index (SSIM) for map comparisons.The original SSIM develop to compare two image mimicking the human visual system, which is presented in the [paper](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1284395). *SSIMmap* package applied this method to maps (polygon and raster). The package includes four key functions for calculating the similarities between maps: **ssim_bandwidth**, **ssim_constant**, **ssim_polygon**,**ssim_raster**. Users who want to compare two maps and quantify their similarities can utilize this package and visualize the results using other tools(e.g.,[*tmap*](https://cran.r-project.org/web/packages/tmap/index.html) and [*ggplot*](https://cran.r-project.org/web/packages/ggplot2/index.html)).

## Installation

You can install the development version of SSIMmap from [GitHub](https://github.com/Hailyee-Ha/SSIMmap.git) with the following commands:

``` r
install.packages("devtools")
devtools::install_github("Hailyee-Ha/SSIMmap")
```

## Usage

To use the SSIMmap package, first load it into your R environment:

``` r
library(SSIMmap)
```

## Functions

SSIMmap includes the following key functions:

1)  **ssim_bandwidth**: This function calculates the bandwidth size for the computation of the SSIM on polygon maps. It offers two options for selecting the bandwidth size from two methods:1. square root of N and 2. the best trade-off between bias and variance. The function takes as input a shape file in *sf* class including columns for the SSIM calculation and returns the two options for selecting bandwidth sizes on the basis of two methods. It also provides a plot illustrating the bias and variance against the size of bandwidth.

2)  **ssim_constant**: This function calculates constants(k1 and k2) for the computation of the SSIM on polygon maps. It takes as input a shape filein *sf* class including columns for the SSIM calculation and returns the constants on the console window.

3)  **ssim_polygon**: This function calculates the SSIM index for a given polygon. It takes as input a shape file in *sf* class including columns for the SSIM calculation and returns either the global SSIM values (global=TRUE) or the SSIM values for each given polygon as the local SSIM (global=FALSE).

4)  **ssim_raster**: This function calculates the SSIM index for raster images. It takes as input a image file importing from the [*terra*](https://cran.r-project.org/web/packages/terra/index.html) and returns either the global SSIM values (global=TRUE) or the SSIM values for each given cell as the local SSIM (global=FALSE).

## Example

Here is a basic example which shows you how to use the package:

``` r
library(SSIMmap)

## basic example code
##example polygon map of Toronto including columns(1. Pampalon Index("PP_SDD"), 2. CIMD Index("CIMD_SDD"), and 3. the percent of the housholds who commute within the same census subdivision("P_commute"))

shape<-SSIMmap::polygon 

##finding the options for selecting the bandwidth sizes from two methods
ssim_bandwidth(shape, "PP_SDD","CIMD_SDD", max_bandwidth=100)
ssim_bandwidth(shape, "P_commute","CIMD_SDD", max_bandwidth=100)
````


##finding the constants(k1 and k2) for the SSIM calculation
ssim_constant(shape,"PP_SDD","CIMD_SDD")
ssim_constant(shape,"P_commute","CIMD_SDD")

##finding the global SSIM for the polygon map
ssim_polygon(shape,"PP_SDD","CIMD_SDD")
ssim_polygon(shape,"P_commute","CIMD_SDD")

##finding the local SSIM for the polygon map and saving the results to the result_polygon1 and the result_polygon2 respectively
result_polygon1<-ssim_polygon(shape,"PP_SDD","CIMD_SDD",global=FALSE)
result_polygon2<-ssim_polygon(shape,"P_commute","CIMD_SDD",global=FALSE)

##visualization of the local SSIM using the mapview package
library(mapview)

# Define the palfunc function that creates a color palette to be used for the map visualization
 palfunc <- function (n, alpha = 1, begin = 1, end = 0, direction = -1)
  {
   # Create a 11-color palette from the RColorBrewer package called "RdBu"
    colors <- RColorBrewer::brewer.pal(11, "RdBu")
    if (direction < 0) colors <- rev(colors)
    colorRampPalette(colors, alpha = alpha)(n)
 }

mapview(result_polygon1,zcol="SSIM",col.regions = palfunc)
mapview(result_polygon2,zcol="SSIM",col.regions = palfunc)




##example raster maps for singletons sperm whales and groups sperm whales
image_path1 <- system.file("data", "groups2nm.tif", package = "SSIMmap")
image_path2 <- system.file("data", "single2nm.tif", package = "SSIMmap")

img1<-rast(image_path1)
img2<-rast(image_path2)

##finding the global SSIM for the raster maps
ssim_raster(img1,img2)

##finding the local SSIM for the raster maps and saving the results to result_raster
result_raster<-ssim_raster(img1,img2,global=FALSE)

##visualizationg the local SSIM for the raster maps
plot(result_raster)

```
