
# SSIMmap

<!-- badges: start -->
<!-- badges: end -->

The goal of the _SSIMmap_ package is to offer a comprehensive set of tools for working with Structural Similarity Index (SSIM) for map comparisons.The original SSIM develop to compare two image mimicking the human visual system, which is presented in the [paper](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1284395). _SSIMmap_ package applied this method to maps (polygon and raster). The package includes four key functions for calculating the similarities between maps: **ssim_bandwidth**, **ssim_constant**, **ssim_polygon**,**ssim_raster**. Users who want to compare two maps and quantify their similarities can utilize this package and visualize the results using other tools(e.g.,[_tmap_](https://cran.r-project.org/web/packages/tmap/index.html) and [_ggplot_](https://cran.r-project.org/web/packages/ggplot2/index.html)).  

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

1) **ssim_bandwidth**: This function calculates the bandwidth size for the computation of the SSIM on polygon maps. It offers two options for selecting the bandwidth size from two methods:1.square root of N and 2. the best trade-off between bias and variance. The function takes as input a shape file in *sf* class including columns for the SSIM calculation and returns the two options for selecting bandwidth sizes on the basis of two methods. It also provides a plot illustrating the bias and variance against the size of bandwidth. 

2) **ssim_constant**: This function calculates constants(k1 and k2) for the computation of the SSIM on polygon maps. It takes as input a shape filein *sf* class including columns for the SSIM calculation and returns the constants on the console window.

3) **ssim_polygon**: This function calculates the SSIM index for a given polygon. It takes as input a shape file in *sf* class including columns for the SSIM calculation and returns either the global SSIM values  (global=TRUE) or the SSIM values for each given polygon as the local SSIM  (global=FALSE).

4) **ssim_raster**: This function calculates the SSIM index for raster images. It takes as input a image file importing from the [_terra_](https://cran.r-project.org/web/packages/terra/index.html) and returns either the global SSIM values (global=TRUE) or the SSIM values for each given cell as the local SSIM (global false)

## Example
Here is a basic example which shows you how to use the package:

``` r
library(SSIMmap)

## basic example code
##example polygon map of Toronto including columns(1. Pampalon Index("PP_SDD"), 2. CIMD Index("CIMD_SDD"), and 3. the percent of the housholds who commute within the same census subdivision("P_commute"))
shape<-SSIMmap::polygon 

##finding the options for selecting the bandwidth sizes from two methods
ssim_bandwidth(shape, "PP_SDD","CIMD_SDD", max_bandwidth=500)

##finding the constants(k1 and k2) for the SSIM calculation
ssim_constant(shape,"PP_SDD","CIMD_SDD")

##finding the global SSIM for the polygon map
ssim_polygon(shape,"PP_SDD","CIMD_SDD")

##finding the local SSIM for the polygon map and saving the results to result_polygon
result_polygon<-ssim_polygon(shape,"PP_SDD","CIMD_SDD",global=FALSE)

##example raster maps for singletons sperm whales and groups sperm whales
img1<-SSIMmap::single2nm
img2<-SSIMmap::groups2nm

##finding the global SSIM for the raster maps
ssim_raster(img1,img2)

##finding the local SSIM for the raster maps and saving the results to result_raster
result_raster<-ssim_raster(img1,img2,global=FALSE)
```

