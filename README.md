# Introduction to SSIMmap

The goal of the _SSIMmap_ package  extend [the classical SSIM method](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1284395) for irregular lattice-based maps and raster images.The original SSIM method was developed to compare two image mimicking the human visual system. The _SSIMmap_ package applies this method to two types of maps (polygon and raster). The geographical SSIM method incorporates well-developed [geographically weighted summary statistics](https://www.sciencedirect.com/science/article/pii/S0198971501000096?via%3Dihub) with an adaptive bandwidth kernel function for irregular lattice-based maps. This package includes four key functions: **ssim_bandwidth**, **ssim_constant**, **ssim_polygon**,**ssim_raster**. Users who want to compare two maps and quantify their similarities can utilize this package and visualize the results using other R packages(e.g.,[_tmap_](https://CRAN.R-project.org/package=tmap) and [_ggplot_](https://CRAN.R-project.org/package=ggplot2)).


## Installation: SSIMmap

You can install the development version of _SSIMmap_ from [_GitHub_](https://github.com/Hailyee-Ha/SSIMmap) with the following commands:

```
devtools::install_github("Hailyee-Ha/SSIMmap")
library(SSIMmap)
```

## Data: Toronto, groups2nm.tif, and single2nm.tif

The package has three example data: 1) Toronto, 2) groups2nm (*terra raster object*), and 3) single2nm (*terra raster object*). First, **Toronto**(*an object of sf class*) is an example data, which stands for Toronto, ON. This contains neighborhood deprivation indices (the Canadian Index of Multiple Deprivation and Pampalon index) and the census variable (the percentage of households who commute within the census subdivision of residence) for 2016. The second and third data are *terra raster objects*. These raster objects represent the geographical location from 38 to 41◦N and 0.5–5◦E, centering on the islands of Ibiza, Mallorca and Menorca. The raster objects, labeled as 'groups2nm' and 'single2nm',  contain information on sperm whales. They represent the predicted probability of presence for groups and singletons at a spatial resolution of 2 nautical miles (NM) on a regular grid. For further information about 'groups2nm' and 'single2nm' please refer to the paper ["Novel application of a quantitative spatial comparison tool to species distribution data"](https://www.sciencedirect.com/science/article/pii/S1470160X16302990).

```
library(SSIMmap)
library(sf)
library(terra)
shape<- SSIMmap::Toronto
plot(shape, border=NA)
img1<-terra::rast(system.file("data", "groups2nm.tif", package = "SSIMmap"))
img2<-terra::rast(system.file("data", "single2nm.tif", package = "SSIMmap"))
plot(img1)
plot(img2)
```

## Functions: ssim_bandwidth

This function calculates the bandwidth size for the computation of the Structural Similarity Index (SSIM) on polygon maps. A user can decide whether or not to standardize the variables of maps. If the maps have different ranges of variables or include negative values, it is necessary to use **standardize =TRUE**. This ensures that the conditions for SSIM (symmetry and boundedness (a range from -1 to 1)) are not violated. If the maps share a similar range, there is no need to use **standardize= TRUE**.The function provides two options for selecting the bandwidth size, which will be displayed in the console window:

1) The square root of N.
2) The trade-off between the bias and variance of the two maps.

The first method for the bandwidth selection calculates the square root of the sample size n, which is based on the relative rate of convergence that could minimize the mean integrated squared error ([Abramson,1982](https://projecteuclid.org/journals/annals-of-statistics/volume-10/issue-4/On-Bandwidth-Variation-in-Kernel-Estimates-A-Square-Root-Law/10.1214/aos/1176345986.full)). The second method aims to balance local variances and bias in spatial data (geographically weighted mean), with an intent to minimize the difference between them. The primary goal of this method is to find a bandwidth size that balances the bias-variance trade-off and therefore effectively captures the scale of the underlying spatial structure. More information about the bias-variance trade-off can be found in the paper ["On the measurement of bias in geographically weighted regression models"](https://www.sciencedirect.com/science/article/pii/S2211675320300476). By default for the second method, the function chooses the mid-point of the optimal trade-off between the bias and variance for map1 and map2. A user can choose the upper (the larger one) or the lower (the smaller one) value between the two of the optimal trade-off between the bias and variance.

The function takes as input an object of  *sf* class, which includes columns necessary for SSIM calculation. It outputs two possible selections for the bandwidth size based on the selected method. Additionally, the function generates a plot illustrating the relationship between bias, variance, and bandwidth size with a vertical line of the square root of N result.

### How to execute

```
ssim_bandwidth(shape,"CIMD_SDD","PP_SDD",max_bandwidth=80)
ssim_bandwidth(shape,"CIMD_SDD","P_commute",max_bandwidth=80)

```

## Functions: ssim_constant

This function calculates constants (k1 and k2) for the computation of the SSIM index on polygon maps based on the ranges of values. It considers a range of an object of *sf* class including columns for the SSIM calculation and returns the rescaled constants on the console window.

### How to execute

```
ssim_constant(shape,"CIMD_SDD","PP_SDD")
```

## Functions: ssim_polygon

This function calculates the SSIM index for a given polygon. It takes an object of the *sf* class as input, which includes the necessary columns for SSIM calculation. The function then returns either the global SSIM values (when global = TRUE) or the local SSIM values for each given polygon (when global = FALSE). By default, the bandwidth size is determined by the square root of N. A user has the flexibility to select their own bandwidth size or use the suggestions from the *ssim_bandwidth* function. 

### How to execute

```
ssim_polygon(shape,"CIMD_SDD","PP_SDD") 
ssim_polygon(shape,"CIMD_SDD","P_commute") 
df<-ssim_polygon(shape,"CIMD_SDD","PP_SDD",global = FALSE) 
df_2<-ssim_polygon(shape,"CIMD_SDD","P_commute",global = FALSE) 
plot(df,border=NA)
plot(df_2,border=NA)
```

## Functions: ssim_raster

This function calculates the SSIM index for raster images. It takes as input a image file importing from the [_terra_](https://CRAN.R-project.org/package=terra) and returns either the global SSIM values (global=TRUE) or the SSIM values for each given cell as the local SSIM (global=FALSE). Default of the window size is 3\*3 and a user can use own window size.


### How to execute

```
library(rasterVis)
ssim_raster(img1,img2)
result_raster<-ssim_raster(img1,img2,global = FALSE)
levelplot(result_raster)
```

## Example of visualizing the SSIM for irregular lattice maps using ggplot2

```
library(ggplot2)
library(RColorBrewer)

# Visualize the SSIM result (the CIMD vs. Pampalon) on the map
ggplot() +
  geom_sf(data = df, aes(fill = SSIM), color = NA) +
  scale_fill_gradientn(colors = brewer.pal(8, "PRGn"), limits = c(-1, 1)) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

# Visualize the SSIM result (the CIMD vs. commuting pattern) on the map
ggplot() +
  geom_sf(data = df_2, aes(fill = SSIM), color = NA) +
  scale_fill_gradientn(colors =brewer.pal(8, "PRGn"), limits = c(-1, 1)) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

```

### Reference

* Brunsdon, C., Fotheringham, A.S. and Charlton, M. (2002). Geographically weighted summary statistics---a framework for localised exploratory data analysis. *Computers, Environment and Urban Systems*, 26(6), pp.501-524. https://doi.org/10.1016/S0198-9715(01)00009-6
* Jones, E.L., Rendell, L., Pirotta, E. and Long, J.A. (2016). Novel application of a quantitative spatial comparison tool to species distribution data. *Ecological Indicators*, 70, pp.67-76. https://doi.org/10.1016/j.ecolind.2016.05.051
* Lu, B., Harris, P., Charlton, M. and Brunsdon, C. (2014). The GWmodel R package: further topics for exploring spatial heterogeneity using geographically weighted models. *Geo-spatial Information Science*, 17(2), pp.85-101. https://doi.org/10.1080/10095020.2014.917453
* Wang, Z., Bovik, A.C., Sheikh, H.R. and Simoncelli, E.P. (2004). Image quality assessment: from error visibility to structural similarity. *IEEE transactions on image processing*, 13(4), pp.600-612. https://doi.org/10.1109/TIP.2003.819861
