# Data, code, and Figures for *Soil water depth distribution and species interactions mediate plant functional type dominance in drylands*

#### Renne, Rachel R., Burke, Ingrid C., & Lauenroth, William K.

#### 2025-07-22

#### Repository created by Rachel R. Renne

---

This repository includes the full analysis for our study on how soil
water amount and depth distribution determine shrub-grass dominance and
interactions in sagebrush ecosystems in the western United States. We
visited 51 relatively undisturbed sites and collected data on plant
cover and perennial grass abundance in shrub canopy and interspace
microsites, then used these data to test a new conceptual framework for
predicting plant functional type dominance in drylands. The manuscript
is in preparation and the DOI will be made available after publication.

Raw field data, summarized simulation output results, and metadata can be found in
`01_Data`.

Scripts (`02_Code`) are numbered in order of the analysis:

- `01_Summarize_vegetation&climate.R` 
- `02_Exploring_dominance_clim&soilwater.R`
- `03_Exploring_RII_clim&soilwater.R`
- `04_Vegetation sensitivity_analysis.code.R`
- `05_Shallow_depth_test.R`


Additional scripts containing functions called in the analysis scripts
include:

- `lpi.R`
- `plantcommunity.R`

Figures generated for the publication and supplemental information are
found in `03_Figures`.

## System Information and Package Versions

<details>

<summary>Click here to expand</summary>

```         
R version 4.5.1 (2025-06-13 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Denver
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] daymetr_1.7       viridis_0.6.2     viridisLite_0.4.2 Ternary_2.1.3     tidyr_1.3.0      
[6] terra_1.8-54      dplyr_1.1.2      

loaded via a namespace (and not attached):
 [1] gtable_0.3.3       compiler_4.5.1     promises_1.2.0.1   maps_3.4.3         tidyselect_1.2.1  
 [6] Rcpp_1.0.10        gridExtra_2.3      later_1.3.0        scales_1.2.1       fastmap_1.1.1     
[11] mime_0.12          lattice_0.21-8     ggplot2_3.4.2      R6_2.5.1           generics_0.1.3    
[16] MASS_7.3-59        tibble_3.2.1       maptools_1.1-6     munsell_0.5.0      shiny_1.7.4       
[21] pillar_1.9.0       RColorBrewer_1.1-3 GISTools_0.7-4     rlang_1.1.0        utf8_1.2.3        
[26] sp_2.2-0           httpuv_1.6.9       cli_3.6.1          withr_2.5.0        magrittr_2.0.3    
[31] ncdf4_1.21         digest_0.6.31      grid_4.5.1         xtable_1.8-4       rstudioapi_0.14   
[36] lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         codetools_0.2-19   rgeos_0.6-2       
[41] colorspace_2.1-1   fansi_1.0.4        foreign_0.8-84     purrr_1.0.1        ellipsis_0.3.2    
[46] tools_4.5.1        pkgconfig_2.0.3    htmltools_0.5.5   
```

</details>
