# Data, code, and Figures for *Soil water depth distribution and species interactions mediate plant functional type dominance in drylands*

#### Renne, Rachel R., Burke, Ingrid C., & Lauenroth, William K.

#### 2025-07-22

#### Repository created by Rachel R. Renne

---

This repository includes the full analysis for our study on how soil
water amount and depth distribution determines shrub-grass dominance and
interactions in sagebrush ecosystems in the western United States. We
visited 51 relatively undisturbed sites and collected data on plant
cover and perennial grass abundance in shrub canopy and interspace
microsites, then used these data to test a new conceptual framework for
predicting plant functional type dominance in drylands. The manuscript
is in preparation and the DOI will be made available after publication.

Raw field data and summarized simulation output results can be found in
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
[1] viridis_0.6.2     viridisLite_0.4.2 tidyr_1.3.0      

loaded via a namespace (and not attached):
 [1] rgeos_0.6-2        utf8_1.2.3         generics_0.1.3     lattice_0.21-8     digest_0.6.31     
 [6] magrittr_2.0.3     grid_4.5.1         RColorBrewer_1.1-3 fastmap_1.1.1      maps_3.4.3        
[11] gridExtra_2.3      promises_1.2.0.1   httr_1.4.7         purrr_1.0.1        fansi_1.0.4       
[16] scales_1.2.1       codetools_0.2-19   cli_3.6.1          shiny_1.7.4        rlang_1.1.0       
[21] ellipsis_0.3.2     munsell_0.5.0      withr_2.5.0        tools_4.5.1        dplyr_1.1.2       
[26] colorspace_2.1-1   ncdf4_1.21         ggplot2_3.4.2      httpuv_1.6.9       curl_5.2.1        
[31] vctrs_0.6.5        R6_2.5.1           mime_0.12          lifecycle_1.0.4    maptools_1.1-6    
[36] MASS_7.3-59        foreign_0.8-84     pkgconfig_2.0.3    terra_1.8-54       pillar_1.9.0      
[41] later_1.3.0        gtable_0.3.3       glue_1.8.0         Rcpp_1.0.10        GISTools_0.7-4    
[46] tibble_3.2.1       tidyselect_1.2.1   rstudioapi_0.14    xtable_1.8-4       htmltools_0.5.5   
[51] compiler_4.5.1     sp_2.2-0    
```

</details>
