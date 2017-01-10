# SST-and-EIS-slopes
These two metrics measure the sensitivities of low cloud cover to interannual variations in estimated inversion strength (EIS) 
and sea surface tempeature (SST) in five tropical oceanic regions dominanted by low clouds.

# References 
Qu, X., A. Hall, S. A. Klein, P. M. Caldwell, 2014: On the spread of changes in marine low cloud cover in climate model 
simulations of the 21st century. 42: 2603-2626. DOI 10.1007/s00382-013-1945-z.

Wood, R. and C. Bretherton, 2006: On the relationship between stratiform low cloud cover and lower-tropospheric 
stability. J. Clim. 19:6425-6432.
# Input
------------

| Frequency | Variable |  CMOR lables |  Unit  |  File Format |
|:----------  |:--------------------|:----------------|:---------------|:------------|
| annual mean |Cloud cover  | cl | %  | nc
|             |Skin temperature | ts | K |nc  
|             |Surface temperature | tas  | K | nc
|             |air temperature     | ta   |K |nc
|             |surface pressure    | ps   |Pa  |nc   


# Output

EIS and SST slope averaged over five regions from 36 CMIP3 and CMIP5 models  [figure and nc file, %K-1]

Is a script to draw a figure in the paper included ?: No
