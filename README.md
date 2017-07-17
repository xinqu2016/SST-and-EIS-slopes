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
| annual mean |cloud cover  | cl | %  | nc
|             |skin temperature | ts | K |nc  
|             |Surface temperature | tas  | K | nc
|             |air temperature     | ta   |K |nc
|             |surface pressure    | ps   |Pa  |nc   
| constant    |land cover          | sftlf |  |nc

These variables were downloaded from https://esgf-node.llnl.gov/search/esgf-llnl/.
# Output

Intermediate output: time series of low cloud cover, EIS and SST during the period of 1900-1999

Final output: EIS and SST slope averaged over five regions from 36 CMIP models [figure and nc file, which only contains values for 18 CMIP5 models.

Is a script to draw a figure in the paper included ?: No
