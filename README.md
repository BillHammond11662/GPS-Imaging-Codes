# README for GPS Imaging codes

The codes presented here are used to generate a robust interpolation of vertical land motion velocities from GPS station data.    

The method is also a very effective interpolator of randomly distributed point data for more general applications. The idea is similar to Kriging but uses weighted medians instead of least squares inversions to determine the field. Uncertainties and network geometry are accounted for by incorporating a spatial structure function that is based on the data and station locations. The results are therfore insensitive to outlier data and can tolerate heterogenous station distibution.  There is no spatial smoothing applied so edges in the field are preserved as long as the data are capable of constraining them.   

The codes are written in Matlab.  There is brief help at the beginning of each program. e.g., in Matlab use “help code.m” or use a text editor to read the comments in the file header.  There is an example Matlab script provided that invokes the codes (Test.m), this should help people to use the codes to study places other than the example given for California and Nevada.  Try to open Matlab and execute “Test.m”.  Everything was tested on in Matlab version R2022b.  

The scripts should work without any special Matlab toolboxes.  However, it is possible to speed execution time by using the parallel computing toolbox.  There is an option within the example to set a flag (“optpar”) to use this feature.  

The algorithm and results for vertical motions in California and Nevada are presented in:
Hammond, W.C., G. Blewitt, C. Kreemer, 2016, GPS Imaging of vertical land motion in California and Nevada: Implications for Sierra Nevada uplift, Journal of Geophysical Research - Solid Earth, v. 121 (10),  https://onlinelibrary.wiley.com/doi/10.1002/2016JB013458/full,
which is available with open access.  A .pdf of the manuscript is included in these files.

The codes have been used by a few students and other scientists who have provided some feedback on their stability and use, but some bugs may persist...  

Bill Hammond
Nevada Geodetic Laboratory
University of Nevada, Reno
