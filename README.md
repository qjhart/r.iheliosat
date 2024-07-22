# DESCRIPTION

r.iheliosat computes integrated beam (direct), diffuse and ground reflected solar irradiation
raster maps for given day, region and atmospheric conditions.

If users select a time, and r.iheliosat will compute the integrated value at that time. 
Two calls to r.heliosat can then be used to calculate time-slices of irradiance that can be combined 
with cloud cover estimates to calculate non-clear sky estimations.

No shadowing from topography is included in the calculations.

The model computes all three components of global radiation (beam, diffuse and reflected) for the 
clear sky conditions, i.e. not taking into consideration the spatial and temporal variation of clouds.

# NOTES
Solar energy is an important input parameter in different models concerning energy industry, landscape, 
vegetation, evapotranspiration, snowmelt or remote sensing. Solar rays incidence angle maps can be 
effectively used in radiometric and topographic corrections in mountainous and hilly terrain where very 
accurate investigations should be performed.

The clear-sky solar radiation model applied in the r.iheliosat is based on the work undertaken for 
development of European Solar Radiation Atlas (Rigollier 2001). The clear sky model estimates the global 
radiation from the sum of its beam, diffuse and reflected components.

# EXAMPLES

# SEE ALSO
[r.sun](https://grass.osgeo.org/grass77/manuals/r.sun.html)

# REFERENCES
Kasten, F. (1996). The Linke turbidity factor based on improved values of the integral Rayleigh optical thickness. Solar Energy, 56 (3), 239-244.
Rigollier, Ch., Bauer, O., Wald, L. (2000). On the clear sky model of the ESRA - European Solar radiation Atlas - with respect to the Heliosat method. Solar energy, 68, 33-48.

# AUTHORS
Quinn Hart, University of California, Davis
© 2004-2024 Quinn Hart. This program is free software under the MIT License
qjhart@ucdavis.edu
