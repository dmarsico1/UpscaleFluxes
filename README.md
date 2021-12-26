# UpscaleFluxes

Code used in the paper "Upscale momentum fluxes due to the Coriolis force acting on  Poloidal flows."

The folder "Matlab" contains the Matlab files used for figure plotting and post-processing.

The fortran code inverts the elliptic operators using a parallel conjugate gradient method.  The explicit forms of these elliptic operators can be found in the aforementioned paper, and they differ from the standard Laplacian operator due to the presence of cross terms in the y-z plane. 
