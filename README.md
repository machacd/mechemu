# Mechanism-based emulator for hydrology

Fortran implementation of a mechanistic emulator applied to hydrological models based on the shallow-water equations.

## Quick intro into the included minimal example

The mechanistic emulator constits of a sum of two parts -- the linear model and
 a Gaussian process conditioned on previous runs of the full model. We deal only
with the implementation and refer to the second article mentioned in the 
[References](https://github.com/machacd/mechemu#references) when it comes to
equations.

### Linear model

The linear model currently implemented (file **core.f90**) consists of *m* coupled 
linear reservoirs, connected in series. One can specifiy, how many of them are
surface reservoirs (they recieve rainfall), which then have their outflow governed
by equation (5).

All the other reservoirs are subsurface linear reservoirs, with the release
coefficient as per equation (10). The output of the emulator is then always read
from the very last reservoir, or the reservoir before that, if the user specifies
*dim_obs=2*.

This linear model is hardcoded into the file **core.f90**, but the user can easily
change it, e.g.\ if she wishes to incorporate a more complex model structure. The 
relevant functions are
```
:w


![equation]().



### Design data

Design data consist of pairs parameter--model output. The user needs to generate them externaly using her model.

## References 

Machac, D.; Reichert, P.; Rieckermann, J. & Albert, C. **Fast emulator of a a slow urban drainage simulator**. *Environmental Modelling & Software*, 2015 (submitted)

Machac, D.; Reichert, P. & Albert, C. **Emulation of dynamic simulators with application to hydrology**. *Journal of Computational Physics*, 2015 (in review)

Albert, C. [**A mechanistic dynamic emulator**](http://arxiv.org/abs/1112.5304). *Nonlinear Analysis: Real World Applications*, 2012, 13, 2747 -- 2754

Reichert, P.; White, G.; Bayarri, M. J. & Pitman, E. B. [**Mechanism-based emulastion of dynamic simulation models: Concept and application in hydrology**](http://dl.acm.org/citation.cfm?id=1923145). *Computational Statistics & Data Analysis*, 2011, 55, 1638--1655


## Acknowledgements

This work is part of the project “Using Commercial Microwave Links and Computer Model Emulation to Reduce Uncertainties in Urban Drainage Simulations” (COMCORDE) funded by
the Swiss National Science Foundation, grants no. CR22I2 135551 and CR22I2 152824.

![Eawag](http://www.fishecology.ch/layout/ealogo-print.gif)

## License
