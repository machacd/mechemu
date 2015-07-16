# Mechanism-based emulator for hydrology

Fortran implementation of a mechanistic emulator applied to hydrological models based on the shallow-water equations.

## Quick intro into setup for own problems

For a more detailed information on the theory of mechanistic emulators, see the
referenced articles.

The mechanistic emulator constits of a sum of two parts -- the linear model and
 a Gaussian process conditioned on previous runs of the full model.

### Linear model

The linear model currently implemented (file **core.f90**) consists of *m* coupled 
linear reservoirs, connected in series. One can specifiy, how many of them are
surface reservoirs (they recieve rainfall), which then have their outflow governed
by equation

![equation](http%3A%2F%2Fwww.sciweavers.org%2Fupload%2FTex2Img_1437057429%2Feqn.png).

All the other reservoirs are subsurface reservoirs, governed by equation

![equation](http%3A%2F%2Fwww.sciweavers.org%2Fupload%2FTex2Img_1437057739%2Feqn.png).



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
