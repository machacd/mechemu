# Mechanism-based emulator for hydrology

Fortran implementation of a mechanistic emulator applied to hydrological models based on the shallow-water equations.

## Purpose

Purpose of this software package is to approximate slow hydrological models by a fast mechanistic emulator. Mechanism-base emulator is a subclass of Gaussian process emulators, which also uses our knowledge of the original simulator.

## Quick intro

A mechanistic emulator consists of a sum of two parts -- the linear model and
 a Gaussian process conditioned on previous runs of the full model. In this readme, we deal only
with the implementation and refer to the second article mentioned in the 
[References](https://github.com/machacd/mechemu#references) when it comes to
equations.

## Folder structure

Source code of the emulator itself is stored in the folder *core*, compiled executables and libraries are stored in *interface* and the example data files are in *example*.

## Prerequisites

TBD

## Minimal included example tutorial

### Linear model set-up

First, we need to decide on which linear model we want to use. This model is then written, as a function of hyperparameters *h* and parameters *p* in the file *example/config.cfg*. After this is done, cd to *interface* and run

```make
make CONFIG_FILE="../example/config.cfg"
```

which compiles the emulator with the specified linear model. 

### Design data

Furthemore, we need the *design data* -- input parameters and corresponding outputs. The example inputs are stored in the file *design_pars.dat*, where each column represents a certain parameter *p*, which need to be in the same order as used in the config file. The file *design_outputs.dat* has the corresponding outputs.In the provided example, the file has 72 lines -- 64 design data sets generated with *Latin hypercube sampling* and 8 test sets.

### Running the emulator

The emulator is then run from the file *main.py*, which includes the bare minimum in order to produce output plot. In time, this file will include more functions, which are already in the classes files in the *interface folder*.


## References 

Machac, D.; Reichert, P.; Rieckermann, J.; Del Giudice, D. & Albert, C. **Accelerating Bayesian inference in hydrological modelling with a mechanistic emulator**, 2016
[**Source data for the article**](https://dl.dropboxusercontent.com/u/23101549/results.zip)

Machac, D.; Reichert, P.; Rieckermann, J. & Albert, C. **Fast emulator of a a slow urban drainage simulator**. *Environmental Modelling & Software*, 2016, 78, 54 -- 67

Machac, D.; Reichert, P. & Albert, C. **Emulation of dynamic simulators with application to hydrology**. *Journal of Computational Physics*, 2016, 313, 352 -- 366

Albert, C. [**A mechanistic dynamic emulator**](http://arxiv.org/abs/1112.5304). *Nonlinear Analysis: Real World Applications*, 2012, 13, 2747 -- 2754

Reichert, P.; White, G.; Bayarri, M. J. & Pitman, E. B. [**Mechanism-based emulastion of dynamic simulation models: Concept and application in hydrology**](http://dl.acm.org/citation.cfm?id=1923145). *Computational Statistics & Data Analysis*, 2011, 55, 1638--1655


## Acknowledgements

This work is part of the project “Using Commercial Microwave Links and Computer Model Emulation to Reduce Uncertainties in Urban Drainage Simulations” (COMCORDE) funded by
the Swiss National Science Foundation, grants no. CR22I2 135551 and CR22I2 152824.

![Eawag](http://www.fishecology.ch/layout/ealogo-print.gif)

## License

TBD
