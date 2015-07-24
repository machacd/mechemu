import numpy as np
import common_functions as cf
import run_swwm_for_calib_output as rsfco

class swmm(object):
    def __init__(self):
        self.typ="swmm"
        initrun=rsfco.run_swmm([1]*10,"4","own")
        self.t=initrun.shape[0]

    def run(self,pars):
        self.pars=pars
        self.result=rsfco.run_swmm(self.pars,"4","own")

