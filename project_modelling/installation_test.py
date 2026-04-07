"""
This file offers a simple test of the installation of HPVsim, to ensure it is working correctly.
"""
import hpvsim as hpv

pars = {}
sim = hpv.Sim(pars=pars)
sim.run()
sim.plot()