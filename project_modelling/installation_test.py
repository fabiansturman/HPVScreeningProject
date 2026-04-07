"""
This file offers a simple test of the installation of HPVsim, to ensure it is working correctly.
"""
import hpvsim as hpv

pars = {
    'n_agents':200_000,
    #'debut':dict(f=dict(dist='normal', par1=1.0, par2=3.1), m=dict(dist='normal', par1=1.0, par2=4.1)),
    'debut':dict(f=dict(dist='exponnorm', par1=0.05, par2=0, par3=0.5), m=dict(dist='exponnorm',par1=0.05, par2=0, par3=0.5))
    #'debut':dict(f=dict(dist='exponnorm', par1=2.0408916163272464, par2=14.86159360812454, par3=1.2446420329339398), m=dict(dist='exponnorm',  par1=2.0408916163272464, par2=14.86159360812454, par3=1.2446420329339398))
    }
sim = hpv.Sim(pars=pars)
sim.run()
sim.plot()