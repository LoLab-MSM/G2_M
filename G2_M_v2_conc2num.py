"""A stochastic model based on:

Mathematical modeling and sensitivity analysis of G2/S phase in the cell cycle
involving the p53/Mdm2 oscillation system.  Yoshihiko Tashima, Yu Kisaka, 
Taizo Hanai, Hiroyuki Hamada, Yukihiro Eguchi, Masahiro Okamoto. 
Proc. Int. Fed. Med. Biomed. Eng., 14(2006), pp.195-198
doi:10.1016/j.biosystems.2008.05.016 "Update"


http://www.researchgate.net/publication/226126852_Mathematical_modeling_of_G2M_phase_in_the_cell_cycle_with_involving_the_p53Mdm2_oscillation_system

Implemented by: Corey Hayford
"""
from G2_M_v2_ssa_params import *
from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify
from scipy import constants


alias_model_components(model)


def set_volume(vol):
    global Na_V
    Na_V = constants.N_A * vol              #1/[]
    
    # ***2nd Order Reactions*** (divide by Avogadro's Number * Volume)
    model.parameters['k3'].value /= Na_V
    model.parameters['k4'].value /= Na_V
    model.parameters['k5'].value /= Na_V
    model.parameters['k6'].value /= Na_V
    model.parameters['k7'].value /= Na_V
    model.parameters['k8'].value /= Na_V
    model.parameters['k10'].value /= Na_V
    model.parameters['km10'].value /= Na_V
    model.parameters['k11'].value /= Na_V
    model.parameters['k12'].value /= Na_V
    model.parameters['k17'].value /= Na_V

    # ***1st Order Reactions*** - LEAVE ALONE!
#     model.parameters['k1'].value
#     model.parameters['k2'].value
#     model.parameters['km3'].value
#     model.parameters['km4'].value
#     model.parameters['km5'].value
#     model.parameters['km6'].value
#     model.parameters['km11'].value
#     model.parameters['k13'].value
#     model.parameters['k15'].value
#     model.parameters['km17'].value
#     model.parameters['k18'].value
#     model.parameters['k19'].value
#     model.parameters['k21'].value
#     model.parameters['k23'].value
#     model.parameters['k25'].value
#     model.parameters['k30'].value
#     model.parameters['k32'].value
#     model.parameters['k33'].value
#     model.parameters['k34'].value
#     model.parameters['k_ex'].value
#     model.parameters['n'].value
    
    # ***0th Order Reactions*** (multiply by Avogadro's Number * Volume)
    model.parameters['k14'].value *= Na_V
    model.parameters['k16'].value *= Na_V
    model.parameters['k20'].value *= Na_V
    model.parameters['k22'].value *= Na_V
    model.parameters['k28'].value *= Na_V
    model.parameters['v_in'].value *= Na_V
    
    # ***In Functions *** (act differently)
    model.parameters['k27'].value /= Na_V
    model.parameters['k31'].value /= Na_V
    model.parameters['k_damp'].value /= Na_V
    model.parameters['Deg_0'].value /= Na_V
    model.parameters['k9'].value *= Na_V
    model.parameters['k24'].value *= Na_V
    model.parameters['k_m'].value *= Na_V
    model.parameters['k26'].value /= (Na_V * Na_V)
    model.parameters['k_deg'].value /= (Na_V * Na_V)
    
    # ***Initial Condition Parameters*** (multiply by Avogadro's Number * Volume)
    model.parameters['X1_0'].value *= Na_V
    model.parameters['X1pre_0'].value *= Na_V
    model.parameters['X2_0'].value *= Na_V
    model.parameters['X3_0'].value *= Na_V
    model.parameters['X4_0'].value *= Na_V
    model.parameters['X5_0'].value *= Na_V
    model.parameters['X6_0'].value *= Na_V
    model.parameters['X7_0'].value *= Na_V
    model.parameters['X8_0'].value *= Na_V
    model.parameters['X9_0'].value *= Na_V
    model.parameters['X10_0'].value *= Na_V
    model.parameters['X11_0'].value *= Na_V
    model.parameters['X12_0'].value *= Na_V
    model.parameters['X13_0'].value *= Na_V
    model.parameters['X14_0'].value *= Na_V
    model.parameters['X15_0'].value *= Na_V
    model.parameters['X16_0'].value *= Na_V
    model.parameters['X17_0'].value *= Na_V
    model.parameters['DDS_0'].value *= Na_V
        
# ***Generate ODEs and Plot***

declare_monomers()
declare_parameters()
declare_initial_conditions()
declare_observables()
declare_functions()
declare_rules()

set_volume(1.0e-20)
     
generate_equations(model, verbose=True)
  
# for monomers in model.monomers:
#     print monomers
# print
#   
# for parameters in model.parameters:
#     print parameters
# print
#   
# for initial_conditions in model.initial_conditions:
#     print initial_conditions
# print
  
# for obs in model.observables:
#     print obs, ":", obs.species, ",", obs.coefficients
#     obs_string = ''
#     for i in range(len(obs.coefficients)):
#         if i > 0: obs_string += " + "
#         obs_string += "__s"+str(obs.species[i])
#         if obs.coefficients[i] > 1:
#              obs_string += "*"+str(obs.coefficients[i])
#     print obs_string
#   
# for rules in model.rules:
#     print rules
# print
#   
# for i in range(len(model.species)):
#     print str(i)+":", model.species[i]
# print
#  
# for i in range(len(model.odes)):
#     print str(i)+":", model.odes[i]
# print
#   
# from pysb.generator.bng import BngGenerator
# print BngGenerator(model).get_content()
     
# quit()

t = linspace(0,4000,4000)

## ** Set No DNA Damage

set_dna_damage(0.0 * Na_V) # Because stochastic, multiply by Na_V to get molecule numbers
y = run_ssa(model,t,verbose=True)

pl.figure()
for obs in ["OBS_MPF", "OBS_p53", "OBS_Wee1"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Protein Dynamics (No DNA Damage)")
pl.savefig("Stochastic G2-M Cell Cycle No DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_aCdc25"], label="OBS_aCdc25")
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level (Active Cdc25)")
pl.title("Protein Dynamics (No DNA Damage)")
pl.savefig("Stochastic G2-M Cell Cycle No DNA Damage2.png", format= "png")

#####
set_dna_damage(0.005 * Na_V) # Because stochastic, multiply by Na_V to get molecule numbers
y = run_ssa(model,t,verbose=True)

pl.figure()
for obs in ["OBS_MPF", "OBS_p53", "OBS_Wee1"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Protein Dynamics (DNA Damage = 0.005)")
pl.savefig("Stochastic G2-M Cell Cycle DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_aCdc25"], label="OBS_aCdc25")
pl.legend(loc=0)
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level (Active Cdc25)")
pl.title("Protein Dynamics (DNA Damage = 0.005)")
pl.savefig("Stochastic G2-M Cell Cycle DNA Damage2.png", format= "png")

pl.show()