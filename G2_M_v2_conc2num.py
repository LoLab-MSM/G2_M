"""A stochastic model based on:

Mathematical modeling and sensitivity analysis of G2/S phase in the cell cycle
involving the p53/Mdm2 oscillation system.  Yoshihiko Tashima, Yu Kisaka, 
Taizo Hanai, Hiroyuki Hamada, Yukihiro Eguchi, Masahiro Okamoto. 
Proc. Int. Fed. Med. Biomed. Eng., 14(2006), pp.195-198
doi:10.1016/j.biosystems.2008.05.016 "Update"


http://www.researchgate.net/publication/226126852_Mathematical_modeling_of_G2M_phase_in_the_cell_cycle_with_involving_the_p53Mdm2_oscillation_system

Implemented by: Corey Hayford
"""
from hayford_model_modified_parameters import *
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
    
#     alias_model_components(model)
#     Parameter("Vol", vol)                #Volume
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

###
#     Parameter("k1", 0.2)        #1/t OK
#     Parameter("k2", 1.0)        #1/t OK
#     Parameter("k3", 1.0)        #1/([]*t) /1/[]
#     Parameter("km3", 1.0)       #1/t OK
#     Parameter("k4", 0.01)        #1/([]*t) /1/[]
#     Parameter("km4", 0.01)        #1/t OK
#     Parameter("k5", 1.0)          #1/([]*t) /1/[]
#     Parameter("km5", 0.01)        #1/t OK
#     Parameter("k6", 1.0)          #1/([]*t) /1/[]
#     Parameter("km6", 0.01)        #1/t OK
#     Parameter("k7", 0.01)        #1/([]*t) /1/[]
#     Parameter("k8", 100.0)       #1/([]*t) /1/[]
#     Parameter("k9", 0.0005)       #[]/t * 1/[]
#     Parameter("k10", 1.0)        #1/([]*t) /1/[]
#     Parameter("km10", 1.0)       #1/([]*t) /1/[]
#     Parameter("k11", 0.1)        #1/([]*t) /1/[]
#     Parameter("km11", 1.0)       #1/t OK
#     Parameter("k12", 0.01)       #1/([]*t) /1/[]
#     Parameter("k13", 1.0)        #1/t OK
#     Parameter("k14", 0.01)       #[]/t * 1/[]
#     Parameter("k15", 0.1)        #1/t OK
#     Parameter("k16", 2.0e-4)     #[]/t *1/[]
#     Parameter("k17", 0.1)        #1/([]*t) /1/[]
#     Parameter("km17", 1.0)        #1/t OK
#     Parameter("k18", 1.0)        #1/t OK
#     Parameter("k19", 1.0)        #1/t OK
#     Parameter("k20", 1.0)        #[]/t * 1/[]
#     Parameter("k21", 0.01)       #1/t OK
#     Parameter("k22", 0.00094)    #[]/t * 1/[]
#     Parameter("k23", 0.02)        #1/t OK 
#     Parameter("k24", 10.0)       #[]/t * 1/[]
#     Parameter("k25", 0.005)      #1/t OK
#     Parameter("k26", 0.004)     #1/[]^2 /1/[]^2 
#     Parameter("k27", 6.0)       #1/([]*t) /1/[]
#     Parameter("k28", 0.0001)    #[]/t * 1/[]
#     Parameter("k30", 0.001)     #1/t OK
#     Parameter("k31", 1.0)        #1/[] /1/[]
#     Parameter("k32", 0.0001)     #1/t OK
#     Parameter("k33", 1e-8)      #1/t OK
#     Parameter("k34", 1.5)       #1/t OK
#     Parameter("k_ex", 1.0)        #1/t OK
#     Parameter("v_in", 1.0e-5)   #[]/t * 1/[]
#     Parameter("k_m", 9.5)       #[] * 1/[]
#     Parameter("n", 9.0)        #unitless OK
#     Parameter("k_damp", 0.02)     #1/([]*t) /1/[]
#     Parameter("k_deg", 0.772)     #1/([]^2 * t) /1/[]^2
#     Parameter("Deg_0", 0.0566)    #1/([]*t) /1/[]
    
    
        
# ***Generate ODEs and Plot***

declare_monomers()
declare_parameters()
declare_initial_conditions()
declare_observables()
declare_functions()
declare_rules()
print model.parameters['k3'].value

set_volume(1.0e-20)
print model.parameters['k3'].value
print model.rules['Chk1_Phos']
print model.rules['Chk1_Phos'].rate_forward.value

# quit()
# print model.expressions['sig_deg'].value


# set_dna_damage(0.005) 
# set_volume(1.0)
# print model.parameters['Vol'].value
# print model.parameters['X1_0'].value
# model.parameters['X1_0'].value = 1.0e-6 * constants.N_A * Vol.value
# print model.parameters['X1_0'].value
# quit()

     
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

#####
set_dna_damage(0.0 * Na_V)
y = run_ssa(model,t,verbose=True)

pl.figure()
for obs in ["OBS_MPF", "OBS_p53", "OBS_Wee1"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc='upper right')
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Protein Dynamics (No DNA Damage)")
pl.savefig("Stochastic G2-M Cell Cycle No DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_aCdc25"], label="OBS_aCdc25")
pl.legend(loc='upper left')
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level (Active Cdc25)")
pl.title("Protein Dynamics (No DNA Damage)")
pl.savefig("Stochastic G2-M Cell Cycle No DNA Damage2.png", format= "png")

#####
set_dna_damage(0.005 * Na_V)
y = run_ssa(model,t,verbose=True)

pl.figure()
for obs in ["OBS_MPF", "OBS_p53", "OBS_Wee1"]:
    pl.plot(t, y[obs], label=obs)
pl.legend(loc='upper right')
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level")
pl.title("Protein Dynamics (DNA Damage = 0.005)")
pl.savefig("Stochastic G2-M Cell Cycle DNA Damage1.png", format= "png")

pl.figure()
pl.plot(t, y["OBS_aCdc25"], label="OBS_aCdc25")
pl.legend(loc='upper left')
pl.xlabel("Time (arbitrary units)")
pl.ylabel("Protein Level (Active Cdc25)")
pl.title("Protein Dynamics (DNA Damage = 0.005)")
pl.savefig("Stochastic G2-M Cell Cycle DNA Damage2.png", format= "png")

pl.show()    

