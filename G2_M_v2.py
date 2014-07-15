"""An implementation of the model from:

Mathematical modeling and sensitivity analysis of G2/S phase in the cell cycle
involving the p53/Mdm2 oscillation system.  Yoshihiko Tashima, Yu Kisaka, 
Taizo Hanai, Hiroyuki Hamada, Yukihiro Eguchi, Masahiro Okamoto. 
Proc. Int. Fed. Med. Biomed. Eng., 14(2006), pp.195-198
doi:10.1016/j.biosystems.2008.05.016 "Update"


http://www.researchgate.net/publication/226126852_Mathematical_modeling_of_G2M_phase_in_the_cell_cycle_with_involving_the_p53Mdm2_oscillation_system

Implemented by: Corey Hayford
"""
from pysb import *
from pysb.util import *
from pysb.macros import *
from pysb.bng import *
from pysb.integrate import odesolve
import pylab as pl
from numpy import linspace
from sympy import sympify

Model()

def set_dna_damage(damage):
    model.parameters['DDS_0'].value = damage

def declare_monomers():
    """Declare the monomers in the Tashima model
    'state' is the activity state between the monomers
    'a' is the active conformation binding state
    'i' is the inactive conformation binding state
    'phos' is the phosphorylation state
    'u' is the unphosphorylated monomer
    'p' is the phosphorylated monomer
    'b' is the binding site between p## and other proteins"""
    
    Monomer("Signal")       #State can be on or off
    Monomer("SignalDamp")   #Dampens signal in Deg(t) function
        
    # **Regulatory proteins**
     
    Monomer("p53")
    Monomer("p21", ['b'])
    
    Monomer("x14_3_3", ['b'])
    Monomer("ATR", ['b'])
    Monomer("Mdm2")
    Monomer("I")
    
    Monomer('Cdc25', ['b', 'state', 'phos'], {'state':['i','a'], 'phos':['u','p']})
    Monomer('Wee1', ['phos'], {'phos':['u','p']})
    Monomer('Chk1', ['phos'], {'phos':['u','p']})
    
#     Monomer('MPF' , ['b', 'state'], {'state':['i','a']})
      #MPF split into Cyclin B and CDK1

    # **Cyclins** #MPF Complex
    
    Monomer("CycB",[           'c'])
    
    # **Cyclin-dependent kinases** """ #MPF Complex
    
    Monomer("CDK1",['phos','b','c'], {'phos':['u','p']})
    

def declare_parameters():
    
# ***Declare Initial Conditions***
    
    Parameter("X1_0", 1.0e-6)         #Chk1p
    Parameter("X1pre_0", 0.9999999)   #Chk1
    Parameter("X2_0", 0.2)            #ATR
    Parameter("X3_0", 0.0265)         #p53
    Parameter("X4_0", 1.0e-6)         #CycB/CDK1 inactive by Wee1
    Parameter("X5_0", 1.0e-8)         #CycB/CDK1 active by Cdc25
    Parameter("X6_0", 0.0)            #p21
    Parameter("X7_0", 0.0)            #p21/CycB/CDK1 inactive by p21 (CKI)
    Parameter("X8_0", 1.0e-6)         #iCdc25
    Parameter("X9_0", 2.0e-5)         #iCdc25Ps216
    Parameter("X10_0", 0.03)          #iCdc25Ps216/14-3-3
    Parameter("X11_0", 1.0e-6)        #aCdc25
    Parameter("X12_0", 0.0)           #aCdc25Ps216
    Parameter("X13_0", 2.0)           #x14_3_3
    Parameter("X14_0", 0.001)         #Wee1
    Parameter("X15_0", 0.0)           #Wee1p
    Parameter("X16_0", 2.35e-4)       #Mdm2
    Parameter("X17_0", 0.0)           #I(Intermediate)
    Parameter("DDS_0")                #DDS
    
# ***Declare Kinetic Parameters***
   
    Parameter("k1", 0.2)    
    Parameter("k2", 1.0)    
    Parameter("k3", 1.0)    
    Parameter("km3", 1.0)   
    Parameter("k4", 0.01)  
    Parameter("km4", 0.01)
    Parameter("k5", 1.0)  
    Parameter("km5", 0.01)
    Parameter("k6", 1.0)  
    Parameter("km6", 0.01)    
    Parameter("k7", 0.01)    
    Parameter("k8", 100.0)   
    Parameter("k9", 0.0005)   
    Parameter("k10", 1.0)    
    Parameter("km10", 1.0)   
    Parameter("k11", 0.1)    
    Parameter("km11", 1.0)   
    Parameter("k12", 0.01)   
    Parameter("k13", 1.0)    
    Parameter("k14", 0.01)   
    Parameter("k15", 0.1)    
    Parameter("k16", 2.0e-4) 
    Parameter("k17", 0.1)    
    Parameter("km17", 1.0)    
    Parameter("k18", 1.0)    
    Parameter("k19", 1.0)    
    Parameter("k20", 1.0)    
    Parameter("k21", 0.01)   
    Parameter("k22", 0.00094)
    Parameter("k23", 0.02)    
    Parameter("k24", 10.0)   
    Parameter("k25", 0.005)  
    Parameter("k26", 0.004) 
    Parameter("k27", 6.0)   
    Parameter("k28", 0.0001)    
    Parameter("k30", 0.001)     
    Parameter("k31", 1.0)    
    Parameter("k32", 0.0001) 
    Parameter("k33", 1e-8)  
    Parameter("k34", 1.5)   
    Parameter("k_ex", 1.0)    
    Parameter("v_in", 1.0e-5)   
    Parameter("k_m", 9.5)    
    Parameter("n", 9.0)    
    Parameter("k_damp", 0.02) 
    Parameter("k_deg", 0.772) 
    Parameter("Deg_0", 0.0566)
    
# ***Initial Conditions***
def declare_initial_conditions():
    
    Initial(Chk1(phos='p'), X1_0)                                       #0
    Initial(Chk1(phos='u'), X1pre_0)                                    #1
    Initial(ATR(b=None), X2_0)                                          #2
    Initial(p53(), X3_0)                                                #3
    
    Initial(CycB(c=1) % CDK1(phos='u',b=None,c=1), X4_0)                #4    preMPF (phos=u by convention)
    Initial(CycB(c=1) % CDK1(phos='p',b=None,c=1), X5_0)                #5    MPF (phos=p by convention)
    
#     Initial(MPF(state='i', b=None), X4_0)                               #4
#     Initial(MPF(state='a', b=None), X5_0)                               #5

    Initial(p21(b=None), X6_0)                                          #6
    
    Initial(CycB(c=2) % CDK1(phos='p',b=1,c=2) % p21(b=1), X7_0)        #7    Inactive MPF (sequestered by p21)
    
#     Initial(p21(b=1) % MPF(b=1, state='a'), X7_0)                       #7

    Initial(Cdc25(b=None, state='i', phos= 'u'), X8_0)                  #8
    Initial(Cdc25(b=None, state='i', phos='p'), X9_0)                   #9
    Initial(Cdc25(b=1, state= 'i', phos= 'p') % x14_3_3(b=1), X10_0)    #10
    Initial(Cdc25(b=None, state= 'a', phos= 'u'), X11_0)                #11
    Initial(Cdc25(b=None, state= 'a', phos= 'p'), X12_0)                #12
    Initial(x14_3_3(b=None), X13_0)                                     #13
    Initial(Wee1(phos= 'u'), X14_0)                                     #14
    Initial(Wee1(phos= 'p'), X15_0)                                     #15
    Initial(Mdm2(), X16_0)                                              #16
    Initial(I(), X17_0)                                                 #17
    Initial(Signal(), DDS_0)                                            #18
    Initial(SignalDamp(), DDS_0)                                        #19
    
# ***Observables***
def declare_observables():
    
    Observable("OBS_aCdc25", Cdc25(b=None, state= 'a', phos= 'u'))
    Observable("OBS_MPF", CycB(c=1) % CDK1(phos='p',b=None,c=1))
    Observable("OBS_p53", p53())
    Observable("OBS_Wee1", Wee1(phos= 'u'))
    Observable("OBS_Int", I())
    Observable("OBS_Mdm2", Mdm2())
    Observable("signal", Signal())
    Observable("signal_damp", SignalDamp())
    
# ***Functions***
#    Developed by Leonard Harris - modified pysb in lharris branch to process functions
def declare_functions():
    
    Expression("create_preMPF", sympify("k9/(1 + k31*OBS_p53)"))
    Expression("Hill_Mdm2", sympify("k24*OBS_Int^n / (k_m^n + OBS_Int^n)"))
    Expression("create_intermediate", sympify("k27*OBS_p53/ (1 + k26*OBS_p53*OBS_Mdm2)"))
    Expression("sig_deg", sympify("Deg_0 - k_deg*(signal-signal_damp)"))
    Expression("kdamp_DDS0", sympify("k_damp*DDS_0"))
    
# ***Rules***
def declare_rules():
    ## ** Functions **
    Rule('Signal_Damp', SignalDamp() >> None, kdamp_DDS0)
    Rule('p53_Create_Mdm2', p53() + Mdm2() >> Mdm2(), sig_deg)
    
    Rule('Create_preMPF', None >> CycB(c=1) % CDK1(phos='u',b=None,c=1), create_preMPF)
#     Rule('Create_preMPF', None >> MPF(b=None, state='i'), create_preMPF) 

    Rule('Create_Mdm2_Hill', None >> Mdm2(), Hill_Mdm2)
    Rule('Signal_3', Signal() >> Signal() + I(), create_intermediate)
    ## ** End Functions **
    
    Rule('Signal_1', Signal() >> Signal() + ATR(b=None), k1)
    Rule('Signal_2', Signal() >> Signal() + p53(), k34)
    Rule('Signal_Degrade', Signal() >> None, k33)
    Rule('Chk1_Dephos', Chk1(phos= 'p') >> Chk1(phos= 'u'), km3)
    Rule('Chk1_Phos', Chk1(phos= 'u') + ATR(b=None) >> Chk1(phos= 'p') + ATR(b=None), k3)
    Rule('iCdc25_Phos', Chk1(phos= 'p') + Cdc25(b=None, state= 'i', phos='u') >> Chk1(phos= 'p') + Cdc25(b=None, state= 'i', phos= 'p'), k7)
    Rule('aCdc25_Phos', Chk1(phos= 'p') + Cdc25(b=None, state= 'a', phos='u') >> Chk1(phos= 'p') + Cdc25(b=None, state= 'a', phos= 'p'), k4)
    Rule('ATR_Degrade', ATR(b=None) >> None, k2)
    Rule('p53_Create', None >> p53(), k28)
    Rule('p53_Degrade', p53() >> None, k30)
    Rule('p53_Create_p21', p53() >> p53() + p21(b=None), k15)
    Rule('p53_Create_x14_3_3', p53() >> p53() + x14_3_3(b=None), k21)
    
    Rule('Activate_MPF', CycB(c=1) % CDK1(phos='u',b=None,c=1) + Cdc25(b=None, state='a') >> CycB(c=1) % CDK1(phos='p',b=None,c=1) + Cdc25(b=None, state='a'), k10) 
    Rule('Inactivate_MPF', CycB(c=1) % CDK1(phos='p',b=None,c=1) + Wee1(phos='u') >> CycB(c=1) % CDK1(phos='u',b=None,c=1) + Wee1(phos='u'), km10)
    Rule('Degrade_MPF', CycB(c=1) % CDK1(phos='p',b=None,c=1) + CycB(c=1) % CDK1(phos='p',b=None,c=1) >> None, k12)
    Rule('Complex_MPF_p21', CycB(c=1) % CDK1(phos='p',b=None,c=1) + p21(b=None) <> CycB(c=2) % CDK1(phos='u',b=1,c=2) % p21(b=1), k11, km11)  
    Rule('Activate_Cdc25', CycB(c=1) % CDK1(phos='p',b=None,c=1) + Cdc25(b=None, state= 'i', phos= 'u') >> CycB(c=1) % CDK1(phos='p',b=None,c=1) + Cdc25(b=None, state= 'a', phos= 'u'), k5)
    Rule('Acitvate_Cdc25Ps216', CycB(c=1) % CDK1(phos='p',b=None,c=1) + Cdc25(b=None, state= 'i', phos= 'p') >> CycB(c=1) % CDK1(phos='p',b=None,c=1) + Cdc25(b=None, state= 'a', phos= 'p'), k6)
    Rule('Wee1_Phos', CycB(c=1) % CDK1(phos='p',b=None,c=1) + Wee1(phos= 'u') >> CycB(c=1) % CDK1(phos='p',b=None,c=1) + Wee1(phos= 'p'), k17)
    
#     Rule('Activate_MPF', MPF(b=None, state='i') + Cdc25(b=None, state='a') >> MPF(b=None, state='a') + Cdc25(b=None, state='a'), k10) #Removed 'phos' from Activate_MPF1/2
#     Rule('Inactivate_MPF', MPF(b=None, state='a') + Wee1(phos='u') >> MPF(b=None, state='i') + Wee1(phos='u'), km10)
#     Rule('Degrade_MPF', MPF(b=None, state='a') + MPF(b=None, state='a') >> None, k12)
#     Rule('Complex_MPF_p21', MPF(b=None, state='a') + p21(b=None) <> MPF(b=1, state= 'a') % p21(b=1), k11, km11)  
#     Rule('Activate_Cdc25', MPF(b=None, state='a') + Cdc25(b=None, state= 'i', phos= 'u') >> MPF(b=None, state= 'a') + Cdc25(b=None, state= 'a', phos= 'u'), k5)
#     Rule('Acitvate_Cdc25Ps216', MPF(b=None, state='a') + Cdc25(b=None, state= 'i', phos= 'p') >> MPF(b=None, state= 'a') + Cdc25(b=None, state= 'a', phos= 'p'), k6)
#     Rule('Wee1_Phos', MPF(b=None, state= 'a') + Wee1(phos= 'u') >> MPF(b=None, state= 'a') + Wee1(phos= 'p'), k17)

    Rule('Create_p21', None >> p21(b=None), k14)
    Rule('Degrade_p21', p21(b=None) >> None, k13)
    Rule('Create_iCdc25', None >> Cdc25(b=None, state= 'i', phos= 'u'), v_in)
    Rule('Inactivate_aCdc25', Cdc25(b=None, state= 'a', phos= 'u') >> Cdc25(b=None, state= 'i', phos= 'u'), km5)
    Rule('Complex_iCdc25Ps216_x14_3_3', Cdc25(b=None, state= 'i', phos= 'p') + x14_3_3(b=None) >> Cdc25(b=1, state= 'i', phos= 'p') % x14_3_3(b=1), k8)
    Rule('Inactive_Cdc25Ps216', Cdc25(b=None, state= 'a', phos= 'p') >> Cdc25(b=None, state= 'i', phos= 'p'), km6)
    Rule('Degrade_Complex_iCdc25Ps216_x14_3_3', Cdc25(b=1, state= 'i', phos= 'p') % x14_3_3(b=1) >> None, k_ex)
    Rule('Degrade_aCdc25', Cdc25(b=None, state= 'a', phos= 'u') >> None, k32)
    Rule('aCdc25Ps216_Dephos', Cdc25(b=None, state= 'a', phos= 'p') >> Cdc25(b=None, state= 'a', phos= 'u'), km4)
    Rule('Create_x14_3_3', None >> x14_3_3(b=None), k20)
    Rule('Degrade_x14_3_3', x14_3_3(b=None) >> None, k19)
    Rule('Create_Wee1', None >> Wee1(phos='u'), k16)
    Rule('Wee1_Dephos', Wee1(phos= 'p') >> Wee1(phos= 'u'), km17)
    Rule('Degrade_Wee1p', Wee1(phos= 'p') >> None, k18)
    Rule('Create_Mdm2', None >> Mdm2(), k22)
    Rule('Degrade_Mdm2', Mdm2() >> None, k23)
    Rule('Degrade_Intermediate', I() >> None, k25)