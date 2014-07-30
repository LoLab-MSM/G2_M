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

# Model()

# def set_dna_damage(damage):
#     model.parameters['DDS_0'].value = damage

def declare_monomers():
    """Declare the monomers in the Tashima model
    'state' is the activity state between the monomers
    'a' is the active conformation binding state
    'i' is the inactive conformation binding state
    'phos' is the phosphorylation state
    'u' is the unphosphorylated monomer
    'p' is the phosphorylated monomer
    'b' is the binding site between p## and other proteins"""
    
#     Monomer("Signal")       #State can be on or off
#     Monomer("SignalDamp")   #Dampens signal in Deg(t) function
        
    # **Regulatory proteins**
     
#     Monomer("p53")
#     Monomer("p21", ['b'])
#     
#     Monomer("x14_3_3", ['b'])
#     Monomer("ATM_ATR", ['b'])
#     Monomer("Mdm2" ['b'])
#     Monomer("I")
#     
#     Monomer('Cdc25', ['b', 'state','state1', 'phos'], {'state':['i','a'], 'state1':['A','C'], 'phos':['u','p']})
#     Monomer('Wee1', ['phos'], {'phos':['u','p']})
#     Monomer('Chk1', ['phos'], {'phos':['u','p']})
#     
# #     Monomer('MPF' , ['b', 'state'], {'state':['i','a']})
#       #MPF split into Cyclin B and CDK1_nuc
# 
#     # **Cyclins** #MPF Complex
#     
#     Monomer("CycB",[           'c'])
#     
#     # **Cyclin-dependent kinases** """ #MPF Complex
#     
#     Monomer("CDK1_nuc",['phos','b','c'], {'phos':['u','p']})
    

def declare_parameters():
    
# ***Declare Initial Conditions***
    
    Parameter("X1_0", 1.0e-6)         #Chk1p
    Parameter("X1pre_0", 0.9999999)   #Chk1
    Parameter("X2_0", 0.2)            #ATM_ATR
    Parameter("X3_0", 0.0265)         #p53
    Parameter("X4_0", 1.0e-6)         #CycB/CDK1_nuc inactive by Wee1
    Parameter("X5_0", 1.0e-8)         #CycB/CDK1_nuc active by Cdc25
    Parameter("X6_0", 0.0)            #p21
    Parameter("X7_0", 0.0)            #p21/CycB/CDK1_nuc inactive by p21 (CKI)
    Parameter("X8_0", 1.0e-6)         #iCdc25C
    Parameter("X9_0", 2.0e-5)         #iCdc25CPs216
    Parameter("X10_0", 0.03)          #iCdc25CPs216/14_3_3
    Parameter("X11_0", 1.0e-6)        #aCdc25C
    Parameter("X12_0", 0.0)           #aCdc25CPs216
    Parameter("X13_0", 2.0)           #x14_3_3
    Parameter("X14_0", 0.001)         #Wee1
    Parameter("X15_0", 0.0)           #Wee1p
    Parameter("X16_0", 2.35e-4)       #Mdm2
    Parameter("X17_0", 0.0)           #I(Intermediate)
    Parameter("DDS_0")                #DDS
    
# ***Declare Kinetic Parameters***
   
    Parameter("G2_M_k1", 0.2)    
    Parameter("G2_M_k2", 1.0)    
    Parameter("G2_M_k3", 1.0)    
    Parameter("G2_M_km3", 1.0)   
    Parameter("G2_M_k4", 0.01)  
    Parameter("G2_M_km4", 0.01)
    Parameter("G2_M_k5", 1.0)  
    Parameter("G2_M_km5", 0.01)
    Parameter("G2_M_k6", 1.0)  
    Parameter("G2_M_km6", 0.01)    
    Parameter("G2_M_k7", 0.01)    
    Parameter("G2_M_k8", 100.0)   
    Parameter("G2_M_k9", 0.0005)   
    Parameter("G2_M_k10", 1.0)    
    Parameter("G2_M_km10", 1.0)   
    Parameter("G2_M_k11", 0.1)    
    Parameter("G2_M_km11", 1.0)   
    Parameter("G2_M_k12", 0.01)   
    Parameter("G2_M_k13", 1.0)    
    Parameter("G2_M_k14", 0.01)   
    Parameter("G2_M_k15", 0.1)    
    Parameter("G2_M_k16", 2.0e-4) 
    Parameter("G2_M_k17", 0.1)    
    Parameter("G2_M_km17", 1.0)    
    Parameter("G2_M_k18", 1.0)    
    Parameter("G2_M_k19", 1.0)    
    Parameter("G2_M_k20", 1.0)    
    Parameter("G2_M_k21", 0.01)   
    Parameter("G2_M_k22", 0.00094)
    Parameter("G2_M_k23", 0.02)    
    Parameter("G2_M_k24", 10.0)   
    Parameter("G2_M_k25", 0.005)  
    Parameter("G2_M_k26", 0.004) 
    Parameter("G2_M_k27", 6.0)   
    Parameter("G2_M_k28", 0.0001)    
    Parameter("G2_M_k30", 0.001)     
    Parameter("G2_M_k31", 1.0)    
    Parameter("G2_M_k32", 0.0001) 
    Parameter("G2_M_k33", 1e-8)  
    Parameter("G2_M_k34", 1.5)   
    Parameter("k_ex", 1.0)    
    Parameter("v_in", 1.0e-5)   
    Parameter("G2_M_k_m", 9.5)    
#     Parameter("n", 9.0)    
    Parameter("k_damp", 0.02) 
    Parameter("k_deg", 0.772) 
    Parameter("Deg_0", 0.0566)
    
    alias_model_components()
    
# ***Initial Conditions***
def declare_initial_conditions():
    x=1
#     Initial(Chk1(phos='p'), X1_0)                                       #0
#     Initial(Chk1(phos='u'), X1pre_0)                                    #1
#     Initial(ATM_ATR(b=None), X2_0)                                          #2
#     Initial(p53(b=None), X3_0)                                                #3
#     
#     Initial(CycB(c=1) % CDK1_nuc(phos='u',b=None,c=1), X4_0)                #4    preMPF (phos=u by convention)
#     Initial(CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1), X5_0)                #5    MPF (phos=p by convention)
#     
# #     Initial(MPF(state='i', b=None), X4_0)                               #4
# #     Initial(MPF(state='a', b=None), X5_0)                               #5
# 
#     Initial(p21(b=None), X6_0)                                          #6
#     
#     Initial(CycB(c=2) % CDK1_nuc(phos='p',b=1,c=2) % p21(b=1), X7_0)        #7    Inactive MPF (sequestered by p21)
#     
# #     Initial(p21(b=1) % MPF(b=1, state='a'), X7_0)                       #7
# 
#     Initial(Cdc25(b=None, state='i', state1='C', phos= 'u'), X8_0)                  #8
#     Initial(Cdc25(b=None, state='i', state1='C', phos='p'), X9_0)                   #9
#     Initial(Cdc25(b=1, state= 'i', state1='C', phos= 'p') % x14_3_3(b=1), X10_0)    #10
#     Initial(Cdc25(b=None, state= 'a', state1='C', phos= 'u'), X11_0)                #11
#     Initial(Cdc25(b=None, state= 'a', state1='C', phos= 'p'), X12_0)                #12
#     Initial(x14_3_3(b=None), X13_0)                                     #13
#     Initial(Wee1(phos= 'u'), X14_0)                                     #14
#     Initial(Wee1(phos= 'p'), X15_0)                                     #15
#     Initial(Mdm2(b=None), X16_0)                                              #16
#     Initial(I(), X17_0)                                                 #17
#     Initial(Signal(), DDS_0)                                            #18
#     Initial(SignalDamp(), DDS_0)                                        #19
    
# ***Observables***
def declare_observables():
    x=1
#     Observable("OBS_aCdc25C", Cdc25(b=None, state= 'a', state1='C', phos= 'u'))
#     Observable("OBS_MPF", CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1))
#     Observable("OBS_p53", p53(b=None))
#     Observable("OBS_Wee1", Wee1(phos= 'u'))
#     Observable("OBS_Int", I())
#     Observable("OBS_Mdm2", Mdm2(b=None))
#     Observable("signal", Signal())
#     Observable("signal_damp", SignalDamp())
  
# ***Functions***
#    Developed by Leonard Harris - modified pysb in lharris branch to process functions
def declare_functions():
    x=1
    
#     Expression("create_preMPF", sympify("G2_M_k9/(1 + G2_M_k31*OBS_p53)"))
#     Expression("create_Mdm2", sympify("G2_M_k24*OBS_Int^n / (G2_M_k_m^n + OBS_Int^n)"))
#     Expression("create_intermediate", sympify("G2_M_k27*OBS_p53/ (1 + G2_M_k26*OBS_p53*OBS_Mdm2)"))
#     Expression("sig_deg", sympify("Deg_0 - k_deg*(signal-signal_damp)"))
#     Expression("kdamp_DDS0", sympify("k_damp*DDS_0"))
#     
# ***Rules***
def declare_rules():
    ## ** Functions **
    Rule('Signal_Damp', SignalDamp() >> None, kdamp_DDS0)
    Rule('Mdm2_Degrade_p53', p53(b=None) + Mdm2(b=None) >> Mdm2(b=None), sig_deg)
      
    Rule('Create_preMPF', None >> CycB(c=1) % CDK1_nuc(phos='u',b=None,c=1), create_preMPF)
#     Rule('Create_preMPF', None >> MPF(b=None, state='i'), create_preMPF) 
  
    Rule('Create_Mdm2_Hill', None >> Mdm2(b=None), create_Mdm2)
    Rule('Signal_Create_I', Signal() >> Signal() + I(), create_intermediate)
    ## ** End Functions **
    
    Rule('Signal_Create_ATM_ATR', Signal() >> Signal() + ATM_ATR(b=None), G2_M_k1)
    Rule('Signal_Create_p53', Signal() >> Signal() + p53(b=None), G2_M_k34)
    Rule('Signal_Degrade', Signal() >> None, G2_M_k33)
    Rule('Chk1_Dephos', Chk1(phos= 'p') >> Chk1(phos= 'u'), G2_M_km3)
    Rule('Chk1_Phos', Chk1(phos= 'u') + ATM_ATR(b=None) >> Chk1(phos= 'p') + ATM_ATR(b=None), G2_M_k3)
    Rule('iCdc25C_Phos', Chk1(phos= 'p') + Cdc25(b=None, state= 'i', state1='C', phos='u') >> Chk1(phos= 'p') + Cdc25(b=None, state= 'i', state1='C', phos= 'p'), G2_M_k7)
    Rule('aCdc25C_Phos', Chk1(phos= 'p') + Cdc25(b=None, state= 'a', state1='C', phos='u') >> Chk1(phos= 'p') + Cdc25(b=None, state= 'a', state1='C', phos= 'p'), G2_M_k4)
    Rule('ATM_ATR_Degrade', ATM_ATR(b=None) >> None, G2_M_k2)
    Rule('p53_Create', None >> p53(b=None), G2_M_k28)
    Rule('p53_Degrade', p53(b=None) >> None, G2_M_k30)
    Rule('p53_Create_p21', p53(b=None) >> p53(b=None) + p21(b=None), G2_M_k15)
    Rule('p53_Create_x14_3_3', p53(b=None) >> p53(b=None) + x14_3_3(b=None), G2_M_k21)
    
    Rule('Activate_MPF', CycB(c=1) % CDK1_nuc(phos='u',b=None,c=1) + Cdc25(b=None, state='a', state1='C') >> CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Cdc25(b=None, state='a', state1='C'), G2_M_k10) 
    Rule('Inactivate_MPF', CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Wee1(phos='u') >> CycB(c=1) % CDK1_nuc(phos='u',b=None,c=1) + Wee1(phos='u'), G2_M_km10)
    Rule('Degrade_MPF', CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) >> None, G2_M_k12)
    Rule('Complex_MPF_p21', CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + p21(b=None) <> CycB(c=1) % CDK1_nuc(phos='u',b=2,c=1) % p21(b=2), G2_M_k11, G2_M_km11)  
    Rule('Activate_Cdc25C', CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Cdc25(b=None, state= 'i', state1='C', phos= 'u') >> CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Cdc25(b=None, state= 'a', state1='C', phos= 'u'), G2_M_k5)
    Rule('Acitvate_Cdc25CPs216', CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Cdc25(b=None, state= 'i', state1='C', phos= 'p') >> CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Cdc25(b=None, state= 'a', state1='C', phos= 'p'), G2_M_k6)
    
#     Rule('Activate_MPF', MPF(b=None, state='i') + Cdc25(b=None, state='a', state1='C') >> MPF(b=None, state='a') + Cdc25(b=None, state='a', state1='C'), G2_M_k10) #Removed 'phos' from Activate_MPF1/2
#     Rule('Inactivate_MPF', MPF(b=None, state='a') + Wee1(phos='u') >> MPF(b=None, state='i') + Wee1(phos='u'), G2_M_km10)
#     Rule('Degrade_MPF', MPF(b=None, state='a') + MPF(b=None, state='a') >> None, G2_M_k12)
#     Rule('Complex_MPF_p21', MPF(b=None, state='a') + p21(b=None) <> MPF(b=1, state= 'a') % p21(b=1), G2_M_k11, G2_M_km11)  
#     Rule('Activate_Cdc25C', MPF(b=None, state='a') + Cdc25(b=None, state= 'i', state1='C', phos= 'u') >> MPF(b=None, state= 'a') + Cdc25(b=None, state= 'a', state1='C', phos= 'u'), G2_M_k5)
#     Rule('Acitvate_Cdc25CPs216', MPF(b=None, state='a') + Cdc25(b=None, state= 'i', state1='C', phos= 'p') >> MPF(b=None, state= 'a') + Cdc25(b=None, state= 'a', state1='C', phos= 'p'), G2_M_k6)
#     Rule('Wee1_Phos', MPF(b=None, state= 'a') + Wee1(phos= 'u') >> MPF(b=None, state= 'a') + Wee1(phos= 'p'), G2_M_k17)

    Rule('Create_p21', None >> p21(b=None), G2_M_k14)
    Rule('Degrade_p21', p21(b=None) >> None, G2_M_k13)
    Rule('Create_iCdc25C', None >> Cdc25(b=None, state= 'i', state1='C', phos= 'u'), v_in)
    Rule('Deactivate_aCdc25C', Cdc25(b=None, state= 'a', phos= 'u') >> Cdc25(b=None, state= 'i', state1='C', phos= 'u'), G2_M_km5)
    Rule('Complex_iCdc25CPs216_x14_3_3', Cdc25(b=None, state= 'i', phos= 'p') + x14_3_3(b=None) >> Cdc25(b=1, state= 'i', state1='C', phos= 'p') % x14_3_3(b=1), G2_M_k8)
    Rule('Deactivate_Cdc25CPs216', Cdc25(b=None, state= 'a', state1='C', phos= 'p') >> Cdc25(b=None, state= 'i', state1='C', phos= 'p'), G2_M_km6)
    Rule('Degrade_Complex_iCdc25CPs216_x14_3_3', Cdc25(b=1, state= 'i', state1='C', phos= 'p') % x14_3_3(b=1) >> None, k_ex)
    Rule('Degrade_aCdc25C', Cdc25(b=None, state= 'a', state1='C', phos= 'u') >> None, G2_M_k32)
    Rule('aCdc25CPs216_Dephos', Cdc25(b=None, state= 'a', state1='C', phos= 'p') >> Cdc25(b=None, state= 'a', state1='C', phos= 'u'), G2_M_km4)
    Rule('Create_x14_3_3', None >> x14_3_3(b=None), G2_M_k20)
    Rule('Degrade_x14_3_3', x14_3_3(b=None) >> None, G2_M_k19)
    Rule('Create_Wee1', None >> Wee1(phos='u'), G2_M_k16)
    Rule('Wee1_Phos', CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Wee1(phos= 'u') >> CycB(c=1) % CDK1_nuc(phos='p',b=None,c=1) + Wee1(phos= 'p'), G2_M_k17)
    Rule('Wee1_Dephos', Wee1(phos= 'p') >> Wee1(phos= 'u'), G2_M_km17)
    Rule('Degrade_Wee1p', Wee1(phos= 'p') >> None, G2_M_k18)
    Rule('Create_Mdm2', None >> Mdm2(b=None), G2_M_k22)
    Rule('Degrade_Mdm2', Mdm2(b=None) >> None, G2_M_k23)
    Rule('Degrade_Intermediate', I() >> None, G2_M_k25)
    