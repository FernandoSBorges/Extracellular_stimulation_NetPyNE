
"""
netParams.py

High-level specifications for S1-thalamus network model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle, json
import os
import numpy as np
import pandas as pd

netParams = specs.NetParams()   # object of class NetParams to store the network parameters


try:
    from __main__ import cfg  # import SimConfig object with params from parent module
except:
    from cfg import cfg

#------------------------------------------------------------------------------
#
# NETWORK PARAMETERS
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# General network parameters
#------------------------------------------------------------------------------
netParams.scale = cfg.scale # Scale factor for number of cells
netParams.sizeX = cfg.sizeX # x-dimension (horizontal length) size in um
netParams.sizeY = cfg.sizeY # y-dimension (vertical height or cortical depth) size in um
netParams.sizeZ = cfg.sizeZ # z-dimension (horizontal depth) size in um
netParams.shape = 'cylinder' # cylindrical (column-like) volume
   
# Layer	height (um)	height (norma)	from	to
# L1	165		    0.079		    0.000	0.079
# L2	149		    0.072		    0.079	0.151
# L3	353		    0.170		    0.151	0.320
# L4	190		    0.091		    0.320	0.412
# L5	525		    0.252		    0.412	0.664
# L6	700		    0.336		    0.664	1.000
# L23	502		    0.241		    0.079	0.320
# All	2082	    1.000	


cellModels = ['HH_full']
Epops = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']

Ipops = []
for popName in cfg.S1pops:
    if popName not in Epops:
        Ipops.append(popName)

layer = {'1':[0.0, 0.079], '2': [0.079,0.151], '3': [0.151,0.320], '23': [0.079,0.320], '4':[0.320,0.412], '5': [0.412,0.664], '6': [0.664,1.0], 
'longS1': [2.2,2.3], 'longS2': [2.3,2.4]}  # normalized layer boundaries

#------------------------------------------------------------------------------
# General connectivity parameters
#------------------------------------------------------------------------------
netParams.defaultThreshold = -10.0 # spike threshold, 10 mV is NetCon default, lower it for all cells
netParams.defaultDelay = 0.1 # default conn delay (ms)
netParams.propVelocity = 300.0 #  300 μm/ms (Stuart et al., 1997)
netParams.scaleConnWeightNetStims = 0.001  # weight conversion factor (from nS to uS)

#------------------------------------------------------------------------------
# Population parameters
#------------------------------------------------------------------------------
## S1

for cellName in cfg.S1cells:
	layernumber = cellName[1:2]
	if layernumber == '2':
		netParams.popParams[cellName] = {'cellType': cellName, 'cellModel': 'HH_full', 'ynormRange': layer['23'], 
                                        'numCells': int(np.ceil(cfg.scaleDensity*cfg.cellNumber[cellName])), 'diversity': True}
	else:
		netParams.popParams[cellName] = {'cellType': cellName, 'cellModel': 'HH_full', 'ynormRange': layer[layernumber], 
                                        'numCells': int(np.ceil(cfg.scaleDensity*cfg.cellNumber[cellName])), 'diversity': True}

#------------------------------------------------------------------------------
# Cell parameters  # L1 70  L23 215  L4 230 L5 260  L6 260  = 1035
#------------------------------------------------------------------------------
## S1 cell property rules

smaller_number_of_axon_sections = {'L1_DAC_bNA': 2, 'L1_DAC_cNA': 1, 'L1_DLAC_cNA': 1, 'L1_HAC_bNA': 3, 'L1_HAC_cIR': 5, 'L1_HAC_cNA': 4, 'L1_NGC_DA_bNA': 3, 'L1_NGC_DA_cAC': 3, 'L1_NGC_DA_cNA': 5, 'L1_NGC_DA_cST': 3, 'L1_NGC_SA_cNA': 4, 'L1_SLAC_bNA': 4, 'L1_SLAC_cAC': 1, 'L1_SLAC_cNA': 1, 'L23_BP_bAC': 1, 'L23_BP_bIR': 1, 'L23_BP_bNA': 2, 'L23_BP_cAC': 3, 'L23_BP_cNA': 2, 'L23_BP_dST': 1, 'L23_BTC_bAC': 5, 'L23_BTC_bIR': 5, 'L23_BTC_bNA': 5, 'L23_BTC_cAC': 5, 'L23_BTC_cNA': 5, 'L23_ChC_cAC': 4, 'L23_ChC_cNA': 3, 'L23_ChC_dNA': 3, 'L23_DBC_bAC': 1, 'L23_DBC_bIR': 1, 'L23_DBC_bNA': 5, 'L23_DBC_cAC': 5, 'L23_LBC_bAC': 3, 'L23_LBC_bNA': 3, 'L23_LBC_cAC': 3, 'L23_LBC_cNA': 4, 'L23_LBC_cST': 3, 'L23_LBC_dNA': 3, 'L23_MC_bAC': 4, 'L23_MC_bNA': 4, 'L23_MC_cAC': 1, 'L23_MC_cNA': 2, 'L23_MC_dNA': 2, 'L23_NBC_bAC': 5, 'L23_NBC_bNA': 1, 'L23_NBC_cAC': 4, 'L23_NBC_cIR': 5, 'L23_NBC_cNA': 1, 'L23_NBC_dNA': 1, 'L23_NGC_bNA': 2, 'L23_NGC_cAC': 2, 'L23_NGC_cNA': 2, 'L23_NGC_cST': 2, 'L23_PC_cAD': 5, 'L23_SBC_bNA': 3, 'L23_SBC_cAC': 3, 'L23_SBC_dNA': 3, 'L4_BP_bAC': 1, 'L4_BP_bIR': 1, 'L4_BP_bNA': 1, 'L4_BP_cAC': 2, 'L4_BP_cNA': 1, 'L4_BP_dST': 1, 'L4_BTC_bAC': 2, 'L4_BTC_bIR': 2, 'L4_BTC_bST': 2, 'L4_BTC_cAC': 4, 'L4_BTC_cNA': 2, 'L4_BTC_dNA': 2, 'L4_ChC_cAC': 1, 'L4_ChC_cNA': 1, 'L4_ChC_dNA': 1, 'L4_DBC_bAC': 1, 'L4_DBC_bIR': 5, 'L4_DBC_bNA': 1, 'L4_DBC_bST': 1, 'L4_DBC_cAC': 1, 'L4_DBC_cIR': 1, 'L4_DBC_cNA': 1, 'L4_LBC_cAC': 1, 'L4_LBC_cNA': 2, 'L4_LBC_cST': 1, 'L4_LBC_dNA': 4, 'L4_LBC_dST': 1, 'L4_MC_bAC': 3, 'L4_MC_bNA': 1, 'L4_MC_cAC': 1, 'L4_MC_cNA': 4, 'L4_MC_dNA': 4, 'L4_NBC_cAC': 1, 'L4_NBC_cIR': 4, 'L4_NBC_cNA': 1, 'L4_NBC_dNA': 2, 'L4_NGC_bNA': 1, 'L4_NGC_cAC': 1, 'L4_NGC_cNA': 2, 'L4_NGC_cST': 1, 'L4_PC_cAD': 4, 'L4_SBC_bNA': 5, 'L4_SBC_cAC': 2, 'L4_SBC_dNA': 2, 'L4_SP_cAD': 2, 'L4_SS_cAD': 3, 'L5_BP_bAC': 1, 'L5_BP_bIR': 1, 'L5_BP_bNA': 1, 'L5_BP_cAC': 3, 'L5_BP_cNA': 5, 'L5_BP_dST': 3, 'L5_BTC_bAC': 1, 'L5_BTC_cAC': 1, 'L5_BTC_cNA': 1, 'L5_ChC_cAC': 1, 'L5_ChC_cNA': 1, 'L5_ChC_dNA': 1, 'L5_DBC_bAC': 2, 'L5_DBC_bIR': 2, 'L5_DBC_bNA': 2, 'L5_DBC_bST': 2, 'L5_DBC_cAC': 2, 'L5_DBC_cIR': 1, 'L5_DBC_cNA': 2, 'L5_LBC_bAC': 4, 'L5_LBC_cAC': 4, 'L5_LBC_cIR': 2, 'L5_LBC_cNA': 5, 'L5_LBC_cST': 1, 'L5_LBC_dNA': 2, 'L5_LBC_dST': 5, 'L5_MC_bAC': 4, 'L5_MC_bIR': 3, 'L5_MC_bST': 5, 'L5_MC_cAC': 5, 'L5_MC_cNA': 5, 'L5_MC_cST': 4, 'L5_MC_dNA': 4, 'L5_NBC_bAC': 1, 'L5_NBC_bIR': 1, 'L5_NBC_bST': 1, 'L5_NBC_cAC': 1, 'L5_NBC_cIR': 1, 'L5_NBC_cNA': 1, 'L5_NBC_cST': 1, 'L5_NBC_dST': 1, 'L5_NGC_bNA': 1, 'L5_NGC_cAC': 1, 'L5_NGC_cNA': 3, 'L5_NGC_cST': 1, 'L5_SBC_bNA': 5, 'L5_SBC_cAC': 1, 'L5_SBC_dNA': 3, 'L5_STPC_cAD': 4, 'L5_TTPC1_cAD': 2, 'L5_TTPC2_cAD': 2, 'L5_UTPC_cAD': 3, 'L6_BPC_cAD': 5, 'L6_BP_bAC': 1, 'L6_BP_bIR': 5, 'L6_BP_bNA': 4, 'L6_BP_cAC': 2, 'L6_BP_cNA': 5, 'L6_BP_dST': 4, 'L6_BTC_bAC': 5, 'L6_BTC_cAC': 3, 'L6_BTC_cNA': 2, 'L6_ChC_cAC': 5, 'L6_ChC_cNA': 4, 'L6_ChC_dNA': 4, 'L6_DBC_bAC': 3, 'L6_DBC_bIR': 3, 'L6_DBC_bNA': 3, 'L6_DBC_bST': 3, 'L6_DBC_cAC': 3, 'L6_DBC_cIR': 2, 'L6_DBC_cNA': 3, 'L6_IPC_cAD': 1, 'L6_LBC_bAC': 5, 'L6_LBC_bIR': 4, 'L6_LBC_bNA': 4, 'L6_LBC_bST': 3, 'L6_LBC_cNA': 3, 'L6_LBC_cST': 3, 'L6_MC_bAC': 1, 'L6_MC_bIR': 2, 'L6_MC_bNA': 2, 'L6_MC_bST': 2, 'L6_MC_cAC': 2, 'L6_MC_cIR': 5, 'L6_MC_cNA': 2, 'L6_NBC_bAC': 2, 'L6_NBC_bIR': 3, 'L6_NBC_bST': 5, 'L6_NBC_cAC': 2, 'L6_NBC_cIR': 2, 'L6_NBC_cNA': 5, 'L6_NBC_cST': 2, 'L6_NBC_dST': 2, 'L6_NGC_bNA': 2, 'L6_NGC_cAC': 2, 'L6_NGC_cNA': 4, 'L6_NGC_cST': 1, 'L6_SBC_bNA': 3, 'L6_SBC_cAC': 5, 'L6_SBC_dNA': 5, 'L6_TPC_L1_cAD': 2, 'L6_TPC_L4_cAD': 3, 'L6_UTPC_cAD': 2}

for cellName in cfg.S1cells:
    
    if int(np.ceil(cfg.scaleDensity*cfg.cellNumber[cellName])) < 5:
        morphoNumbers = int(np.ceil(cfg.scaleDensity*cfg.cellNumber[cellName]))
    else:
        morphoNumbers = 5
        
    cellFraction = 1.0/morphoNumbers
    
    for morphoNumber in range(morphoNumbers):
        cellMe = cfg.cellLabel[cellName] + '_' + str(morphoNumber+1)
        
        # netParams.loadCellParamsRule(label = cellMe, fileName = 'cells/' + cellMe + '_cellParams.json')  
        netParams.loadCellParamsRule(label = cellMe, fileName = 'cells/' + cfg.cellLabel[cellName] + '_' + str(smaller_number_of_axon_sections[cellName]) + '_cellParams.json')  

        netParams.cellParams[cellMe]['diversityFraction'] = cellFraction        
        netParams.cellParams[cellMe]['secLists']['spiny'] = [sec for sec in netParams.cellParams[cellMe]['secLists']['all'] if sec not in netParams.cellParams[cellMe]['secLists']['axonal']]
        netParams.cellParams[cellMe]['secLists']['spinyEE'] = [sec for sec in netParams.cellParams[cellMe]['secLists']['spiny'] if sec not in netParams.cellParams[cellMe]['secLists']['somatic']]
        netParams.cellParams[cellMe]['conds']['cellType'] = cellName
        
        #-----------------------------------------------------------------------------------#
        if cfg.reducedtest:
            cellRule = {'conds': {'cellType': cellName}, 'diversityFraction': cellFraction, 'secs': {}}  # cell rule dict
            cellRule['conds'] = netParams.cellParams[cellMe]['conds']    
            cellRule['secs'] = {}
            cellRule['secs']['soma_0'] = netParams.cellParams[cellMe]['secs']['soma_0']
            cellRule['secLists'] = {}
            cellRule['secLists']['spiny'] = ['soma_0']
            cellRule['secLists']['spinyEE'] = ['soma_0']
            cellRule['secLists']['all'] = ['soma_0']
            cellRule['secLists']['basal'] = ['soma_0']   
            cellRule['secLists']['apical'] = ['soma_0']    
            netParams.cellParams[cellMe] = cellRule   # add dict to list of cell params   
        #-----------------------------------------------------------------------------------#
   
# #------------------------------------------------------------------------------
# #  extracellular mechs
# #------------------------------------------------------------------------------
for celltyp in netParams.cellParams.keys():
    label = []
    for secname in netParams.cellParams[celltyp]['secs'].keys():
        netParams.cellParams[celltyp]['secs'][secname]['mechs']['extracellular'] = {}

# ------------------------------------------------------------------------------
# load data from S1 conn pre-processing file 
#------------------------------------------------------------------------------
with open('conn/conn.pkl', 'rb') as fileObj: connData = pickle.load(fileObj)

lmat = connData['lmat']
a0mat = connData['a0mat']
d0 = connData['d0']

a0e = connData['a0mat_exp']
l0e = connData['lmat_exp']
d0e = connData['d0_exp']

a0g = connData['a0mat_gauss']
x0g = connData['x0_gauss']
l0g = connData['lmat_gauss']
d0g = connData['d0_gauss']

dfinal = connData['dfinal']
pmat = {}
pmat[12.5] = connData['pmat12um']
pmat[25] = connData['pmat25um']
pmat[50] = connData['pmat50um']
pmat[75] = connData['pmat75um']
pmat[100] = connData['pmat100um']
pmat[125] = connData['pmat125um']
pmat[150] = connData['pmat150um']
pmat[175] = connData['pmat175um']
pmat[200] = connData['pmat200um'] #max value for d0=200

synperconnNumber = connData['synperconnNumber']
connNumber = connData['connNumber']
decay = connData['decay']
gsyn = connData['gsyn']
use = connData['use']

ConnTypesNumber = connData['ConnTypesNumber'] 
ConnTypes = connData['ConnTypes']   

connIEtype = connData['connIEtype']  
connEItype = connData['connEItype']
parameters_syn = connData['parameters_syn']

physColumnNames = []
syntypes = []
for name,syntype in parameters_syn:    
    if name not in physColumnNames:
        physColumnNames.append(name) 
    if syntype not in syntypes:
        syntypes.append(syntype)
        
dfS6 = pd.DataFrame(index=syntypes, columns=physColumnNames)
for syntype in syntypes:
    for name in physColumnNames:    
        dfS6.loc[syntype][name] = parameters_syn[name,syntype]

#------------------------------------------------------------------------------
# Synaptic mechanism parameters
#------------------------------------------------------------------------------
#  mods from S1 BBP - deterministic version
for syntype in syntypes:
    if syntype > 50:  # Exc
        
        netParams.synMechParams['S1_EE_STP_Det_' + str(syntype)] = {'mod': 'DetAMPANMDA',
                                         'Use': dfS6['use'][syntype]*cfg.use_frac['EE'], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype] 
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_AMPA': 1.74, # ± 0.18 ms
                                         'tau_r_AMPA': 0.2, 
                                         'tau_r_NMDA': 0.29,
                                         'tau_d_NMDA': 43,   
                                         'NMDA_ratio': 0.8, # ± 0.1 for EE -- experimentally measured for some path?
                                         'mg':1.0, #    0.5mM where exceptionally specified?                                                                
                                            }
        netParams.synMechParams['S1_EIproximal_STP_Det_' + str(syntype)] = {'mod': 'DetAMPANMDA',
                                         'Use': dfS6['use'][syntype]*cfg.use_frac['EIproximal'], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype] 
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_AMPA': 1.74, # ± 0.18 ms
                                         'tau_r_AMPA': 0.2,
                                         'tau_r_NMDA': 0.29,
                                         'tau_d_NMDA': 43,   
                                         'NMDA_ratio': 0.4, # ± 0.1  for EI -- experimentally measured for some path?
                                         'mg':1.0, #    0.5mM where exceptionally specified?                                                                
                                            }
        netParams.synMechParams['S1_EIdistal_STP_Det_' + str(syntype)] = {'mod': 'DetAMPANMDA',
                                         'Use': dfS6['use'][syntype]*cfg.use_frac['EIdistal'], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype] 
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_AMPA': 1.74, # ± 0.18 ms
                                         'tau_r_AMPA': 0.2,
                                         'tau_r_NMDA': 0.29,
                                         'tau_d_NMDA': 43,   
                                         'NMDA_ratio': 0.4, # ± 0.1  for EI -- experimentally measured for some path?
                                         'mg':1.0, #    0.5mM where exceptionally specified?                                                                
                                            }
    else: # Inh
        
        netParams.synMechParams['S1_II_STP_Det_' + str(syntype)] = {'mod': 'DetGABAAB',
                                         'Use': dfS6['use'][syntype]*cfg.use_frac['Inh'], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype]  
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_GABAA': dfS6['decay'][syntype], # ± dfS6['decayStd'][syntype]
                                         'tau_r_GABAA': 0.2,   #rng.lognormal(0.2, 0.1) in synapses.hoc  
                                         'tau_d_GABAB': 260.9,
                                         'tau_r_GABAB': 3.5,
#                                          'GABAB_ratio': 1.0,  #=0(1):The ratio of GABAB to GABAA  ?          
                                            }
        
        netParams.synMechParams['S1_IE_STP_Det_' + str(syntype)] = {'mod': 'DetGABAAB',
                                         'Use': dfS6['use'][syntype]*cfg.use_frac['Inh'], # ± dfS6['useStd'][syntype]
                                         'Dep': dfS6['dep'][syntype], # ± dfS6['depStd'][syntype]  
                                         'Fac': dfS6['fac'][syntype], # ± dfS6['facStd'][syntype]
                                         'tau_d_GABAA': dfS6['decay'][syntype], # ± dfS6['decayStd'][syntype]
                                         'tau_r_GABAA': 0.2,   #rng.lognormal(0.2, 0.1) in synapses.hoc  
                                         'tau_d_GABAB': 260.9,
                                         'tau_r_GABAB': 3.5,
#                                          'GABAB_ratio': 1.0,  #=0(1):The ratio of GABAB to GABAA   ?       
                                            }

# Th NEW
#E2 -> syn 134
netParams.synMechParams['TC:S1'] = {'mod': 'DetAMPANMDA',
                                          'Dep': 227.0,
                                          'Fac': 13.0,
                                          'Use': 0.72,
                                          'tau_r_AMPA': 0.2,
                                          'tau_d_AMPA': 1.74,
                                          'NMDA_ratio': 0.4,
                                          'tau_r_NMDA': 0.29,
                                          'tau_d_NMDA': 43.0}


# Spont and BG
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 1.74, 'e': 0}
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 0.29, 'tau2NMDA': 43, 'e': 0}
netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.2, 'tau2': 8.3, 'e': -80}
netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93} 
ESynMech = ['AMPA', 'NMDA']
ISynMech = ['GABAA', 'GABAB']

# Th
netParams.synMechParams['NMDA_Th']             = {'mod': 'MyExp2SynNMDABB',    'tau1NMDA': 15, 'tau2NMDA': 150,                'e': 0}
netParams.synMechParams['AMPA_Th']             = {'mod': 'MyExp2SynBB',        'tau1': 0.05,   'tau2': 5.3, 'e': 0}
netParams.synMechParams['GABAB_Th']            = {'mod': 'MyExp2SynBB',        'tau1': 3.5,    'tau2': 260.9,                  'e': -93} 
netParams.synMechParams['GABAA_Th']            = {'mod': 'MyExp2SynBB',        'tau1': 0.07,   'tau2': 18.2,                   'e': -80}
ESynMech_Th    = ['AMPA_Th', 'NMDA_Th']
PVSynMech_Th   = ['GABAA_Th']
NGFSynMech_Th  = ['GABAA_Th', 'GABAB_Th']

#------------------------------------------------------------------------------
# S1 Local connectivity parameters 
#------------------------------------------------------------------------------

#             if pre_new == 'L1_NGC_DA':
#                 pre = 'L1_NGC-DA'
#             else:
#                 pre = pre_new

#             if pre_new == 'L1_NGC_SA':
#                 pre = 'L1_NGC-SA'
#             else:
#                 pre = pre_new

#             if post_new == 'L1_NGC_DA':
#                 post = 'L1_NGC-DA'
#             else:
#                 post = pre_new

#             if post_new == 'L1_NGC_SA':
#                 post = 'L1_NGC-SA'
#             else:
#                 post = pre_new

contA = 0

if cfg.addConn:    
    for pre in Ipops+Epops:
        for post in Ipops+Epops:
            if float(connNumber[pre][post]) > 0:           
                # ------------------------------------------------------------------------------    
                #  2D distance prob rules
                # ------------------------------------------------------------------------------ 
                if int(float(d0[pre][post])) < 25:    # single fit
                    if 'exp' in connData['best_fit'][pre][post]:  # exponential            
                        prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s)' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post])      
                    else: # gaussian
                        prob = '%s*exp(-(dist_2D-%s)**2/(2*%s**2))*(dist_2D<%s)' % (a0g[pre][post],x0g[pre][post],l0g[pre][post],dfinal[pre][post])   
                        
                else:
                    if 'expl' in connData['best_fit'][pre][post]:  # exponential + linear interpolation [25:d0]
                        if int(float(d0[pre][post])) == 25:    #d0==25 -> exponential fit when dist_2D>25, else prob[0um:25um] = pmat[12.5]
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],float(pmat[12.5][pre][post]))
                        else:    #d0>25 -> exponential fit when dist_2D>d0, else prob[0um:d0] = linear interpolation [25:d0]
                            d01 = int(float(d0[pre][post]))
                            y1 = float(pmat[25][pre][post])
                            y2 = float(pmat[d01][pre][post])
                            x1 = 25
                            x2 = d01                   
                            angular = (y2 - y1)/(x2 - x1)
                            linear = y2 - x2*angular
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f * dist_2D + %f' % (a0mat[pre][post],lmat[pre][post],dfinal[pre][post],d0[pre][post],angular,linear)
                   
                    elif 'exp' in connData['best_fit'][pre][post]:  # exponential     
                        if float(pmat[12.5][pre][post]) > float(pmat[25][pre][post]):
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s)' % (a0e[pre][post],l0e[pre][post],dfinal[pre][post])
                        else:  
                            prob = '%s*exp(-dist_2D/%s)*(dist_2D<%s) if dist_2D > %s else %f' % (a0e[pre][post],l0e[pre][post],dfinal[pre][post],d0e[pre][post],float(pmat[12.5][pre][post]))      
                    
                    else: # gaussian
                        prob = '%s*exp(-(dist_2D-%s)**2/(2*%s**2))*(dist_2D<%s)' % (a0g[pre][post],x0g[pre][post],l0g[pre][post],dfinal[pre][post])             
                        
                # ------------------------------------------------------------------------------    
                # I -> I
                # ------------------------------------------------------------------------------
                if pre in Ipops:
                    if post in Ipops:                             
                        connID = ConnTypes[pre][post][0]                        
                        synMechType = 'S1_II_STP_Det_' + str(connID)   
                        contA+= 1
                        netParams.connParams['II_' + pre + '_' + post] = { 
                                        'preConds': {'pop': cfg.popLabelEl[pre]}, 
                                        'postConds': {'pop': cfg.popLabelEl[post]},
                                        'synMech': synMechType,
                                        'probability': prob,
                                        'weight': parameters_syn['gsyn',connID] * cfg.IIGain, 
                                        'synMechWeightFactor': cfg.synWeightFractionII,
                                        'delay': 'defaultDelay+dist_3D/propVelocity',
                                        'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                        'sec': 'spiny'}        
                # ------------------------------------------------------------------------------
                #  I -> E  # with ME conn diversity
                # ------------------------------------------------------------------------------
                if pre in Ipops:
                    if post in Epops:                                                       
                        cellpreList_A = []
                        cellpreList_B = []
                        cellpreList_C = []
                        connID_B = -1    
                        connID_C = -1                               
                        if 'SBC' in pre or 'LBC' in pre or 'NBC' in pre:                              
                            cellpost = cfg.popLabelEl[post][0]   
                            for npre,cellpre in enumerate(cfg.popLabelEl[pre]):   
                                premtype = pre[-3:]
                                preetype = cellpre[-3:]                                    
                                connID = connIEtype[premtype][preetype]                                     
                                if connID == ConnTypes[pre][post][0]:
                                    cellpreList_A.append(cellpre)    
                                elif connID == ConnTypes[pre][post][1]:
                                    cellpreList_B.append(cellpre)
                                    connID_B = ConnTypes[pre][post][1]
                                elif connID == ConnTypes[pre][post][2]:
                                    cellpreList_C.append(cellpre)
                                    connID_C = ConnTypes[pre][post][2]
                                else:
                                    print('ERROR')                                    
                        else:   
                            cellpreList_A = cfg.popLabelEl[pre]                              
                            
                        connID = ConnTypes[pre][post][0]                            
                        synMechType = 'S1_IE_STP_Det_' + str(connID)
                        
                        contA+= 1                          
                        netParams.connParams['IE_'+pre+'_'+post] = { 
                                    'preConds': {'pop': cellpreList_A}, 
                                    'postConds': {'pop': cfg.popLabelEl[post]},
                                    'synMech': synMechType,
                                    'probability': prob,
                                    'weight': parameters_syn['gsyn',connID] * cfg.IEGain, 
                                    'synMechWeightFactor': cfg.synWeightFractionIE,
                                    'delay': 'defaultDelay+dist_3D/propVelocity',
                                    'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                    'sec': 'spiny'}  
                

                        if connID_B >= 0:          
                            connID = connID_B
                            synMechType = 'S1_IE_STP_Det_' + str(connID)         
                            netParams.connParams['IE_'+pre+'_'+post+'_B'] = { 
                                        'preConds': {'pop': cellpreList_B}, 
                                        'postConds': {'pop': cfg.popLabelEl[post]},
                                        'synMech': synMechType,
                                        'probability': prob,
                                        'weight': parameters_syn['gsyn',connID] * cfg.IEGain, 
                                        'synMechWeightFactor': cfg.synWeightFractionIE,
                                        'delay': 'defaultDelay+dist_3D/propVelocity',
                                        'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                        'sec': 'spiny'}                       
                
                                
                            if connID_C >= 0:          
                                connID = connID_C
                                synMechType = 'S1_IE_STP_Det_' + str(connID)         
                                netParams.connParams['IE_'+pre+'_'+post+'_C'] = { 
                                            'preConds': {'pop': cellpreList_C}, 
                                            'postConds': {'pop': cfg.popLabelEl[post]},
                                            'synMech': synMechType,
                                            'probability': prob,
                                            'weight': parameters_syn['gsyn',connID] * cfg.IEGain, 
                                            'synMechWeightFactor': cfg.synWeightFractionIE,
                                            'delay': 'defaultDelay+dist_3D/propVelocity',
                                            'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                            'sec': 'spiny'}                       
                                
                #------------------------------------------------------------------------------   
                # E -> E
                #------------------------------------------------------------------------------
                if pre in Epops:
                    if post in Epops:    
                        connID = ConnTypes[pre][post][0]                        
                        synMechType = 'S1_EE_STP_Det_' + str(connID)   
                        contA+= 1   
                        netParams.connParams['EE_'+pre+'_'+post] = { 
                            'preConds': {'pop': cfg.popLabelEl[pre]}, 
                            'postConds': {'pop': cfg.popLabelEl[post]},
                            'synMech': synMechType,
                            'probability': prob, 
                            'weight': parameters_syn['gsyn',connID] * cfg.EEGain, 
                            'synMechWeightFactor': cfg.synWeightFractionEE,
                            'delay': 'defaultDelay+dist_3D/propVelocity',
                            'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                            'sec': 'spinyEE'}    
    
                #------------------------------------------------------------------------------               
                #  E -> I  with ME conn diversity
                #------------------------------------------------------------------------------   
                if pre in Epops:
                    if post in Ipops:                        
                        cellpostList_A = []
                        cellpostList_B = []
                        connID_B = -1                          
                        if ConnTypes[pre][post][0] == 131 or ConnTypes[pre][post][0] == 132: # EXCEPTIONS -> L6_IPC:L6_(DBC-LBC-NBC-SBC) and  L6_TPC_L:L6_(DBC-LBC-NBC-SBC)    
                            cellpostList_A = cfg.popLabelEl[post]     
                        elif 'LBC' in post or 'NBC' in post or 'BP' in post or 'DBC' in post or 'BTC' in post:    
                            cellpre = cfg.popLabelEl[pre][0]
                            for npost,cellpost in enumerate(cfg.popLabelEl[post]):                                
                                postmtype = post[-3:]
                                postetype = cellpost[-3:]
                                if 'BP' in postmtype:
                                    postmtype = post[-2:]       
                                connID = connEItype[postmtype][postetype]                                
                                if connID == ConnTypes[pre][post][0]:
                                    cellpostList_A.append(cellpost)    
                                elif connID == ConnTypes[pre][post][1]:
                                    cellpostList_B.append(cellpost)
                                    connID_B = ConnTypes[pre][post][1]
                                else:
                                    print('ERROR')                                
                        else:                           
                            cellpostList_A = cfg.popLabelEl[post]         
                             
                        connID = ConnTypes[pre][post][0]  

                        if 'DBC' in post or 'BTC' in post or 'MC' in post or 'BP' in post:  # steep Ca2+ dependence for connections between PC-distal targeting cell types (DBC, BTC, MC, BP)
                            synMechType = 'S1_EIdistal_STP_Det_' + str(connID)
                        else: # shallow dependence between PC-proximal targeting cell types (LBCs, NBCs, SBCs, ChC) + L1s and NGCs ????
                            synMechType = 'S1_EIproximal_STP_Det_' + str(connID)  

                        contA+= 1                                                              
                        netParams.connParams['EI_'+pre+'_'+post] = { 
                                        'preConds': {'pop': cfg.popLabelEl[pre]}, 
                                        'postConds': {'pop': cellpostList_A},
                                        'synMech': synMechType,
                                        'probability': prob, 
                                        'weight': parameters_syn['gsyn',connID] * cfg.EIGain, 
                                        'synMechWeightFactor': cfg.synWeightFractionEI,
                                        'delay': 'defaultDelay+dist_3D/propVelocity',
                                        'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                        'sec': 'spiny'}   

                        if connID_B >= 0:      

                            connID = connID_B

                            if 'DBC' in post or 'BTC' in post or 'MC' in post or 'BP' in post:  # steep Ca2+ dependence for connections between PC-distal targeting cell types (DBC, BTC, MC, BP)
                                synMechType = 'S1_EIdistal_STP_Det_' + str(connID)
                            else: # shallow dependence between PC-proximal targeting cell types (LBCs, NBCs, SBCs, ChC) + L1s and NGCs ????
                                synMechType = 'S1_EIproximal_STP_Det_' + str(connID)  

                            netParams.connParams['EI_'+pre+'_'+post+'_B'] = { 
                                            'preConds': {'pop': cfg.popLabelEl[pre]}, 
                                            'postConds': {'pop': cellpostList_B},
                                            'synMech': synMechType,
                                            'probability': prob, 
                                            'weight': parameters_syn['gsyn',connID] * cfg.EIGain, 
                                            'synMechWeightFactor': cfg.synWeightFractionEI,
                                            'delay': 'defaultDelay+dist_3D/propVelocity',
                                            'synsPerConn': int(synperconnNumber[pre][post]+0.5),
                                            'sec': 'spiny'}   


#------------------------------------------------------------------------------
# NetStim inputs to simulate Spontaneous synapses + background in S1 neurons - data from Rat
#------------------------------------------------------------------------------
SourcesNumber = 5 # for each post Mtype - sec distribution
synperNeuronStimI = connData['synperNeuronStimI']
synperNeuronStimE = connData['synperNeuronStimE']
GsynStimI = connData['GsynStimI']
GsynStimE = connData['GsynStimE']
   
if cfg.addStimSynS1:      
    for post in Ipops + Epops:

        synperNeuron = synperNeuronStimI[post]
        ratespontaneous = cfg.rateStimI
        for qSnum in range(SourcesNumber):
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
            netParams.stimSourceParams['StimSynS1_S_all_INH->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}

        synperNeuron = synperNeuronStimE[post]
        ratespontaneous = cfg.rateStimE
        for qSnum in range(SourcesNumber):
            ratesdifferentiation = (0.8 + 0.4*qSnum/(SourcesNumber-1)) * (synperNeuron*ratespontaneous)/SourcesNumber
            netParams.stimSourceParams['StimSynS1_S_all_EXC->' + post + '_' + str(qSnum)] = {'type': 'NetStim', 'rate': ratesdifferentiation, 'noise': 1.0}
            
    #------------------------------------------------------------------------------
    for post in Epops:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
                'conds': {'cellType': cfg.popLabelEl[post]}, 
                'synMech': 'AMPA', 
                'sec': 'spinyEE', 
                'weight': GsynStimE[post],
                'delay': 0.1}

    for post in Ipops:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_EXC->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_EXC->' + post + '_' + str(qSnum), 
                'synMech': 'AMPA', 
                'conds': {'cellType': cfg.popLabelEl[post]}, 
                'sec': 'spiny', 
                'weight': GsynStimE[post],
                'delay': 0.1}

    for post in Epops+Ipops:
        for qSnum in range(SourcesNumber):
            netParams.stimTargetParams['StimSynS1_T_all_INH->' + post + '_' + str(qSnum)] = {
                'source': 'StimSynS1_S_all_INH->' + post + '_' + str(qSnum), 
                'conds': {'cellType': cfg.popLabelEl[post]}, 
                'synMech': 'GABAA', 
                'sec': 'spiny', 
                'weight': GsynStimI[post],
                'delay': 0.1}


#------------------------------------------------------------------------------
# Description
#------------------------------------------------------------------------------
netParams.description = """ 
- Code based: M1 net, 6 layers, 7 cell types - v103
- v0 - insert cell diversity
- v1 - insert connection rules
- v2 - insert phys conn parameters
- v3 - ajust conn number
- v4 - NetStim inputs to simulate Spontaneous synapses + background in S1 neurons - data from Rat
"""
