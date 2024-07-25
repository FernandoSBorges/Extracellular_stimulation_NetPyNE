"""
cfg.py 

Simulation configuration for S1-thalamus model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle
import os
import numpy as np

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

cfg.simType='S1_stim'
cfg.coreneuron = False

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 5000.0 ## Duration of the sim, in ms  
cfg.dt = 0.05
cfg.seeds = {'cell': 4322, 'conn': 4322, 'stim': 4322, 'loc': 4322} 
cfg.hParams = {'celsius': 34, 'v_init': -71}  
cfg.verbose = False
cfg.createNEURONObj = True
cfg.createPyStruct = True  
cfg.cvode_active = False
cfg.cvode_atol = 1e-6
cfg.cache_efficient = True
cfg.printRunTime = 0.1

cfg.includeParamsLabel = False
cfg.printPopAvgRates = True
cfg.checkErrors = False

#------------------------------------------------------------------------------
# Cells
#------------------------------------------------------------------------------
# TO DEBUG - import and simulate only the Cell soma (to study only the Net)
cfg.reducedtest = False    

# TO DEBUG - Create only 5 Cells for each MEtype in S1
cfg.oneCellperMEtypeS1 = False 


cfg.rootFolder = os.getcwd()

# Load cells info from previously saved using netpyne (False: load from HOC BBP files, slower)
cfg.loadcellsfromJSON = True

cfg.poptypeNumber = 55
cfg.celltypeNumber = 207

cfg.cao_secs = 1.2

cfg.use_frac = {} # use[invivo] = cfg.use_frac * use[invitro]

cfg.use_frac['EIproximal'] = 0.75 # shallow dependence between PC-proximal targeting cell types (LBCs, NBCs, SBCs, ChC)
cfg.use_frac['Inh'] = 0.50 # Pathways that had not been studied experimentally were assumed to have an intermediate level of dependence
cfg.use_frac['EE'] = 0.25 # steep Ca2+ dependence for connections between PC-PC and PC-distal targeting cell types (DBC, BTC, MC, BP)
cfg.use_frac['EIdistal'] = 0.25 

#------------------------------------------------------------------------------  
#------------------------------------------------------------------------------  
# S1 Cells
# Load 55 Morphological Names and Cell pop numbers -> L1:6 L23:10 L4:12 L5:13 L6:14
# Load 207 Morpho-electrical Names used to import the cells from 'cell_data/' -> L1:14 L23:43 L4:46 L5:52 L6:52
# Create [Morphological,Electrical] = number of cell metype in the sub-pop

with open('cells/S1-cells-distributions-Human.txt') as mtype_file:
    mtype_content = mtype_file.read()       

cfg.popNumber = {}
cfg.cellNumber = {} 
cfg.popLabel = {} 
popParam = []
cellParam = []
cfg.meParamLabels = {} 
cfg.popLabelEl = {} 
cfg.cellLabel = {}

for line in mtype_content.split('\n')[:-1]:
    cellname, mtype, etype, n, m = line.split()
    metype = mtype + '_' + etype[0:3]
    cfg.cellNumber[metype] = int(n)
    cfg.popLabel[metype] = mtype
    cfg.popNumber[mtype] = int(m)
    cfg.cellLabel[metype] = cellname[0:-3]

    if mtype not in popParam:
        popParam.append(mtype)
        cfg.popLabelEl[mtype] = [] 
               
    cfg.popLabelEl[mtype].append(metype)
    
    cellParam.append(mtype + '_' + etype[0:3])
    
cfg.S1pops = popParam[0:55]
cfg.S1cells = cellParam[0:207]

#------------------------------------------------------------------------------  
cfg.popParamLabels = popParam
cfg.cellParamLabels = cellParam


Epops = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']

cfg.Ecells = [] 
for metype in cfg.S1cells: # metype      
    mtype = cfg.popLabel[metype]            
    if mtype in Epops:  
        cfg.Ecells.append(metype)      

# subPopLabels = ['L1_DAC','L1_DLAC', 'L1_HAC','L1_SLAC', #'L1_NGC_DA','L1_NGC_SA',
#  'L23_PC','L23_MC','L23_SBC', 'L23_BP','L23_BTC','L23_ChC','L23_DBC','L23_LBC','L23_NBC','L23_NGC',
#  'L4_PC','L4_SP','L4_SS','L4_SBC','L4_MC','L4_BP','L4_BTC','L4_ChC','L4_DBC','L4_LBC','L4_NBC','L4_NGC',
#  'L5_TTPC2','L5_STPC','L5_TTPC1','L5_UTPC','L5_SBC','L5_MC', 'L5_BP','L5_BTC','L5_ChC','L5_DBC','L5_LBC','L5_NBC','L5_NGC',
#  'L6_TPC_L4','L6_TPC_L1','L6_UTPC','L6_IPC','L6_BPC','L6_SBC','L6_MC','L6_BP','L6_BTC','L6_ChC','L6_DBC','L6_LBC','L6_NBC','L6_NGC']

# subPopLabels = cfg.S1pops[6:28] # from 0 to 55 is full S1 -> L1:6 L23:10 L4:12 L5:13 L6:14
# subPopLabels = Epops[0:12]

# subPopLabels = ['L1_DAC',
#  'L23_PC','L23_SBC',
#  'L4_PC','L4_SBC',
#  'L5_TTPC2','L5_SBC', 
#  'L6_TPC_L4','L6_SBC']

# subPopLabels = ['L1_DAC','L1_DLAC',
#  'L23_PC','L23_SBC','L23_LBC','L23_MC',
#  'L4_PC','L4_SP','L4_SS','L4_SBC','L4_LBC','L4_MC',
#  'L5_TTPC2','L5_STPC','L5_TTPC1','L5_UTPC','L5_SBC','L5_LBC','L5_MC', 
#  'L6_TPC_L4','L6_TPC_L1','L6_UTPC','L6_IPC','L6_BPC','L6_SBC','L6_LBC','L6_MC'] 

subPopLabels = ['L1_DAC','L1_DLAC',#'L1_HAC','L1_NGC_DA','L1_NGC_SA','L1_SLAC',
 'L23_PC','L23_MC','L23_SBC', #'L23_BP','L23_BTC','L23_ChC','L23_DBC','L23_LBC','L23_NBC','L23_NGC',
 'L4_PC','L4_SBC','L4_MC', # 'L4_BP','L4_BTC','L4_ChC','L4_DBC','L4_LBC','L4_NBC','L4_NGC','L4_SP','L4_SS',
 'L5_TTPC2','L5_SBC','L5_MC', #'L5_BP','L5_BTC','L5_ChC','L5_DBC','L5_LBC','L5_NBC','L5_NGC','L5_STPC','L5_TTPC1','L5_UTPC',
 'L6_TPC_L4','L6_SBC','L6_MC'] #,'L6_BPC','L6_BP','L6_BTC','L6_ChC','L6_DBC','L6_IPC','L6_LBC','L6_NBC','L6_NGC','L6_TPC_L1','L6_UTPC']

#------------------------------------------------------------------------------  
cfg.S1pops = subPopLabels
cfg.S1cells = []
for metype in cfg.cellParamLabels:
    if cfg.popLabel[metype] in subPopLabels:        
        cfg.S1cells.append(metype)
        
cfg.thalamicpops = []

cfg.popParamLabels = cfg.S1pops
cfg.cellParamLabels = cfg.S1cells

cfg.cellNumber[metype]
cfg.popNumber[cfg.popLabel[metype]]

#------------------------------------------------------------------------------  
# for metype in cfg.cellParamLabels:
#     print(metype,cfg.cellNumber[metype],cfg.popLabel[metype],cfg.popNumber[cfg.popLabel[metype]])   
    
#------------------------------------------------------------------------------  
## Change popNumber
#------------------------------------------------------------------------------  
# cfg.cellNumber[metype] = 25
# cfg.popNumber[cfg.popLabel[metype]] = 25
# #------------------------------------------------------------------------------  
# for metype in cfg.cellParamLabels:
#     print(metype,cfg.cellNumber[metype],cfg.popLabel[metype],cfg.popNumber[cfg.popLabel[metype]])   

#--------------------------------------------------------------------------
# Recording 
#--------------------------------------------------------------------------
cfg.allpops = cfg.cellParamLabels
cfg.cellsrec = 1
if cfg.cellsrec == 0:  cfg.recordCells = cfg.allpops # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in cfg.allpops] # record one cell of each pop
elif cfg.cellsrec == 2: # record one cell of each cellMEtype for Epops
    cfg.recordCells = []
    for metype in cfg.cellParamLabels:
        if metype in cfg.Ecells:
            for numberME in range(5):
                cfg.recordCells.append((metype,numberME))
        else:
            cfg.recordCells.append((metype,0))
            cfg.recordCells.append((metype,1))

cfg.recordTraces = {'V_soma': {'sec':'soma_0', 'loc':0.5, 'var':'v'},
                    'V_axon_0': {'sec':'axon_0', 'loc':0.5, 'var':'v'},
                    'V_Myelin_0': {'sec':'Myelin_0', 'loc':0.5, 'var':'v'},
                    'V_Myelin_10': {'sec':'Myelin_10', 'loc':0.5, 'var':'v'},
                    'V_Node_0': {'sec':'Node_0', 'loc':0.5, 'var':'v'},
                    'V_Node_10': {'sec':'Node_10', 'loc':0.5, 'var':'v'},
                    # 'V_Unmyelin_0': {'sec':'Unmyelin_0', 'loc':0.5, 'var':'v'},
                    # 'V_Unmyelin_10': {'sec':'Unmyelin_10', 'loc':0.5, 'var':'v'},
                    # 'V_apic_0': {'sec':'apic_0', 'loc':0.5, 'var':'v'},
                    # 'V_apic_5': {'sec':'apic_5', 'loc':0.5, 'var':'v'},
                    # 'V_apic_95': {'sec':'apic_95', 'loc':0.5, 'var':'v'},                
                    # 'V_dend_5': {'sec':'dend_5', 'loc':0.5, 'var':'v'},
                    # 'V_dend_65': {'sec':'dend_65', 'loc':0.5, 'var':'v'},
                    'V_dend_10': {'sec':'dend_10', 'loc':0.5, 'var':'v'},
                    }

cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 0.1           

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'tms_test_batch0'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = True         	## Save pkl file
cfg.saveJson = False	           	## Save json file
cfg.saveDataInclude = ['simData', 'simConfig', 'netParams', 'net'] ## , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			
cfg.saveCellConns = False	

#------------------------------------------------------------------------------
# Analysis and plotting 
# ------------------------------------------------------------------------------
cfg.analysis['plotRaster'] = {'include': cfg.allpops, 'saveFig': True, 'showFig': False, 'orderInverse': True, 'timeRange': [0,cfg.duration], 'figSize': (24,24), 'fontSize':6, 'lw': 4, 'markerSize':4, 'marker': '.', 'dpi': 300} 
cfg.analysis['plot2Dnet']   = {'include': cfg.allpops, 'saveFig': True, 'showConns': False, 'figSize': (18,18), 'fontSize':8}   # Plot 2D cells xy
cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'oneFigPer': 'cell', 'overlay': True, 'timeRange': [0,cfg.duration], 'saveFig': True, 'showFig': False, 'figSize':(18,12)}
# cfg.analysis['plot2Dfiring']={'saveFig': True, 'figSize': (24,24), 'fontSize':16}
# cfg.analysis['plotConn'] = {'includePre': cfg.allpops, 'includePost': cfg.allpops, 'feature': 'numConns', 'groupBy': 'pop', 'figSize': (24,24), 'saveFig': True, 'orderBy': 'gid', 'graphType': 'matrix', 'saveData':'../data/v5_batch0/v5_batch0_matrix_numConn.json', 'fontSize': 18}
# cfg.analysis['plotConn'] = {'includePre': ['L1_DAC_cNA','L23_MC_cAC','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], 'includePost': ['L1_DAC_cNA','L23_MC_cAC','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], 'feature': 'convergence', 'groupBy': 'pop', 'figSize': (24,24), 'saveFig': True, 'orderBy': 'gid', 'graphType': 'matrix', 'fontSize': 18}
# cfg.analysis['plot2Dnet']   = {'include': ['L5_LBC', 'VPM_sTC', 'POm_sTC_s1'], 'saveFig': True, 'showConns': True, 'figSize': (24,24), 'fontSize':16}   # Plot 2D net cells and connections
# cfg.analysis['plotShape'] = {'includePre': cfg.recordCells, 'includePost': cfg.recordCells, 'showFig': False, 'includeAxon': False, 
                            # 'showSyns': False, 'saveFig': True, 'dist': 0.55, 'cvar': 'voltage', 'figSize': (24,12), 'dpi': 600}

#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
soma_area_scaling_factor = 2.453
axon_diameter_scaling_factor = 2.453
main_axon_diameter_scaling_factor = 1
apic_diameter_scaling_factor = 1.876
dend_diameter_scaling_factor = 1.946
dend_length_scaling_factor = 1.17

L25_human = 950 + 380 + 700
L25_Rat = 502 + 190 + 525
# print(L25_human,L25_Rat,L25_human/L25_Rat)
# print("Human_height",2082*L25_human/L25_Rat)
Human_Rat_height_ratio = 1.668
Human_height = 3472.85

cfg.scale = 1.0 # reduce size
# cfg.sizeY = 2082.0*dend_length_scaling_factor
cfg.sizeY = 3500.0

cfg.scaleDensity = 0.5 # Number of cells = 31346 if 1.0

cfg.keepdensity = True

if cfg.keepdensity:
    cfg.sizeX = 420.0*np.sqrt(cfg.scaleDensity) # r = 210 um and hexagonal side length = 230.9 um -> 
    cfg.sizeZ = 420.0*np.sqrt(cfg.scaleDensity)
else:
    cfg.sizeX = 420.0 # keep the original radius
    cfg.sizeZ = 420.0

#------------------------------------------------------------------------------
# Spontaneous synapses + background - data from Rat
#------------------------------------------------------------------------------
cfg.addStimSynS1 = True
cfg.rateStimE = 9.0
cfg.rateStimI = 9.0

#------------------------------------------------------------------------------
# Connectivity
#------------------------------------------------------------------------------
## S1->S1
cfg.addConn = True

cfg.synWeightFractionEE = [1.0, 1.0] # E -> E AMPA to NMDA ratio
cfg.synWeightFractionEI = [1.0, 1.0] # E -> I AMPA to NMDA ratio
cfg.synWeightFractionII = [1.0, 1.0]  # I -> I GABAA to GABAB ratio
cfg.synWeightFractionIE = [1.0, 1.0]  # I -> E GABAA to GABAB ratio
cfg.EEGain = 1.0
cfg.EIGain = 1.0
cfg.IIGain = 1.0
cfg.IEGain = 1.0


#------------------------------------------------------------------------------
# External Stimulation
#------------------------------------------------------------------------------

cfg.addExternalStimulation = True

# The parameters of the extracellular point current source
cfg.acs_params = {'position': [0.0, -1710.0, 0.0],  # um # y = [pia, bone]
              'amp': -1250.,  # uA,
              'stimstart': 300,  # ms
              'stimend': 400.0,  # ms
              'frequency': 5,  # Hz
              'sigma': 0.57  # decay constant S/m
              }

cfg.tms_params = dict(
    freq_Hz=30.,
    duration_ms=cfg.duration,
    pulse_resolution_ms=cfg.dt,
    stim_start_ms=2000.,
    stim_end_ms=3000.,
    ef_amp_V_per_m=80.,
    width_ms=1.,
    pshape="Sine",
    decay_rate_percent_per_mm=10,
    E_field_dir=[-1, -1, -1],
    decay_dir=[0, 0, -1],
    ref_point_um=[0, 0, 0],
)