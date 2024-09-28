"""
cfg.py 

Simulation configuration for human cortex model (using NetPyNE)
This file has sim configs as well as specification for parameterized values in netParams.py 

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

from netpyne import specs
import pickle
import os
import numpy as np
import pandas as pd

cfg = specs.SimConfig()  

#------------------------------------------------------------------------------
#
# SIMULATION CONFIGURATION
#
#------------------------------------------------------------------------------

cfg.simType='h01_stim'
cfg.coreneuron = False

#------------------------------------------------------------------------------
# Run parameters
#------------------------------------------------------------------------------
cfg.duration = 5000.0 ## Duration of the sim, in ms  
cfg.dt = 0.05
cfg.seeds = {'cell': 4321, 'conn': 4321, 'stim': 4321, 'loc': 4321} 
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

cfg.rootFolder = os.getcwd()

# Load cells info from previously saved using netpyne (False: load from HOC BBP files, slower)
cfg.loadcellsfromJSON = True

cfg.poptypeNumber = 55
cfg.celltypeNumber = 134

cfg.cao_secs = 1.2

cfg.use_frac = {} # use[invivo] = cfg.use_frac * use[invitro]

cfg.use_frac['EIproximal'] = 0.75 # shallow dependence between PC-proximal targeting cell types (LBCs, NBCs, SBCs, ChC)
cfg.use_frac['Inh'] = 0.50 # Pathways that had not been studied experimentally were assumed to have an intermediate level of dependence
cfg.use_frac['EE'] = 0.25 # steep Ca2+ dependence for connections between PC-PC and PC-distal targeting cell types (DBC, BTC, MC, BP)
cfg.use_frac['EIdistal'] = 0.25 

#------------------------------------------------------------------------------  
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

    # print(cellname, mtype, etype, n, m)
    
cfg.S1pops = popParam
cfg.S1cells = cellParam

cfg.popParamLabels = popParam
cfg.cellParamLabels = cellParam


#------------------------------------------------------------------------------  
Epops = ['L23_PC', 'L4_PC', 'L4_SS', 'L4_SP', 
             'L5_TTPC1', 'L5_TTPC2', 'L5_STPC', 'L5_UTPC',
             'L6_TPC_L1', 'L6_TPC_L4', 'L6_BPC', 'L6_IPC', 'L6_UTPC']

cfg.Ecells = [] 
cfg.h01_mtype = {}

for metype in cfg.S1cells: # metype    

    mtype = cfg.popLabel[metype]

    if mtype[1] == '2':
        layer = 'L23'      
    else:
        layer = mtype[0:2]         

    if mtype in Epops:  
        cfg.h01_mtype[layer+'E'] = 0
        cfg.Ecells.append(metype)      
    else:
        cfg.h01_mtype[layer+'I'] = 0
        
for metype in cfg.S1cells: # metype    

    mtype = cfg.popLabel[metype]

    if mtype[1] == '2':
        layer = 'L23'      
    else:
        layer = mtype[0:2]         

    if mtype in Epops:  
        cfg.h01_mtype[layer+'E'] = cfg.h01_mtype[layer+'E'] + cfg.cellNumber[metype]
    else:
        cfg.h01_mtype[layer+'I'] = cfg.h01_mtype[layer+'I'] + cfg.cellNumber[metype]
        
#------------------------------------------------------------------------------
# Network 
#------------------------------------------------------------------------------
# Aberra neuron model
soma_area_scaling_factor = 2.453
axon_diameter_scaling_factor = 2.453
main_axon_diameter_scaling_factor = 1
apic_diameter_scaling_factor = 1.876
dend_diameter_scaling_factor = 1.946
dend_length_scaling_factor = 1.17

# h01 L1-L6 rotated
Human_height = 2556.35

cfg.cylinderRadius_h01 = 300.0 # 

cfg.scale = 1.0 # reduce size
cfg.sizeX = 600.0 
cfg.sizeY = 3200.00
cfg.sizeZ = 600.0
cfg.scaleDensity = 1.0 # Number of cells = 5527 if cfg.sizeZ = cfg.sizeX = 600.0 um


#------------------------------------------------------------------------------  
nodes_new = pd.read_csv('../data/cell_positions_h01_rotated.csv')

h01type, h01Number = np.unique(nodes_new[nodes_new['distance2Dcenter'] < cfg.cylinderRadius_h01]['mtype'].values, return_counts=True)

print("Cell Number per layer for h01 considering all cells in a cylinder of Radius_um =",cfg.cylinderRadius_h01,h01type, h01Number)

cfg.h01_ratio_number = {}
cfg.List_h01 = {}
# print('layer','\t','h01','\t','S1','\t', "perc\n")
for h01N, h01t in enumerate(['L1I', 'L23E', 'L23I', 'L4E', 'L4I', 'L5E', 'L5I', 'L6E', 'L6I']):
    cfg.h01_ratio_number[h01t] = h01Number[h01N]/cfg.h01_mtype[h01t]
    # print(h01t,'\t',h01Number[h01N],'\t',cfg.h01_mtype[h01t],'\t',"%.2f" % (100.0*cfg.h01_ratio_number[h01t]))
    cfg.List_h01[h01t] = []

for mtype in cfg.S1pops:

    if mtype[1] == '2':
        layer = 'L23'      
    else:
        layer = mtype[0:2]         

    if mtype in Epops:    
        h01t = layer+'E'
    else:
        h01t = layer+'I'    

    cfg.List_h01[h01t].append(mtype)

# cfg.List_h01

#------------------------------------------------------------------------------  
cfg.popNumber_new = {}
number=0
listnew = []

for h01N, h01t in enumerate(['L1I', 'L23E', 'L23I', 'L4E', 'L4I', 'L5E', 'L5I', 'L6E', 'L6I']):
# for h01t in cfg.List_h01.keys():
    number=0
    # print()

    for mtype in cfg.List_h01[h01t]:

        if mtype[1] == '2':
            layer = 'L23'      
        else:
            layer = mtype[0:2]         


        cfg.popNumber_new[mtype] = int(0.5 + cfg.popNumber[mtype] * cfg.h01_ratio_number[h01t]) 
        
        number+=cfg.popNumber_new[mtype]

        # print(mtype, cfg.popNumber[mtype], cfg.popNumber_new[mtype], h01Number[h01N], number)
        
    cfg.popNumber_new[mtype] = cfg.popNumber_new[mtype] + h01Number[h01N] - number
    number+= h01Number[h01N] - number

    # print(" ",mtype, cfg.popNumber[mtype], cfg.popNumber_new[mtype], h01Number[h01N], number)

#------------------------------------------------------------------------------  
cfg.cellNumber_new = {}

for mtype in cfg.S1pops:
    # print()
    number=0
    for cellEl in range(min(np.size(cfg.popLabelEl[mtype]),cfg.popNumber_new[mtype])):

        metype = cfg.popLabelEl[mtype][cellEl]

        if cfg.popNumber_new[mtype] <= np.size(cfg.popLabelEl[mtype]):
            cfg.cellNumber_new[metype] = int(1.0)
        else:
            cfg.cellNumber_new[metype] = int(0.5 + cfg.cellNumber[metype]*cfg.popNumber_new[mtype]/cfg.popNumber[mtype])

        number+=cfg.cellNumber_new[metype]

        # print(int(number), metype, cfg.cellNumber[metype], mtype, cfg.popNumber_new[mtype],cfg.cellNumber_new[metype])

    cfg.cellNumber_new[metype] = cfg.cellNumber_new[metype] + cfg.popNumber_new[mtype] - number
    # print(" ", int(number), metype, cfg.cellNumber[metype], mtype, cfg.popNumber_new[mtype],cfg.cellNumber_new[metype])

    if cfg.cellNumber_new[metype] == 0:

        if 'L6_DBC_cAC' in cfg.cellNumber_new.keys():
            if  cfg.cellNumber_new['L6_DBC_cAC'] == 2 and cfg.cellNumber_new['L6_DBC_cNA'] == 0:
                print(" ########### \n fixing cellNumber_new = 0 for L6_DBC_cNA")
                cfg.cellNumber_new['L6_DBC_cAC'] = 1
                cfg.cellNumber_new['L6_DBC_cNA'] = 1
                print(" ", int(number), metype, cfg.cellNumber[metype], mtype, cfg.popNumber_new[mtype],cfg.cellNumber_new[metype])

        if cfg.cellNumber_new[metype] == 0:
            print("fix this like above")
            break
#------------------------------------------------------------------------------  
number=0
for mtype in cfg.S1pops:
    for cellEl in range(np.size(cfg.popLabelEl[mtype])):

        metype = cfg.popLabelEl[mtype][cellEl]

        try:
            number+=cfg.cellNumber_new[metype]
        except:
            print(metype, "not inclued for this size")

# print('Cell Number =', number)


cfg.popLabelEl = {} 
for metype in cfg.cellNumber_new.keys():

    mtype = cfg.popLabel[metype] 

    if mtype not in cfg.popLabelEl.keys():
        cfg.popLabelEl[mtype] = [] 
               
    cfg.popLabelEl[mtype].append(metype)
    

number=0
for mtype in cfg.S1pops:
    for cellEl in range(np.size(cfg.popLabelEl[mtype])):

        metype = cfg.popLabelEl[mtype][cellEl]

        try:
            number+=cfg.cellNumber_new[metype]
        except:
            print(metype, "not inclued for this size")

# print('Cell Number =', number)

#--------------------------------------------------------------------------
cfg.S1pops = list(cfg.popNumber_new.keys())
cfg.S1cells = list(cfg.cellNumber_new.keys())

cfg.popParamLabels = list(cfg.popNumber_new.keys())
cfg.cellParamLabels = list(cfg.cellNumber_new.keys())

#--------------------------------------------------------------------------

# TO DEBUG - Create only 5 Cells for each MEtype
cfg.oneCellperMEtypeS1 = False 

#------------------------------------------------------------------------------  
# TO DEBUG - Create only one Cell per MEtype
if cfg.oneCellperMEtypeS1:
    cfg.popNumber_new = {}
    for mtype in cfg.S1pops:
        cfg.popNumber_new[mtype] = 0

    for metype in cfg.cellNumber_new.keys():

        mtype = cfg.popLabel[metype] 

        if mtype in Epops:
            cfg.popNumber_new[mtype] = cfg.popNumber_new[mtype] + 5
            cfg.cellNumber_new[metype] = 5
        else:
            cfg.popNumber_new[mtype] = cfg.popNumber_new[mtype] + 1
            cfg.cellNumber_new[metype] = 1

#--------------------------------------------------------------------------
# Recording 
#--------------------------------------------------------------------------
cfg.allpops = cfg.cellParamLabels
cfg.cellsrec = 2
if cfg.cellsrec == 0:  cfg.recordCells = cfg.allpops # record all cells
elif cfg.cellsrec == 1: cfg.recordCells = [(pop,0) for pop in cfg.allpops] # record one cell of each pop
elif cfg.cellsrec == 2: # record one cell of each cellMEtype for Epops
    cfg.recordCells = []
    for metype in cfg.cellParamLabels:
        if metype not in cfg.Ecells:
            cfg.recordCells.append((metype,0))
            cfg.recordCells.append((metype,1))
        else:
            numberME = 0
            diference = cfg.cellNumber_new[metype] - 5.0*int(cfg.cellNumber_new[metype]/5.0)
            
            for number in range(5):            
                cfg.recordCells.append((metype,numberME))
                
                if number < diference:              
                    numberME+=int(np.ceil(cfg.cellNumber_new[metype]/5.0))  
                else:
                    numberME+=int(cfg.cellNumber_new[metype]/5.0)


cfg.recordTraces = {'V_soma': {'sec':'soma_0', 'loc':0.5, 'var':'v'},
                    # 'V_axon_0': {'sec':'axon_0', 'loc':0.5, 'var':'v'},
                    # 'V_Myelin_0': {'sec':'Myelin_0', 'loc':0.5, 'var':'v'},
                    # 'V_Myelin_10': {'sec':'Myelin_10', 'loc':0.5, 'var':'v'},
                    # 'V_Node_0': {'sec':'Node_0', 'loc':0.5, 'var':'v'},
                    # 'V_Node_10': {'sec':'Node_10', 'loc':0.5, 'var':'v'},
                    # 'V_Unmyelin_0': {'sec':'Unmyelin_0', 'loc':0.5, 'var':'v'},
                    # 'V_Unmyelin_10': {'sec':'Unmyelin_10', 'loc':0.5, 'var':'v'},
                    # 'V_apic_0': {'sec':'apic_0', 'loc':0.5, 'var':'v'},
                    # 'V_apic_5': {'sec':'apic_5', 'loc':0.5, 'var':'v'},
                    # 'V_apic_95': {'sec':'apic_95', 'loc':0.5, 'var':'v'},                
                    # 'V_dend_5': {'sec':'dend_5', 'loc':0.5, 'var':'v'},
                    # 'V_dend_65': {'sec':'dend_65', 'loc':0.5, 'var':'v'},
                    # 'V_dend_10': {'sec':'dend_10', 'loc':0.5, 'var':'v'},
                    }

cfg.recordStim = False			
cfg.recordTime = False  		
cfg.recordStep = 0.05           

#------------------------------------------------------------------------------
# Saving
#------------------------------------------------------------------------------
cfg.simLabel = 'v0_batch0'
cfg.saveFolder = '../data/'+cfg.simLabel
# cfg.filename =                	## Set file output name
cfg.savePickle = False         	## Save pkl file
cfg.saveJson = False	           	## Save json file
cfg.saveDataInclude = ['simConfig'] # ['simData', 'simConfig', 'netParams', 'net'] ## , 'simConfig', 'netParams'
cfg.backupCfgFile = None 		##  
cfg.gatherOnlySimData = False	##  
cfg.saveCellSecs = False			
cfg.saveCellConns = False	

#------------------------------------------------------------------------------
# Analysis and plotting 
# ------------------------------------------------------------------------------
cfg.analysis['plotRaster'] = {'include': cfg.allpops, 'saveFig': True, 'showFig': False, 'orderInverse': True, 'timeRange': [0.0,cfg.duration], 'saveData': True, 'figSize': (24,24), 'fontSize':12, 'lw': 8, 'markerSize':8, 'marker': '.', 'dpi': 300} 
cfg.analysis['plot2Dnet']   = {'include': cfg.allpops, 'saveFig': True, 'showConns': False, 'figSize': (12,18), 'fontSize':8}   # Plot 2D cells xy
cfg.analysis['plotTraces'] = {'include': cfg.recordCells, 'oneFigPer': 'trace', 'overlay': True, 'timeRange': [0.0,cfg.duration], 'saveData': True, 'saveFig': True, 'showFig': False, 'figSize':(18,12)}
# cfg.analysis['plot2Dfiring']={'saveFig': True, 'figSize': (24,24), 'fontSize':16}
# cfg.analysis['plotConn'] = {'includePre': cfg.allpops, 'includePost': cfg.allpops, 'feature': 'numConns', 'groupBy': 'pop', 'figSize': (48,48), 'saveFig': True, 'orderBy': 'gid', 'graphType': 'matrix', 'saveData':True, 'fontSize': 18}
# cfg.analysis['plot2Dnet']   = {'include': ['L5_LBC', 'VPM_sTC', 'POm_sTC_s1'], 'saveFig': True, 'showConns': True, 'figSize': (24,24), 'fontSize':16}   # Plot 2D net cells and connections
# cfg.analysis['plotShape'] = {'includePre': cfg.recordCells, 'includePost': cfg.recordCells, 'showFig': False, 'includeAxon': False, 
                            # 'showSyns': False, 'saveFig': True, 'dist': 0.55, 'cvar': 'voltage', 'figSize': (24,12), 'dpi': 600}

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
    ef_amp_V_per_m=60.,
    width_ms=1.0,
    pshape="Sine",
    decay_rate_percent_per_mm=10,
    E_field_dir=[-1, -1, -1],
    decay_dir=[0, 0, -1],
    ref_point_um=[0, 0, 0],
)

#------------------------------------------------------------------------------
# Current inputs 
#------------------------------------------------------------------------------

cfg.addIClamp = False
cfg.IClamp_nA = 0.310
