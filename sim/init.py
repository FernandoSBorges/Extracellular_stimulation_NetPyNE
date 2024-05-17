"""
init.py

Starting script to run NetPyNE-basedS1 model.

Usage:
    python init.py # Run simulation, optionally plot a raster

MPI usage:
    mpiexec -n 4 nrniv -python -mpi init.py

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""

import matplotlib; matplotlib.use('Agg')  # to avoid graphics error in servers
from netpyne import sim
import neuron
import pickle, json
import numpy as np
import os


# cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')
cfg, netParams = sim.readCmdLineArgs()

sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations

## Load cells positions
with open('cells/spkTimes_v9_batch8_highgsynCT.pkl', 'rb') as fileObj: simData = pickle.load(fileObj)

cellsTags = simData['cellsTags']

# print(sim.rank,sim.net.cells[0].tags)

for i,metype in enumerate(sim.net.cells):

    if 'presyn' in metype.tags['pop']:
        ii = int(metype.tags['cellLabel'])        
        metype.tags['xnorm'] = cellsTags[ii]['xnorm']
        metype.tags['ynorm'] = cellsTags[ii]['ynorm']
        metype.tags['znorm'] = cellsTags[ii]['znorm']
        metype.tags['x'] = cellsTags[ii]['x']
        metype.tags['y'] = cellsTags[ii]['y']
        metype.tags['z'] = cellsTags[ii]['z']   

    else:
        ii2 = int(0.000001+(metype.tags['fraction']/(1/cfg.Nmorpho[metype.tags['pop']])))  

        ii = cfg.listmorphonumber[metype.tags['pop']][ii2]

        metype.tags['xnorm'] = cellsTags[ii]['xnorm']
        metype.tags['ynorm'] = cellsTags[ii]['ynorm']
        metype.tags['znorm'] = cellsTags[ii]['znorm']
        metype.tags['x'] = cellsTags[ii]['x']
        metype.tags['y'] = cellsTags[ii]['y']
        metype.tags['z'] = cellsTags[ii]['z']   

sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.net.defineCellShapes()

def collect_pt3d(self, section):
        """        collect the pt3d info, for each section
        """
        n3dsec = 0
        r3dsec = np.zeros(3)
        for sec in [sec for secName, sec in self.secs.items() if section in secName]:
            sec['hObj'].push()
            n3d = int(neuron.h.n3d())  # get number of n3d points in each section
            # print("get number of n3d points in each section",n3d)
            r3d = np.zeros((3, n3d))  # to hold locations of 3D morphology for the current section
            n3dsec += n3d

            for i in range(n3d):
                r3dsec[0] += neuron.h.x3d(i)
                r3dsec[1] += neuron.h.y3d(i)
                r3dsec[2] += neuron.h.z3d(i)
            
            neuron.h.pop_section()

        r3dsec /= n3dsec
        
        return r3dsec
    
def getSecsPos(self, secList):
        """        Get Secs position
        """
        x3d, y3d, z3d = [], [], []
        
        for secName in secList:
            # print(secName)
            r3dsec = collect_pt3d(self, secName)
            # print(secName, r3dsec)
            
            x3d.append(r3dsec[0])
            y3d.append(r3dsec[1])
            z3d.append(r3dsec[2])
            
        return x3d, y3d, z3d
    
def insert_v_ext(self, v_ext, t_ext):

    self.t_ext = neuron.h.Vector(t_ext)
    self.v_ext = []
    for v in v_ext:
        self.v_ext.append(neuron.h.Vector(v))
    
    # play v_ext into e_extracellular reference
    i = 0
    for secName, sec in self.secs.items(): 
        # print(secName,i)
        for seg in sec['hObj']:
            self.v_ext[i].play(seg._ref_e_extracellular, self.t_ext)
        i += 1
  
def make_extracellular_stimuli(acs_params, self, secList):
    """ Function to calculate and apply external potential """
    x0, y0, z0 = acs_params['position']
    ext_field = np.vectorize(lambda x, y, z: 1 / (4 * np.pi *
                                                  (acs_params['sigma'] * 
                                                   np.sqrt((x0 - x)**2 + (y0 - y)**2 + (z0 - z)**2) + skull_attenuation)))

    stimstart = acs_params['stimstart']
    stimend = acs_params['stimend']
    stimdif = stimend-stimstart

    # MAKING THE EXTERNAL FIELD
    n_tsteps = int(stimdif / cfg.dt + 1)
    n_start = int(stimstart/cfg.dt)
    n_end = int(stimend/cfg.dt + 1)
    t = np.arange(start=n_start, stop=n_end) * cfg.dt
    pulse = acs_params['amp'] * 1000. * \
        np.sin(2 * np.pi * acs_params['frequency'] * t / 1000)

    totnsegs = len(secList)    
    v_cell_ext = np.zeros((totnsegs, n_tsteps))    
    v_cell_ext[:, :] = ext_field(getSecsPos(metype,secList)[0], -1*np.array(getSecsPos(metype, secList)[1]), getSecsPos(metype,secList)[2]).reshape(totnsegs, 1) * pulse.reshape(1, n_tsteps)
    
    insert_v_ext(self, v_cell_ext, t)

    return v_cell_ext, self


# The parameters of the extracellular point current source
acs_params = {'position': [0.0, -1710.0, 0.0],  # um # y = [pia, bone]
              'amp': 400.,  # uA,
              'stimstart': 1000,  # ms
              'stimend': 3000,  # ms
              'frequency': 5,  # Hz
              'sigma': 0.57  # decay constant S/m
              }

skull_attenuation = 0.01*710 #conductivity of bone(S/m) * thickness of rat skull um

#Add extracellular stim
for c,metype in enumerate(sim.net.cells):
    if 'presyn' not in metype.tags['pop']:
        # print("\n", metype.tags)
        secList = [secs for secs in metype.secs.keys() if "pt3d" in metype.secs[secs]['geom']]
        # print(secList)
        v_cell_ext, cell = make_extracellular_stimuli(acs_params, metype,secList)


sim.runSim()                      			# run parallel Neuron simulation  
# Gather/save data option 1: standard
# sim.gatherData()                 			# gather spiking data and cell info from each node
# Gather/save data option 2: distributed saving across nodes 
sim.saveDataInNodes()
sim.gatherDataFromFiles()
sim.analysis.plotData()         			# plot spike raster etc

if sim.rank == 0:
    file_path = cfg.saveFolder+'/'+cfg.simLabel+'_node_data'
    print(file_path)

    if os.path.exists(file_path):
        print("The files exist")
        filelist = os.listdir(file_path)
        for f in filelist:
            os.remove(os.path.join(file_path, f))
            print(f"The file {f} has been deleted.")
    
sim.saveData()                    			# save params, cell info and sim output to file (pickle,mat,txt,etc)#
# sim.analysis.plotData()         			# plot spike raster etc


# sim.analysis.plotRaster(include=cfg.S1cells, timeRange=[1000,cfg.duration], orderBy='gid', orderInverse=True, figSize=(18, 12), fontSize=9, dpi=300, saveFig=True, showFig=False)
# sim.analysis.plotRaster(timeRange=[0,cfg.duration], orderBy='gid', orderInverse=True, labels=None, popRates=False, lw=1, marker='.', markerSize=2, figSize=(18, 12), fontSize=9, dpi=300, saveFig=True, showFig=False)
# sim.analysis.plotTraces(include=cfg.recordCells, overlay=True, oneFigPer='cell', figSize=(12, 4), fontSize=7, saveFig=True)
#sim.analysis.plotTraces(include=cfg.recordCells, overlay=False, oneFigPer='trace', figSize=(18, 12), fontSize=9, saveFig=True)
# features = ['numConns','convergence']
# groups =['pop']
# for feat in features:
#    for group in groups:
#        sim.analysis.plotConn(includePre=['L1_DAC_cNA','L23_MC_cAC','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], includePost=['L1_DAC_cNA','L23_MC_cAC','L4_SS_cAD','L4_NBC_cNA','L5_TTPC2_cAD', 'L5_LBC_cNA', 'L6_TPC_L4_cAD', 'L6_LBC_cNA', 'ss_RTN_o', 'ss_RTN_m', 'ss_RTN_i', 'VPL_sTC', 'VPM_sTC', 'POm_sTC_s1'], feature=feat, groupBy=group, figSize=(24,24), saveFig=True, orderBy='gid', graphType='matrix', fontSize=18, saveData='../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' + group + '_' + feat+ '_matrix.json')

# sim.analysis.plotLFP(**{'plots': ['timeSeries'], 'electrodes': [0,1,2,3], 'timeRange': [150, cfg.duration], 'maxFreq':80, 'figSize': (16,8), 'saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' +'LFP1.png', 'showFig': False})
# sim.analysis.plotLFP(**{'plots': ['timeSeries'], 'electrodes': [4,5,6,7], 'timeRange': [0, 300], 'maxFreq':80, 'figSize': (16,8), 'saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' +'LFP2', 'showFig': False})
# sim.analysis.plotLFP(**{'plots': ['timeSeries'], 'electrodes': [8,9,10,11], 'timeRange': [0, 300], 'maxFreq':80, 'figSize': (16,8), 'saveFig': '../data/'+cfg.simLabel[0:9]+'/'+cfg.simLabel + '_' +'LFP3', 'showFig': False})
