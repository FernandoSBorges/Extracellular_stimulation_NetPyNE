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

from stimulation import make_extracellular_stimuli
from tms_tools import apply_tms

tms = True

# cfg, netParams = sim.readCmdLineArgs(simConfigDefault='cfg.py', netParamsDefault='netParams.py')
cfg, netParams = sim.readCmdLineArgs()

sim.initialize(
    simConfig = cfg, 	
    netParams = netParams)  				# create network object and set cfg and net params
sim.net.createPops()               			# instantiate network populations
sim.net.createCells()              			# instantiate network cells based on defined populations

for i,metype in enumerate(sim.net.cells):

    if 'presyn' not in metype.tags['pop']:

        ii = int(0.000001+(metype.tags['fraction']/(1/cfg.cellNumber_new[metype.tags['pop']])))  

        metype.tags['xnorm'] = cfg.pyr_positions[metype.tags['pop']][ii][0]/cfg.sizeX
        metype.tags['ynorm'] = (3200.929881191896 - cfg.pyr_positions[metype.tags['pop']][ii][1])/cfg.sizeY
        metype.tags['znorm'] = cfg.pyr_positions[metype.tags['pop']][ii][2]/cfg.sizeZ
        metype.tags['x'] = cfg.pyr_positions[metype.tags['pop']][ii][0]
        metype.tags['y'] = 3200.929881191896 - cfg.pyr_positions[metype.tags['pop']][ii][1]
        metype.tags['z'] = cfg.pyr_positions[metype.tags['pop']][ii][2]   

    else: # include external cells as vecStims # not implemented yet
        
        ii = int(metype.tags['cellLabel'])      

        metype.tags['xnorm'] = cfg.pyr_positions[metype.tags['pop']][ii][0]/cfg.sizeX
        metype.tags['ynorm'] = (3200.929881191896 - cfg.pyr_positions[metype.tags['pop']][ii][1])/cfg.sizeY
        metype.tags['znorm'] = cfg.pyr_positions[metype.tags['pop']][ii][2]/cfg.sizeZ
        metype.tags['x'] = cfg.pyr_positions[metype.tags['pop']][ii][0]
        metype.tags['y'] = 3200.929881191896 - cfg.pyr_positions[metype.tags['pop']][ii][1]
        metype.tags['z'] = cfg.pyr_positions[metype.tags['pop']][ii][2]   


sim.net.connectCells()            			# create connections between cells based on params
sim.net.addStims() 							# add network stimulation
sim.setupRecording()              			# setup variables to record for each cell (spikes, V traces, etc)
sim.net.defineCellShapes()


if cfg.addExternalStimulation:
    if tms:
        apply_tms(sim.net, **cfg.tms_params)
    else:
        #Add extracellular stim
        for c,metype in enumerate(sim.net.cells):
            if metype.tags['cellModel'] == 'HH_full':
                secList = [secs for secs in metype.secs.keys() if "pt3d" in metype.secs[secs]['geom']]
                print("\n", metype.tags, "nsec =",len(secList))
                # print(secList)
                v_cell_ext, cell = make_extracellular_stimuli(cfg.acs_params, metype, secList)

sim.runSim()                      			# run parallel Neuron simulation  

if cfg.addExternalStimulation:
    if tms:
        for cell in sim.net.cells:
            cell.t_ext.clear()
            cell.v_ext.clear()
    else:
        for c,metype in enumerate(sim.net.cells):
            if metype.tags['cellModel'] == 'HH_full':
                # print("\n", metype.tags)
                metype.t_ext.clear()
                metype.v_ext.clear()

sim.gatherData()                 			# gather spiking data and cell info from each node
sim.saveData()
sim.analysis.plotData()         			# plot spike raster etc

# sim.analysis.plotTraces(overlay=False, oneFigPer='cell', timeRange = [626.0, 630.0], ylim=[-85,50], figSize=(36,26), fontSize=15, saveFig=True);