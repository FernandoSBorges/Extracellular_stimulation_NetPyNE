"""
batch.py 

Batch simulation for S1 model using NetPyNE

Contributors: salvadordura@gmail.com, fernandodasilvaborges@gmail.com
"""
from netpyne.batch import Batch
from netpyne import specs
import numpy as np

# ----------------------------------------------------------------------------------------------
# Custom
# ----------------------------------------------------------------------------------------------
def custom():
    params = specs.ODict()
    
    params[('seeds', 'stim')] =  [1000]

    # params[('fracmorphoradius')] = [1.0/2.0]

    b = Batch(params=params, netParamsFile='netParams.py', cfgFile='cfg.py')

    return b

# ----------------------------------------------------------------------------------------------
# Run configurations
# ----------------------------------------------------------------------------------------------
def setRunCfg(b, type='mpi_bulletin'):
    if type=='mpi_bulletin' or type=='mpi':
        b.runCfg = {'type': 'mpi_bulletin', 
            'script': 'init.py', 
            'skip': True}

    elif type=='mpi_direct':
        b.runCfg = {'type': 'mpi_direct',
            'cores': 12,
            'script': 'init.py',
            'mpiCommand': 'mpiexec', # --use-hwthread-cpus
            'skip': True}

    elif type=='mpi_direct2':
        b.runCfg = {'type': 'mpi_direct',
            'mpiCommand': 'mpirun -n 80 ./x86_64/special -mpi -python init.py', # --use-hwthread-cpus
            'skip': True}
            
-----------------------------------------------------------------------------
# Main code
# ----------------------------------------------------------------------------------------------
if __name__ == '__main__': 
    b = custom() #

    b.batchLabel = 'v1_batch2'  
    b.saveFolder = '../data/'+b.batchLabel
    b.method = 'grid'
    setRunCfg(b, 'mpi_direct')
    b.run() # run batch
