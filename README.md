# NetPyNE implementation of Extracellular stimulation
Extracellular stimulation in NetPyNE


## Description
This code reproduces the simulations of Fig 4 (poster_TMS_SfN_2024.pdf)

Fernando S Borges, Jacob Tajchman, Tarek Khashan, Eugenio Urdapilleta, Emiliano Santarnecchi, Salvador Dura-Bernal. Towards a detailed mechanistic model of human cortical microcircuits that accurately predicts the cellular- and circuit-level effects of TMS. Society for Neuroscience, 2024.

The "v1_batch5" reproduces "data/v1_batch5_TMS_Human.png"

### Branches
1. `main`: files needed to run the code (4.5 GB)

## Setup and execution
Requires NEURON with Python and MPI support. 

### NEURON libraries 
1. From /sim run `nrnivmodl mod`. This should create a directory called x86_64. 
2. In cfg.py make sure cfg.coreneuron = False
3. To run type: `python batch.py` or `mpiexec -n [num_proc] nrniv -python -mpi init.py`

### CoreNEURON libraries
1. From /sim run `nrnivmodl -coreneuron mod`. This should create a directory called x86_64. 
2. In cfg.py make sure cfg.coreneuron = True
3. To run type: `python batch.py` or `mpirun -n [num_proc] ./x86_64/special -mpi -python init.py`


## Overview of file structure:

* `/sim/init.py`: Main executable; calls functions from other modules. Sets what parameter file to use.

* `/sim/netParams.py`: Network parameters

* `/sim/cfg.py`: Simulation configuration

* `/sim/batch.py`: Run multiple simulations

* `/sim/cells`: source files for the different cell types used in the model; these will be imported into netpyne

* `/sim/mod`: NMODL files containing the ionic channel and synaptic mechanisms used in the model 

* `/data`: where the model and simulation data is stored 

* `/test`: tests, comparations, and validations



## Extra description see:
 
Fernando da Silva Borges,  Joao V.S. Moreira,  Lavinia M. Takarabe,  William W. Lytton,  Salvador Dura-Bernal. **Large-scale biophysically detailed model of somatosensory thalamocortical circuits in NetPyNE**. Frontiers in Neuroinformatics. https://doi.org/10.1101/2022.02.03.479029

https://github.com/suny-downstate-medical-center/S1_Thal_NetPyNE_Frontiers_2022



For further information please contact: fernandodasilvaborges@gmail.com


