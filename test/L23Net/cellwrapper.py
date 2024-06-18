import os
import sys

# cellName = 'HL23PYR'

    
def loadCell_Net(cellName):

    templatepath = 'models/NeuronTemplate.hoc'
    biophysics = 'models/biophys_' + cellName + '.hoc'
    morphpath = 'morphologies/' + cellName + '.swc'

    from neuron import h

    h.load_file("stdrun.hoc")
    h.load_file('import3d.hoc')

    h.xopen(biophysics)
        
    try:
       h.xopen(templatepath)
    except:
        pass

    cell = getattr(h, 'NeuronTemplate')(morphpath)
    
    print (cell)

    return cell
    
