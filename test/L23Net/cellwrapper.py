import sys
import os

def loadCell_HL23PYR(cellName):

    templatepath = 'models/NeuronTemplate_HL23PYR.hoc'
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

    cell = getattr(h, 'NeuronTemplate_HL23PYR')(morphpath)
    
    print (cell)

    h.biophys_HL23PYR(cell)

    return cell


def loadCell_HL23VIP(cellName):

    templatepath = 'models/NeuronTemplate_HL23VIP.hoc'
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

    cell = getattr(h, 'NeuronTemplate_HL23VIP')(morphpath)
    
    print (cell)

    h.biophys_HL23VIP(cell)

    return cell


def loadCell_HL23PV(cellName):

    templatepath = 'models/NeuronTemplate_HL23PV.hoc'
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

    cell = getattr(h, 'NeuronTemplate_HL23PV')(morphpath)
    
    print (cell)

    h.biophys_HL23PV(cell)

    return cell


def loadCell_HL23SST(cellName):

    templatepath = 'models/NeuronTemplate_HL23SST.hoc'
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

    cell = getattr(h, 'NeuronTemplate_HL23SST')(morphpath)
    
    print (cell)

    h.biophys_HL23SST(cell)

    return cell
