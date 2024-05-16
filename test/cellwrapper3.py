from neuron import h
import os
import sys

def loadCell(cellName, cellTemplateName):
    origDir = os.getcwd()
    os.chdir('WeiseEtAl2023/cells/'+cellName+'/')
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")    
    h.load_file("template.hoc")
    # h.xopen('template.hoc')
    # Instantiate the cell from the template
    
    add_synapses=False    
    cell = getattr(h, cellTemplateName)(1 if add_synapses else 0)
    
    os.chdir(origDir)
    
    return cell
    
def loadCell_Net(cellName, cellTemplateName):
    origDir = os.getcwd()
    os.chdir('WeiseEtAl2023/cells/'+cellName+'/')
    h.load_file("stdrun.hoc")
    h.load_file('import3d.hoc')
    try:
        h.xopen("morphology.hoc")
    except:
        pass
    try:
        h.xopen("biophysics.hoc")
    except:
        pass
    try:
        h.xopen("synapses/synapses.hoc")
    except:
        pass
    h.xopen('template.hoc')
    cell = getattr(h, cellTemplateName)(0)
    
    print (cell)
    os.chdir(origDir)
    return cell

# from neuron import h
# import os

# h.load_file("stdrun.hoc")
# h.load_file("import3d.hoc")    
# os.chdir(rootFolder)
# os.chdir('WeiseEtAl2023/')
# # h.load_file("template.hoc")
# h.load_file("nrngui.hoc")
# h.load_file("interpCoordinates.hoc")
# h.load_file("setPointers.hoc")
# h.load_file("calcVe.hoc")
# h.load_file("stimWaveform.hoc")
# h.load_file("cellChooser.hoc")
# h.load_file("setParams.hoc")
# h.load_file("editMorphology.hoc")
# os.chdir('cells/'+cellName+'/')
# h.load_file("createsimulation.hoc")
# # Instantiate the cell from the template
# add_synapses=False
# cell = getattr(h, cellTemplateName)(1 if add_synapses else 0)
# print (cell)
