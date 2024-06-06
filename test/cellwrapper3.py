import os
import sys

# sys.path.append(os.path.join(os.path.dirname(__file__), '../'))

def loadCell(cellName, cellTemplateName):
    from neuron import h
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
    from neuron import h
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
    
def loadCell_L5_TTPC2_cADpyr(cellNumber):

    from cellwrapperTMSlike import L5_TTPC2_cADpyr

    # from neuron import h 

    # h.load_file("stdrun.hoc")

    cell = L5_TTPC2_cADpyr(cellNumber)
    print(cell.loaded)

    cell.load()
    print(cell.loaded)
    
    return cell