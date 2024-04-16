from neuron import h
import os

def loadCell(cellName, cellTemplateName):
    origDir = os.getcwd()
    os.chdir('WeiseEtAl2023/cells/'+cellName+'/')
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")    
    h.load_file("template.hoc")
    # Instantiate the cell from the template
    
    add_synapses=False    
    cell = getattr(h, cellTemplateName)(1 if add_synapses else 0)
    
    os.chdir(origDir)
    
    return cell

def loadCell2(hocName, MorphoName):

    h.load_file('import3d.hoc')
    h.load_file('stdrun.hoc')
    MorphologyPath = '/home/fernando/Documents/SCx_model/O1_data_physiology/morphologies/ascii/'
    gid = 1
    h.load_file('/home/fernando/Documents/SCx_model/O1_data_physiology/emodels_hoc/' + hocName + '.hoc')

    cell = getattr(h,  hocName)(gid,MorphologyPath,MorphoName)

    return cell


def loadCell3(cellName, cellTemplateName):
    origDir = os.getcwd()
    os.chdir('cells/'+cellName+'/')
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
    
    os.chdir(origDir)
    
    return cell


def loadCell4(cellName, cellTemplateName):
    os.chdir('cell_data/'+cellName+'/')
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")    
    h.load_file("template.hoc")
    # Instantiate the cell from the template
    add_synapses=False
    cell = getattr(h, cellTemplateName)(1 if add_synapses else 0)
    
    i=0
    for secs in cell.somatic:
        sec = cell.soma[i]
        listmech = list(cell.soma[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.soma[i].gmax_StochKv_deterministic = 1e-4 * cell.soma[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.soma[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1

    i=0
    for secs in cell.basal:
        sec = cell.dend[i]
        listmech = list(cell.dend[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.dend[i].gmax_StochKv_deterministic = 1e-4 * cell.dend[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.dend[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1

    i=0
    for secs in cell.apical:
        sec = cell.apic[i]
        listmech = list(cell.apic[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.apic[i].gmax_StochKv_deterministic = 1e-4 * cell.apic[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.apic[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1

    i=0
    for secs in cell.axonal:
        sec = cell.axon[i]
        listmech = list(cell.axon[i](0.5))      
        for mech in listmech:
            if str(mech) == 'StochKv':
                sec.insert('StochKv_deterministic')
                cell.axon[i].gmax_StochKv_deterministic = 1e-4 * cell.axon[i].gkbar_StochKv
                # ~ print (sec, mech, i, cell.axon[i].gmax_StochKv_deterministic)
                sec.uninsert('StochKv')
        i=i+1     
    
    print (cell)
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
