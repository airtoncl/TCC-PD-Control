## pd: 0 ou 1 - saudavel ou PD
## t_sim: tempo da simulacao em ms
## freq: frequencia da estimulacao deterministica no TH ou STN. Deve valer 0 se for estocastica
## cortstim: 0 ou 1 - 1 se cortex esta sendo estimulado para analise do PSTH
## str_stim: True ou False - Se ha estimulacao no Estriado para analise de transmissao do sinal senoidal
## ac_str_stim: Valor da oscilacao sobre o nivel DC na estimulacao do Str
## std: desvio padrao do ruido de estimulacao. Deve valer 0 se for deterministico
## stim_type: None, STN ou TH
## f_range: white, low (4, 20) , medium (40, 80), high (80, 300): faixa do ruido colorido 

from netpyne import specs, sim
import random
import numpy as np
import scipy.io as sio
import scipy.signal as sig
from scipy.fftpack import fft, ifft
from matplotlib import pyplot as plt
import pylab
from scipy.signal import welch
from bisect import bisect_left
from scipy.integrate import simpson
#healthy LFP:
#{'StrD1': 0.1200112728246025, 'StrD2': 0.1200112728246025, 'TH': 0.09821650165842044, 'GPi': 0.06291261066904336, 'GPe': 0.015480569050848833, 'CtxRS': 0.14325787235124335, 'CtxFSI': 0.14325787235124335, 'STN': 0.6177382415160736}  

#pd LFP:
#{'StrD1': 0.20870637896851477, 'StrD2': 0.20870637896851477, 'TH': 0.33296972490482285, 'GPi': 0.21134348090184935, 'GPe': 0.2153055978681896, 'CtxRS': 0.22857895992541208, 'CtxFSI': 0.22857895992541208, 'STN': 0.8609998410768475}   

def calc_power_band(lfp_fft, lfp_f, a, b):

    res = list()
    a_ = bisect_left(lfp_f, a)
    b_ = bisect_left(lfp_f, b)
    # Integrate for each channel
    for ch in lfp_fft:
        res.append( simpson( ch[a_:b_], lfp_f[a_:b_] ) )
    return np.array(res)

def full(pd=0, it_num=1, t_sim=2000, freq=0, cortstim=0, str_stim=False, ac_str_stim=0.25,
         std=0, stim_type ='None', f_range = 'high', Input_STN_amp=2e-3):
    

    electrodesPos = [ [5000, 4900, 4000],  # StrD1
                               [5000, 4900, 4000],  # StrD2
                               [1000, 2600, 1800],  # TH
                               [4500, 1200, 1000],  # GPi
                               [4500, 2200, 2700],  # GPe
                               [6500, 7800, 4000],  # CtxRS
                               [6500, 7800, 4000],  # CtxFSI
                               [2000, 1200, 1200] ] # STN

    ###################### Health / Parkinson ################### 
        # RS-> StrD1 connections 
        # GPe-< GPe connections 
        # Str.mod
    #############################################################

        

    #####################################################################################################
    ################################### 1 - Network parameters ##########################################
    #####################################################################################################
        
    netParams = specs.NetParams()  # object of class NetParams to store the network parameters

    netParams.sizeX = 7500  # x-dimension (horizontal length) size in um
    netParams.sizeY = 8800  # y-dimension (vertical height or cortical depth) size in um
    netParams.sizeZ = 5000  # z-dimension (horizontal length) size in um
    #####################################################################################################
    ## 1.I - Population parameters

    netParams.popParams['StrD1'] = {'cellModel': 'StrD1','cellType': 'StrD1', 'numCells': 10,
                                                 'xRange': [4000, 6000],
                                                 'yRange': [3900, 5900],
                                                 'zRange': [3000, 5000]}
    netParams.popParams['StrD2'] = {'cellModel': 'StrD2','cellType': 'StrD2', 'numCells': 10,
                                                 'xRange': [4000, 6000],
                                                 'yRange': [3900, 5900],
                                                 'zRange': [3000, 5000]}

    netParams.popParams['TH'] = {'cellModel': 'TH','cellType': 'Thal', 'numCells': 10,
                                              'xRange': [0, 2000],
                                              'yRange': [1600, 3600],
                                              'zRange': [800, 2800]}

    netParams.popParams['GPi'] = {'cellModel': 'GPi','cellType': 'GPi', 'numCells': 10,
                                               'xRange': [3500, 5500],
                                               'yRange': [200, 2200],
                                               'zRange': [0, 2000]}

    netParams.popParams['GPe'] = {'cellModel': 'GPe','cellType': 'GPe', 'numCells': 10,
                                               'xRange': [3500, 5500],
                                               'yRange': [1200, 3200],
                                               'zRange': [1700, 3700]}

    netParams.popParams['CTX_RS'] = {'cellModel': 'CTX_RS','cellType': 'CTX_RS', 'numCells': 10,
                                                  'xRange': [5500, 7500],
                                                  'yRange': [6800, 8800],
                                                  'zRange': [3000, 5000]}
    netParams.popParams['CTX_FSI'] = {'cellModel': 'CTX_FSI','cellType': 'CTX_FSI', 'numCells': 10,
                                                   'xRange': [5500, 7500],
                                                   'yRange': [6800, 8800],
                                                   'zRange': [3000, 5000]}

    netParams.popParams['STN'] = {'cellType': 'STN', 'numCells': 10,
                                               'xRange': [1000, 3000],
                                               'yRange': [0, 2000],
                                               'zRange': [200, 2200]} #'cellModel': 'STN',


    #####################################################################################################
    ## 1.II - Cell Params

    ## RS
    cellRule = {'conds': {'cellModel': 'CTX_RS', 'cellType': 'CTX_RS'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'pointps': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1, 'cm': 1}
    cellRule['secs']['soma']['pointps']['Izhi'] = {'mod':'Izhi2003b',  'a': 0.02, 'b': 0.2, 'c': -65, 'd': 8, 'f': 5, 'g': 140, 'thresh': 30}
    cellRule['secs']['soma']['vinit'] = -65
    cellRule['secs']['soma']['threshold'] = 30
    netParams.cellParams['CTX_RS'] = cellRule
    

    ## FSI
    cellRule = {'conds': {'cellModel': 'CTX_FSI', 'cellType': 'CTX_FSI'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'pointps': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1, 'cm': 1}
    cellRule['secs']['soma']['pointps']['Izhi'] = {'mod':'Izhi2003b',  'a': 0.1, 'b': 0.2, 'c': -65, 'd': 2, 'f': 5, 'g': 140, 'thresh': 30}
    cellRule['secs']['soma']['vinit'] = -65
    cellRule['secs']['soma']['threshold'] = 30
    netParams.cellParams['CTX_FSI'] = cellRule
    

    ## StrD1
    cellRule = {'conds': {'cellModel': 'StrD1', 'cellType': 'StrD1'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['soma']['mechs']['Str'] = {'gmbar' : (2.6e-3-pd*1.1e-3)}
    cellRule['secs']['soma']['vinit'] = random.gauss(-63.8,5)
    cellRule['secs']['soma']['threshold'] = -10
    netParams.cellParams['StrD1'] = cellRule
    

    ## StrD2
    cellRule = {'conds': {'cellModel': 'StrD2', 'cellType': 'StrD2'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['soma']['mechs']['Str'] = {'gmbar' : (2.6e-3-pd*1.1e-3)}
    cellRule['secs']['soma']['vinit'] = random.gauss(-63.8,5)
    cellRule['secs']['soma']['threshold'] = -10
    netParams.cellParams['StrD2'] = cellRule
    

    ## TH
    cellRule = {'conds': {'cellModel': 'TH', 'cellType': 'Thal'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['soma']['mechs']['thalamus'] = {}
    cellRule['secs']['soma']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['soma']['threshold'] = -10
    netParams.cellParams['TH'] = cellRule
    

    ## GPi
    cellRule = {'conds': {'cellModel': 'GPi', 'cellType': 'GPi'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['soma']['mechs']['GP'] = {}
    cellRule['secs']['soma']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['soma']['threshold'] = -10
    netParams.cellParams['GPi'] = cellRule
    

    ## GPe
    cellRule = {'conds': {'cellModel': 'GPe', 'cellType': 'GPe'},  'secs': {}} 	
    cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['soma']['mechs']['GP'] = {}
    cellRule['secs']['soma']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['soma']['threshold'] = -10
    netParams.cellParams['GPe'] = cellRule
    

    ## STN
    cellRule = {'conds': { 'cellType': 'STN'},  'secs': {}} #'cellModel': 'STN',	
    cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['soma']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['soma']['mechs']['STN'] = {}
    cellRule['secs']['soma']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['soma']['threshold'] = -10
    netParams.cellParams['STN'] = cellRule
    

    ##########################################################################################################
    ## 1.III - Synaptic mechanism parameters

    ## TH
    netParams.synMechParams['Igith'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': -85} #gpi -<th

    ## GPe
    netParams.synMechParams['Insge,ampa'] = {'mod': 'Exp2Syn','tau1': 0.4, 'tau2': 2.5,'e': 0}  #stn -> gpe
    netParams.synMechParams['Insge,nmda'] = {'mod': 'Exp2Syn','tau1': 2, 'tau2': 67,'e': 0}  #stn -> gpe
    netParams.synMechParams['Igege'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': -85}  #gpe -< gpe
    netParams.synMechParams['Istrgpe'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': -85}  #D2 -> gpe

    ## GPi
    netParams.synMechParams['Igegi'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': -85}  #gpe -< gp 
    netParams.synMechParams['Isngi'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': 0}  #stn -> gpi 
    netParams.synMechParams['Istrgpi'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': -85}  #D1 -> gpi 

    ## STN
    netParams.synMechParams['Igesn'] = {'mod': 'Exp2Syn','tau1': 0.4, 'tau2': 7.7,'e': -85}  #gpe -< stn
    netParams.synMechParams['Icosn,ampa'] = {'mod': 'Exp2Syn','tau1': 0.5, 'tau2': 2.49 , 'e': 0}   #ctx -> gpe
    netParams.synMechParams['Icosn,nmda'] = {'mod': 'Exp2Syn','tau1': 2, 'tau2': 90, 'e': 0}   #ctx -> gpe

    ## Str
    netParams.synMechParams['Igabadr'] = {'mod': 'Exp2Syn', 'tau1' : 0.1, 'tau2' : 13,'e': -80} #str -< str
    netParams.synMechParams['Igabaindr'] = {'mod': 'Exp2Syn', 'tau1' : 0.1, 'tau2' : 13,'e': -80} #str -< str
    netParams.synMechParams['Icostr'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': 0} #ctx -> str

    ## CTX
    netParams.synMechParams['Iei'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5,  'e': 0} #rs->fsi 
    netParams.synMechParams['Iie'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5,  'e': -85}  #fsi<-rs 
    netParams.synMechParams['Ithco'] = {'mod': 'Exp2Syn','tau1': 5,'tau2': 5, 'e': 0} #th->rs 


    ##########################################################################################################
    ## 1.IV - Cell connectivity rules

    #########################   TH   ###########################
        ## GPi-> Th connections 
    netParams.connParams['GPi->th'] = {
        'preConds': {'pop': 'GPi'}, 'postConds': {'pop': 'TH'},  # GPi-> th
            'connList':
       [[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]],
        'weight':  0.0336e-3, 		# synaptic weight (conductance) 
            'delay': 5,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Igith'}   		# target synaptic mechanism 


    #########################   GPe   ###########################
        ## STN->GPe connections
    # Two aleatory GPe cells (index i) receive synapse from cells i and i - 1
    aux = random.sample(list(range(0, 10)), 2)
    connList = [[aux[0]-1,aux[0]],[aux[0],aux[0]],[aux[1]-1,aux[1]],[aux[1],aux[1]]]
    weight = []
    for k in range (0,len(connList)):
        weight = weight + [random.uniform(0,0.3)*0.43e-3]
    netParams.connParams['STN->GPe'] = {
        'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'GPe'},  # STN-> GPe
        'connList': connList,               # AMPA
        'weight':  weight, 		# synaptic weight (conductance)
        'delay': 2,					# transmission delay (ms) 
        'loc': 1,					# location of synapse
        'synMech': 'Insge,ampa'}  		# target synaptic mechanism

        ## STN->GPe connections
    # Two aleatory GPe cells (index i) receive synapse from cells i and i - 1
    aux = random.sample(list(range(0, 10)), 2)
    connList = [[aux[0]-1,aux[0]],[aux[0],aux[0]],[aux[1]-1,aux[1]],[aux[1],aux[1]]]
    weight = []
    for k in range (0,len(connList)):
        weight = weight + [random.uniform(0,0.002)*0.43e-3]
    netParams.connParams['STN->GPe2'] = {
        'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'GPe'},  # STN-> GPe
        'connList': connList,                       # NMDA
        'weight':  weight, 		            # synaptic weight (conductance)
        'delay': 2,					# transmission delay (ms) 
        'loc': 1,				    # location of synapse
        'synMech': 'Insge,nmda'}  		# target synaptic mechanism

        ## GPe-< GPe connections
    connList = [[2,1],[3,2],[4,3],[5,4],[6,5],[7,6],[8,7],[9,8],[0,9],[1,0],
                [8,0],[9,1],[0,2],[1,3],[2,4],[3,5],[4,6],[5,7],[6,8],[7,9]]
    weight = []
    for k in range (0,len(connList)):
        weight = weight + [(0.25+0.75*pd)*random.uniform(0,1)*0.3e-3]
    netParams.connParams['GPe->GPe'] = {
        'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'GPe'},  # GPe-< GPe
            'connList': connList,
        'weight':  weight, 		# synaptic weight (conductance)
            'delay': 1,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Igege'}   		# target synaptic mechanism   

        ## StrD2>GPe connections 
    netParams.connParams['StrD2->GPe'] = {
        'preConds': {'pop': 'StrD2'}, 'postConds': {'pop': 'GPe'},  # StrD2-> GPe
            'connList':
        [[0,0],[1,0],[2,0],[3,0],[4,0],[5,0],[6,0],[7,0],[8,0],[9,0],
        [0,1],[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[7,1],[8,1],[9,1],
        [0,2],[1,2],[2,2],[3,2],[4,2],[5,2],[6,2],[7,2],[8,2],[9,2],
        [0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[7,3],[8,3],[9,3],
        [0,4],[1,4],[2,4],[3,4],[4,4],[5,4],[6,4],[7,4],[8,4],[9,4],
        [0,5],[1,5],[2,5],[3,5],[4,5],[5,5],[6,5],[7,5],[8,5],[9,5],
        [0,6],[1,6],[2,6],[3,6],[4,6],[5,6],[6,6],[7,6],[8,6],[9,6],
        [0,7],[1,7],[2,7],[3,7],[4,7],[5,7],[6,7],[7,7],[8,7],[9,7],
        [0,8],[1,8],[2,8],[3,8],[4,8],[5,8],[6,8],[7,8],[8,8],[9,8],
        [0,9],[1,9],[2,9],[3,9],[4,9],[5,9],[6,9],[7,9],[8,9],[9,9],],
        'weight': 0.15e-3, 		# synaptic weight (conductance)
            'delay': 5,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Istrgpe'}   		# target synaptic mechanism


    #########################   GPi   ###########################        
        ## STN-> GPi connections
    # Five aleatory GPi cells (index i) receive synapse from cells i and i - 1
    aux = random.sample(list(range(0, 10)), 5)
    ## Parte auxiliar para o PSTH:
    gsngi = np.zeros(10)
    for k in range(0,10):
        if (k == aux[0] or k == aux[1] or k == aux[2] or k == aux[3] or k == aux[4]):
            gsngi[k] = 1
        else:
            gsngi[k] = 0
    connList = [[aux[0]-1,aux[0]],[aux[0],aux[0]],[aux[1]-1,aux[1]],[aux[1],aux[1]],[aux[2]-1,aux[2]],[aux[2],aux[2]],
                [aux[3]-1,aux[3]],[aux[3],aux[3]],[aux[4]-1,aux[4]],[aux[4],aux[4]]]
    netParams.connParams['STN->GPi'] = {
        'preConds': {'pop': 'STN'}, 'postConds': {'pop': 'GPi'}, 
        'connList': connList,
        'weight':  0.0645e-3, 		# synaptic weight (conductance) 
        'delay': 1.5,			    # transmission delay (ms) 
        'loc': 1,					# location of synapse
        'synMech': 'Isngi'}   		# target synaptic mechanism

        ## GPe-< GPi connections 
    netParams.connParams['GPe->GPi'] = {
        'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'GPi'},  
            'connList':
       [[9,1],[0,2],[1,3],[2,4],[3,5],[4,6],[5,7],[6,8],[7,9],[8,0],
        [1,0],[2,1],[3,2],[4,3],[5,4],[6,5],[7,6],[8,7],[9,8],[0,9]],
        'weight':  0.15e-3, 		# synaptic weight (conductance) 
            'delay': 3,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Igegi'}   		# target synaptic mechanism

        ## StrD1>GPi connections 
    netParams.connParams['StrD1->GPe'] = {
        'preConds': {'pop': 'StrD1'}, 'postConds': {'pop': 'GPi'},  # StrD1-> GPi
            'connList':
        [[0,0],[1,0],[2,0],[3,0],[4,0],[5,0],[6,0],[7,0],[8,0],[9,0],
        [0,1],[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[7,1],[8,1],[9,1],
        [0,2],[1,2],[2,2],[3,2],[4,2],[5,2],[6,2],[7,2],[8,2],[9,2],
        [0,3],[1,3],[2,3],[3,3],[4,3],[5,3],[6,3],[7,3],[8,3],[9,3],
        [0,4],[1,4],[2,4],[3,4],[4,4],[5,4],[6,4],[7,4],[8,4],[9,4],
        [0,5],[1,5],[2,5],[3,5],[4,5],[5,5],[6,5],[7,5],[8,5],[9,5],
        [0,6],[1,6],[2,6],[3,6],[4,6],[5,6],[6,6],[7,6],[8,6],[9,6],
        [0,7],[1,7],[2,7],[3,7],[4,7],[5,7],[6,7],[7,7],[8,7],[9,7],
        [0,8],[1,8],[2,8],[3,8],[4,8],[5,8],[6,8],[7,8],[8,8],[9,8],
        [0,9],[1,9],[2,9],[3,9],[4,9],[5,9],[6,9],[7,9],[8,9],[9,9],],
        'weight': 0.15e-3, 		# synaptic weight (conductance) 
            'delay': 4,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Istrgpi'}   		# target synaptic mechanism
        
        
    #########################   STN   ###########################
        ## GPe-> STN connections 
    netParams.connParams['GPe->STN'] = {
        'preConds': {'pop': 'GPe'}, 'postConds': {'pop': 'STN'},  # GPe-< STN
            'connList':
       [[2,1],[3,2],[4,3],[5,4],[6,5],[7,6],[8,7],[9,8],[0,9],[1,0],
        [0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]],
        'weight':  0.15e-3, 		# synaptic weight (conductance)
            'delay': 4,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Igesn'}   		# target synaptic mechanism

        ## CTX-> STN connections
    connList = [[2,1],[3,2],[4,3],[5,4],[6,5],[7,6],[8,7],[9,8],[0,9],[1,0],
                [0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]]
    weight = []
    for k in range (0,len(connList)):
        weight = weight + [random.uniform(0,0.3)*0.43e-3]
    netParams.connParams['CTX->STN'] = {
        'preConds': {'pop': 'CTX_RS'}, 'postConds': {'pop': 'STN'},  # CTX-> STN
            'connList': connList,
        'weight':  weight, 		# synaptic weight (conductance)
            'delay': 5.9,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Icosn,ampa'}   		# target synaptic mechanism   

        ## CTX-> STN connections 
    connList = [[2,1],[3,2],[4,3],[5,4],[6,5],[7,6],[8,7],[9,8],[0,9],[1,0],
                [0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]]
    weight = []
    for k in range (0,len(connList)):
        weight = weight + [random.uniform(0,0.003)*0.43e-3]
    netParams.connParams['CTX->STN2'] = {
        'preConds': {'pop': 'CTX_RS'}, 'postConds': {'pop': 'STN'},  # CTX-> STN
            'connList': connList,
        'weight':  weight, 		# synaptic weight (conductance)
            'delay': 5.9,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Icosn,nmda'}   		# target synaptic mechanism   


    #########################   Str   ###########################
        ## StrD2-< StrD2 connections
    # Each StrD2 cell receive synapse from 4 aleatory StrD2 cell (except from itself)
    connList = []
    for i in range (0, 10):
        #cells = range(0, 10)
        cells = list(range(0, 10))
        cells.remove(i)
        aux = []
        for k in range(0,4):
            aux = aux+ [random.choice(cells)]
        connList = connList + [[aux[0],i], [aux[1], i], [aux[2],i], [aux[3],i]]
    netParams.connParams['StrD2-> StrD2'] = {
        'preConds': {'pop': 'StrD2'}, 'postConds': {'pop': 'StrD2'},  #StrD2-< StrD2 
        'connList': connList,
        'weight':  0.1/4*0.5e-3, 		# synaptic weight (conductance) -> mudar essa maluquisse
        'delay': 0,					# transmission delay (ms) 
        'loc': 1,					# location of synapse
        'synMech': 'Igabaindr'}   		# target synaptic mechanism

        ## StrD1-< StrD1 connections
    # Each StrD1 cell receive synapse from 3 aleatory StrD1 cell (except from itself)
    connList = []
    for i in range (0, 10):
        #cells = range(0, 10)
        cells = list(range(0, 10))
        cells.remove(i)
        aux = []
        for k in range(0,3):
            aux = aux+ [random.choice(cells)]
        connList = connList + [[aux[0],i], [aux[1], i], [aux[2],i]]
    netParams.connParams['StrD1-> StrD1'] = {
        'preConds': {'pop': 'StrD1'}, 'postConds': {'pop': 'StrD1'},  #StrD1-< StrD1 
        'connList': connList,
        'weight':  0.1/3*0.5e-3, 		# synaptic weight (conductance) -> mudar aqui tb
        'delay': 0,					# transmission delay (ms) 
        'loc': 1,					# location of synapse
        'synMech': 'Igabadr'}   		# target synaptic mechanism
        
         ## RS-> StrD1 connections 
    netParams.connParams['RS-> StrD1'] = {
        'preConds': {'pop': 'CTX_RS'}, 'postConds': {'pop': 'StrD1'},  # RS-> StrD1
            'connList':
        [[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]],
        'weight':  (0.07-0.044*pd)*0.43e-3, 		# synaptic weight (conductance) 
            'delay': 5.1,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Icostr'}   		# target synaptic mechanism

         ## RS-> StrD2 connections 
    netParams.connParams['RS-> StrD2'] = {
        'preConds': {'pop': 'CTX_RS'}, 'postConds': {'pop': 'StrD2'},  # RS-> StrD2 
            'connList':
        [[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]],
        'weight':  0.07*0.43e-3, 		# synaptic weight (conductance) 
            'delay': 5.1,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Icostr'}   		# target synaptic mechanism

     
    #########################   CTX   ###########################    
        ## RS -> FSI connections
    # Each FSI cell receive synapse from 4 aleatory RS cells
    connList = []
    for i in range (0, 10):
        aux = []
        for k in range(0,4):
            aux = aux+ [random.choice(list(range(0,10)))]
        connList = connList + [[aux[0],i], [aux[1], i], [aux[2],i], [aux[3],i]]
    netParams.connParams['ctx_rs->ctx_fsi'] = {
        'preConds': {'pop': 'CTX_RS'}, 'postConds': {'pop': 'CTX_FSI'},  #  ctx_rs -> ctx_fsi
            'connList': connList,
        'weight': 0.043e-3, 		# synaptic weight (conductance) 
            'delay': 1,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Iei'}   		# target synaptic mechanism

        ## FSI -> RS connections
    # Each RS cell receive synapse from 4 aleatory FSI cells
    connList = []
    for i in range (0, 10):
        aux = []
        for k in range(0,4):
            aux = aux+ [random.choice(list(range(0,10)))]
        connList = connList + [[aux[0],i], [aux[1], i], [aux[2],i], [aux[3],i]]
    netParams.connParams['ctx_fsi->ctx_rs'] = {
        'preConds': {'pop': 'CTX_FSI'}, 'postConds': {'pop': 'CTX_RS'},  #  ctx_fsi -< ctx_rs
            'connList': connList,
        'weight': 0.083e-3, 		# synaptic weight (conductance)
            'delay': 1,					# transmission delay (ms) 
            'loc': 1,					# location of synapse
        'synMech': 'Iie'}   		# target synaptic mechanism

        ## Th -> RS connections 
    netParams.connParams['th->ctx_rs'] = {
        'preConds': {'pop': 'TH'}, 'postConds': {'pop': 'CTX_RS'},  #  th -> ctx_rs
            'connList':
       [[0,0],[1,1],[2,2],[3,3],[4,4],[5,5],[6,6],[7,7],[8,8],[9,9]],
        'weight':  0.0645e-3, 		# synaptic weight (conductance)
            'delay': 5,					# transmission delay (ms)
            'loc': 1,					# location of synapse
        'synMech': 'Ithco'}   		# target synaptic mechanism
        

    #####################################################################################################
    ## 1.V - Stimulation Source Parameters   
    dbs_times = []
    pulse_amp = 2.0e-3  # Amplitude dos pulsos de DBS (pode ajustar)
    pulse_dur = 1.0     # Duração de cada pulso (ms)

    if stim_type == 'STN' and freq > 0:
        period = 1000.0 / freq  # ms entre pulsos
        dbs_times = np.arange(0, t_sim, period).tolist()


    # FS receve a constante 3 density current or 1 during cortical stimulation
    netParams.stimSourceParams['Input_FS'] = {'type': 'IClamp',
                                                    'delay': 0,
                                                    'dur': t_sim,
                                                    'amp': Input_STN_amp} #0
    netParams.stimTargetParams['Input_FS->FS'] = {'source': 'Input_FS',
                                                        'conds': {'pop': 'CTX_FSI'},
                                                        'sec': 'soma',
                                                        'loc': 0}

    # RS receve a constante 3 density current or 1 during cortical stimulation
    netParams.stimSourceParams['Input_RS'] = {'type': 'IClamp',
                                                    'delay': 0,
                                                    'dur': t_sim,
                                                    'amp': 0}
    netParams.stimTargetParams['Input_RS->RS'] = {'source': 'Input_RS',
                                                        'conds': {'pop': 'CTX_RS'},
                                                        'sec': 'soma',
                                                        'loc': 0}

        ## GPe receive a constante 3 density current or 1 during cortical stimulation
    netParams.stimSourceParams['Input_GPe'] = {'type': 'IClamp', 'delay': 0, 'dur': t_sim,
                                               'amp': (3e-3 - cortstim*(not pd)*2e-3)} #Stimulus Amplitude (in nA)
    netParams.stimTargetParams['Input_GPe->GPe'] = {'source': 'Input_GPe',  'conds': {'cellType':'GPe'},'sec':'GPe',
                                                    'loc': 0}

        ## GPi receive a constante 3 density current
    netParams.stimSourceParams['Input_GPi'] = {'type': 'IClamp', 'delay': 0, 'dur': t_sim, 'amp': 3e-3}
    netParams.stimTargetParams['Input_GPi->GPi'] = {'source': 'Input_GPi',  'conds': {'cellType':'GPi'},'sec':'GPi',
                                                    'loc': 0}     
        ## Thalamus receive a constante 1.2 density current
    netParams.stimSourceParams['Input_th'] = {'type': 'IClamp', 'delay': 0, 'dur': t_sim,
                                              'amp': 1.2e-3} 
    netParams.stimTargetParams['Input_th->TH'] = {'source': 'Input_th',  'conds': {'cellType':'Thal'},'sec':'TH',
                                                  'loc': 0}
        ## Striatum receive a constante 1.2 density current
    if(str_stim == True):
        netParams.stimSourceParams['Input_Str'] = {'type': 'IClamp', 'delay': 0, 'dur': t_sim,
                                                   'amp': 2.4e-3} 
        netParams.stimTargetParams['Input_Str->Str'] = {'source': 'Input_Str',  'conds': {'cellType':['StrD1', 'StrD2']},
                                                        'loc': 0}

        ## STN receve a constante 3 density current or 1 during cortical stimulation
    if stim_type == 'STN' and freq > 0:
        # DBS com pulsos periódicos via VecStim
        netParams.stimSourceParams['DBS_STN'] = {
            'type': 'NetStim',
            'start': 0,
            'interval': 1,
            'number': len(dbs_times),
            'noise': 0,
            'spikeTimes': dbs_times
        }
        
        netParams.stimTargetParams['DBS_STN->STN'] = {
            'source': 'DBS_STN',
            'conds': {'pop': 'STN'},
            'sec': 'soma',
            'loc': 0,
            'weight': pulse_amp,   # equivalente à amplitude da corrente
            'synMech': 'dbs_syn'
        }

        # Define o mecanismo sináptico como impulso breve
        netParams.synMechParams['dbs_syn'] = {
            'mod': 'Exp2Syn',
            'tau1': 0.1,
            'tau2': pulse_dur,
            'e': 0  # excitatório
        }

    else:
        # Caso não haja DBS, use corrente direta padrão
        netParams.stimSourceParams['Input_STN'] = {
            'type': 'IClamp',
            'delay': 0,
            'dur': t_sim,
            'amp': Input_STN_amp
        }
        netParams.stimTargetParams['Input_STN->STN'] = {
            'source': 'Input_STN',
            'conds': {'pop': 'STN'},
            'sec': 'soma',
            'loc': 0
        }


    # dStr receve a constante 3 density current
    netParams.stimSourceParams['Input_StrD1'] = {'type': 'IClamp',
                                                        'delay': 0,
                                                        'dur': t_sim,
                                                        'amp': 0}
    netParams.stimTargetParams['Input_StrD1->StrD1'] = {'source': 'Input_StrD1',
                                                                'conds': {'pop': 'StrD1'},
                                                                'sec': 'soma',
                                                                'loc': 0}

    # iStr receve a constante 3 density current
    netParams.stimSourceParams['Input_StrD2'] = {'type': 'IClamp',
                                                        'delay': 0, 'dur': t_sim,
                                                        'amp': 0}
    netParams.stimTargetParams['Input_StrD2->StrD2'] = {'source': 'Input_StrD2',
                                                                'conds': {'pop': 'StrD2'},
                                                                'sec': 'soma',
                                                                'loc': 0}


    #####################################################################################################
    ################################### 2 - Simulation parameters #######################################
    #####################################################################################################


    ## Create simulation
    dt = 0.1
    simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration    
    simConfig.duration = t_sim+dt          # Duration of the simulation, in ms
    simConfig.dt = dt                # Internal integration timestep to use
    simConfig.verbose = False           # Show detailed messages

    simConfig.recordStep = 1 # Step size in ms to save data (eg. V traces, LFP, etc)
    simConfig.recordCells = ['allCells']
    simConfig.recordSpikesGids = True
    simConfig.printPopAvgRates = True

    #simConfig.saveJson = True
    #simConfig.filename = 'sim_output'  # Nome base dos arquivos

    simConfig.recordLFP = [ [5000, 4900, 4000],  # StrD1
                               [5000, 4900, 4000],  # StrD2
                               [1000, 2600, 1800],  # TH
                               [4500, 1200, 1000],  # GPi
                               [4500, 2200, 2700],  # GPe
                               [6500, 7800, 4000],  # CtxRS
                               [6500, 7800, 4000],  # CtxFSI
                               [2000, 1200, 1200] ] # STN
    simConfig.saveLFPCells = True

    sim.createSimulateAnalyze(netParams=netParams, simConfig=simConfig)
        

    nelec = 1

    lfp = sim.allSimData['LFP']
    # [ f, t ]
    lfp = np.transpose(lfp, [1, 0])

    # calculate LFP using Welch method
    lfp_f, lfp_dimensions = sig.welch( lfp[0], 1000, nperseg=1024, detrend=False )
    lfp_fft = np.zeros(( len(electrodesPos)//nelec, lfp_dimensions.shape[0] ))
    for i in range( 0, lfp.shape[0], nelec ):
        reg_fft = list()
        for j in range( nelec ):
            reg_fft.append( sig.welch( lfp[i+j], 1000, nperseg=1024, detrend=False ) )
        lfp_f, lfp_fft[i//nelec, :] = np.mean( reg_fft, axis=0 )

    denominator = calc_power_band( lfp_fft, lfp_f, 0, 51 )
    alpha = calc_power_band( lfp_fft, lfp_f, 8, 14 )
    beta  = calc_power_band( lfp_fft, lfp_f, 13, 31 )
    betaH = calc_power_band( lfp_fft, lfp_f, 30, 51 )
    coefs = ( alpha + beta + betaH ) / denominator

    lfp_coef = dict()
    area_names = ["StrD1", "StrD2", "TH", "GPi", "GPe", "CtxRS", "CtxFSI", "STN"]
    for i in range( len(coefs) ):
        lfp_coef[ area_names[i] ] = ( coefs[i] )
    
    #print("#####################################################################################################")
    #print(lfp_coef)
    #print("#####################################################################################################")


    #####################################################################################################
    ################################### 3 - Recording Data ##############################################
    #####################################################################################################


    #####################################################################################################
    pylab.show()
    
    lfp = np.array(sim.allSimData['LFP'])  # shape: [num_eletrodos, num_tsteps]
    lfp = np.transpose(lfp, [1, 0])
    fs = 1.0 / simConfig.dt  # ex: 20 kHz → 0.05 ms -> fs = 20 kHz

    # Plot do sinal bruto
    #plt.figure(figsize=(10, 3))
    #plt.plot(lfp[0])
    #plt.title("LFP - Eletrodo 0")
    #plt.xlabel("Tempo (passos)")
    #plt.ylabel("Potencial (uV)")
    #plt.grid(True)
    #plt.tight_layout()
    #plt.savefig('teste1_PD.png')

    return lfp_coef, lfp

    # Espectro de potência via Welch
    #f, Pxx = welch(lfp[0], fs=fs, nperseg=1024)
    #plt.figure(figsize=(8, 4))
    #plt.semilogy(f, Pxx)
    #plt.xlim([0, 100])
    #plt.title('Espectro de Potência - LFP (Eletrodo 0)')
    #plt.xlabel('Frequência (Hz)')
    #plt.ylabel('Potência (uV^2/Hz)')
    #plt.grid(True)
    #plt.tight_layout()
    #plt.savefig('teste2_PD.png')

    # Potência Beta (13–30 Hz)
    #beta_power = np.trapz(Pxx[(f >= 13) & (f <= 30)], f[(f >= 13) & (f <= 30)])
    #print(f"Potência beta: {beta_power:.4f} uV^2")

    #print(sim.allSimData.keys())  # deve conter 'LFP'
    #print(len(sim.allSimData['LFP']))