## pd: 0 ou 1 - saudavel ou PD
## t_sim: tempo da simulacao em ms
## freq: frequencia da estimulacao deterministica no TH ou STN. Deve valer 0 se for estocastica
## cortstim: 0 ou 1 - 1 se cortex esta sendo estimulado para analise do PSTH
## str_stim: True ou False - Se ha estimulacao no Estriado para analise de transmissao do sinal senoidal
## ac_str_stim: Valor da oscilacao sobre o nivel DC na estimulacao do Str
## std: desvio padrao do ruido de estimulacao. Deve valer 0 se for deterministico
## stim_type: None, STN ou TH
## f_range: white, low (4, 20) , medium (40, 80), high (80, 300): faixa do ruido colorido 

def calc_power_band(lfp_fft, lfp_f, a, b):
    from bisect import bisect_left
    from scipy.integrate import simpson
    import numpy as np

    res = list()
    a_ = bisect_left(lfp_f, a)
    b_ = bisect_left(lfp_f, b)
    # Integrate for each channel
    for ch in lfp_fft:
        res.append( simpson( ch[a_:b_], lfp_f[a_:b_] ) )
    return np.array(res)

def extractLFP_SP(self):
        import numpy as np
        from scipy import signal
        nelec = 1

        lfp = sim.allSimData['LFP']
        # [ f, t ]
        lfp = np.transpose(lfp, [1, 0])

        # calculate LFP using Welch method
        lfp_f, lfp_dimensions = signal.welch( lfp[0], 1000, nperseg=1024, detrend=False )
        lfp_fft = np.zeros(( len(electrodesPos)//nelec, lfp_dimensions.shape[0] ))
        for i in range( 0, lfp.shape[0], nelec ):
            reg_fft = list()
            for j in range( nelec ):
                reg_fft.append( signal.welch( lfp[i+j], 1000, nperseg=1024, detrend=False ) )
            lfp_f, lfp_fft[i//nelec, :] = np.mean( reg_fft, axis=0 )
        return lfp_f, lfp_fft


def extractLFP_raw(self):
    lfp = sim.allSimData['LFP']
    # [ f, t ]
    lfp = np.transpose(lfp, [1, 0])
    return lfp

def full(pd=0, it_num=1, t_sim=1000, freq=0, cortstim=0, str_stim=False, ac_str_stim=0.25,
         std=0, stim_type ='None', f_range = 'high'):
    from netpyne import specs, sim
    import random
    import numpy as np
    import scipy.io as sio
    import scipy.signal as sig
    from scipy.fftpack import fft, ifft
    import math
    from matplotlib import pyplot as plt
    import os

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
    cellRule['secs']['rs'] = {'geom': {}, 'pointps': {}}  								
    cellRule['secs']['rs']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1, 'cm': 1}
    cellRule['secs']['rs']['pointps']['Izhi'] = {'mod':'Izhi2003b',  'a': 0.02, 'b': 0.2, 'c': -65, 'd': 8, 'f': 5, 'g': 140, 'thresh': 30}
    cellRule['secs']['rs']['vinit'] = -65
    cellRule['secs']['rs']['threshold'] = 30
    netParams.cellParams['CTX_RS'] = cellRule
    

    ## FSI
    cellRule = {'conds': {'cellModel': 'CTX_FSI', 'cellType': 'CTX_FSI'},  'secs': {}} 	
    cellRule['secs']['fsi'] = {'geom': {}, 'pointps': {}}  								
    cellRule['secs']['fsi']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1, 'cm': 1}
    cellRule['secs']['fsi']['pointps']['Izhi'] = {'mod':'Izhi2003b',  'a': 0.1, 'b': 0.2, 'c': -65, 'd': 2, 'f': 5, 'g': 140, 'thresh': 30}
    cellRule['secs']['fsi']['vinit'] = -65
    cellRule['secs']['fsi']['threshold'] = 30
    netParams.cellParams['CTX_FSI'] = cellRule
    

    ## StrD1
    cellRule = {'conds': {'cellModel': 'StrD1', 'cellType': 'StrD1'},  'secs': {}} 	
    cellRule['secs']['StrD1'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['StrD1']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['StrD1']['mechs']['Str'] = {'gmbar' : (2.6e-3-pd*1.1e-3)}
    cellRule['secs']['StrD1']['vinit'] = random.gauss(-63.8,5)
    cellRule['secs']['StrD1']['threshold'] = -10
    netParams.cellParams['StrD1'] = cellRule
    

    ## StrD2
    cellRule = {'conds': {'cellModel': 'StrD2', 'cellType': 'StrD2'},  'secs': {}} 	
    cellRule['secs']['StrD2'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['StrD2']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['StrD2']['mechs']['Str'] = {'gmbar' : (2.6e-3-pd*1.1e-3)}
    cellRule['secs']['StrD2']['vinit'] = random.gauss(-63.8,5)
    cellRule['secs']['StrD2']['threshold'] = -10
    netParams.cellParams['StrD2'] = cellRule
    

    ## TH
    cellRule = {'conds': {'cellModel': 'TH', 'cellType': 'Thal'},  'secs': {}} 	
    cellRule['secs']['th'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['th']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['th']['mechs']['thalamus'] = {}
    cellRule['secs']['th']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['th']['threshold'] = -10
    netParams.cellParams['TH'] = cellRule
    

    ## GPi
    cellRule = {'conds': {'cellModel': 'GPi', 'cellType': 'GPi'},  'secs': {}} 	
    cellRule['secs']['GPi'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['GPi']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['GPi']['mechs']['GP'] = {}
    cellRule['secs']['GPi']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['GPi']['threshold'] = -10
    netParams.cellParams['GPi'] = cellRule
    

    ## GPe
    cellRule = {'conds': {'cellModel': 'GPe', 'cellType': 'GPe'},  'secs': {}} 	
    cellRule['secs']['GPe'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['GPe']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['GPe']['mechs']['GP'] = {}
    cellRule['secs']['GPe']['vinit'] = random.gauss(-62,5)
    cellRule['secs']['GPe']['threshold'] = -10
    netParams.cellParams['GPe'] = cellRule
    

    ## STN
    cellRule = {'conds': { 'cellType': 'STN'},  'secs': {}} #'cellModel': 'STN',	
    cellRule['secs']['STN'] = {'geom': {}, 'mechs': {}}  								
    cellRule['secs']['STN']['geom'] = {'diam': 5.642, 'L': 5.642, 'Ra': 1, 'nseg': 1}
    cellRule['secs']['STN']['mechs']['STN'] = {}
    cellRule['secs']['STN']['vinit'] = random.gauss(-62,5)
    netParams.cellParams['STN'] = cellRule
    cellRule['secs']['STN']['threshold'] = -10
    


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


    #####################################################################################################
    ################################### 2 - Simulation parameters #######################################
    #####################################################################################################


    ## Create simulation
    dt = 0.05
    simConfig = specs.SimConfig()       # object of class SimConfig to store simulation configuration    
    simConfig.duration = t_sim+dt          # Duration of the simulation, in ms
    simConfig.dt = dt                # Internal integration timestep to use
    simConfig.verbose = False           # Show detailed messages

    ##simConfig.recordTraces = {'vsn':{'sec':'STN','loc':0.5,'var':'v'}}
    simConfig.recordStep = dt # Step size in ms to save data (eg. V traces, LFP, etc)
    simConfig.recordCells = ['allCells']
    simConfig.recordSpikesGids = True
    simConfig.printPopAvgRates = True

    simConfig.saveJson = True
    simConfig.filename = 'sim_output'  # Nome base dos arquivos

    simConfig.analysis['plotTraces'] = {'include': [('STN',[1])]}
    simConfig.analysis['plotRaster'] ={'include': ['CTX_RS', 'CTX_FSI', 'TH', 'GPi', 'GPe', 'STN', 'StrD2', 'StrD1']}
    
    # --- Definições da Simulação ---
    # Registrar Potencial de Membrana (Traços)
    #simConfig.analysis['plotTraces'] = {
    #    'include': [('STN', [0])],  # Célula 0 da população STN
    #    'saveFig': True,
    #    'saveFigPath': './results',
    #    'saveFigFormats': ['png']
    #}
#
    ## Registrar Raster Plot
    #simConfig.analysis['plotRaster'] = {
    #    'include': ['CTX_RS', 'CTX_FSI', 'TH', 'GPi', 'GPe', 'STN', 'StrD2', 'StrD1'],
    #    'saveFig': True,
    #    'saveFigPath': './results',
    #    'saveFigFormats': ['png']
    #}
#
    ## Registrar LFP
    #simConfig.recordLFP = [[0, 0, 0]]  # Um eletrodo virtual no centro da rede
#
    ## Garantir que o diretório de resultados exista
    #if not os.path.exists('./results'):
    #    os.makedirs('./results')
    #simConfig.recordLFP = [[0, 0, 0], [100, 0, 0], [200, 0, 0]]  # 3 eletrodos simulados em posições distintas
    simConfig.recordLFP = [ [5000, 4900, 4000],  # StrD1
                               [5000, 4900, 4000],  # StrD2
                               [1000, 2600, 1800],  # TH
                               [4500, 1200, 1000],  # GPi
                               [4500, 2200, 2700],  # GPe
                               [6500, 7800, 4000],  # CtxRS
                               [6500, 7800, 4000],  # CtxFSI
                               [2000, 1200, 1200] ] # STN
    simConfig.saveLFPCells = True
    #simConfig.analysis['plotRaster'] = True
    #simConfig.analysis['plotLFP'] = {'electrodes': ['all'],
    #                                  'includeAxon': False,
    #                                  'timeRange': [0, 2000],
    #                                  'plots': ['timeSeries', 'locations', 'PSD'],
    #                                  'plots': ['locations'],
    #                                  'showFig': True}

    # Atribui um yRange diferente para cada núcleo
    #pop_positions = {
    #    'CTX': [0, 50],
    #    'TH': [100, 150],
    #    'STR': [200, 250],
    #    'STN': [300, 350],
    #    'GPe': [400, 450],
    #    'GPi': [500, 550]
    #}
#
    #for popName, yRange in pop_positions.items():
    #    if popName in netParams.popParams:
    #        netParams.popParams[popName]['yRange'] = yRange

    pops, cells, conns, stims, simData, *rest = sim.create(netParams=netParams, simConfig=simConfig, output=True)
    
    
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
    
    print("#####################################################################################################")
    print(lfp_coef)
    print("#####################################################################################################")

    #####################################################################################################
    ################################### 3 - Injected Currents ###########################################
    #####################################################################################################
    
    t = np.arange(0, t_sim+dt, dt)
    i_est_stn = np.zeros(len(t))
    i_est_th = np.zeros(len(t))
    i_est_str = np.zeros(len(t))
    i_ctx = np.zeros(len(t))
    
    ## i_est_stn is the an arbitrary current stimulation in the STN cells
    if (stim_type == 'STN'):
        if (std != 0) and (freq == 0): ## Stochastic
            if f_range == 'white': # White
                i_est_stn = np.random.normal(0, std, len(t))
            else: # Colored
                i_est_stn = np.random.normal(0, std, len(t)-1);
                filtro = np.zeros(len(t)-1);
                if f_range == 'low':
                    faixa = (4,20);
                elif f_range == 'medium':
                    faixa = (40,80);
                elif f_range == 'high':
                    faixa = (80,300);
                filtro[(faixa[0]*100):(faixa[1]*100+1)] = np.ones(faixa[1]*100-faixa[0]*100+1);
                filtro[(-faixa[1]*100):(-faixa[0]*100+1)] = np.ones(faixa[1]*100-faixa[0]*100+1);
                s = fft(i_est_stn);
                s = np.multiply(s,filtro);
                ruido_colorido = ifft(s);
                i_est_stn = np.real(ruido_colorido);
                i_est_stn = np.concatenate((i_est_stn, [0]))    
        elif (std == 0) and (freq != 0): ## Deterministic 
            i_est_stn = -0.15*sig.square(2*np.pi*freq*t/1000, duty = 0.3*freq/1000) - 0.15*np.ones(len(t))
            ##i_est_stn = -9.0/2000.0*sig.square(2*np.pi*freq*t/1000, duty = 3.59*freq/1000) - 9.0/2000.0*np.ones(len(t))

    ## i_est_th is the an arbitrary current stimulation in the TH cells
    elif (stim_type == 'TH'):
        if (std != 0) and (freq == 0): ## Stochastic
            i_est_th = np.random.normal(0, std, len(t))
        elif (std == 0) and (freq != 0): ## Deterministic
            i_est_th = -0.15*sig.square(2*np.pi*freq*t/1000, duty = 0.3*freq/1000) - 0.15*np.ones(len(t))

    ## i_ctx is an pulsed stimulation current in the RS cells
    if (cortstim == 1):
        i_ctx = 175*sig.square(2*np.pi*t*cortstim/1000, duty = 0.3*cortstim/1000) + 175*np.ones(len(t))
        i_ctx[0:31] = np.zeros(31)
        i_ctx[-1] = 0
    #i_ctx = 3*np.ones(len(t))
    

    ## i_str is an sinusoidal stimulation current in the Str cells
    if str_stim == True:
        amp = (ac_str_stim*1e-3) # uA/cm2
        i_est_str = -amp*np.sin(2*5*np.pi*t/1000)

    
    ## Convert to HOC variables
    t_h = sim.h.Vector(len(t))
    i_est_th_h = sim.h.Vector(len(t))
    i_est_stn_h = sim.h.Vector(len(t))
    i_ctx_h = sim.h.Vector(len(t))
    i_est_str_h = sim.h.Vector(len(t))
    for i in range (0,len(t)):
        t_h.x[i] = t[i]
        i_est_th_h.x[i] = i_est_th[i]
        i_est_stn_h.x[i] = i_est_stn[i]
        i_ctx_h.x[i] = i_ctx[i]
        i_est_str_h.x[i] = i_est_str[i]
    for i in range(0,10): # inject these HOC variables to the variables defined in the MOD files
        #print(cells[i+70].secs.STN.hObj, '\n')
        i_est_stn_h.play(cells[i+70].secs.STN.hObj(0.5)._ref_i_est_STN, t_h, 1) # e.g. i_est_STN is a variable in STN.mod
        i_est_th_h.play(cells[i+20].secs.th.hObj(0.5)._ref_i_est_thalamus, t_h, 1)
        i_est_str_h.play(cells[i+10].secs.StrD2.hObj(0.5)._ref_i_est_Str, t_h, 1)
        i_est_str_h.play(cells[i].secs.StrD1.hObj(0.5)._ref_i_est_Str, t_h, 1)
        i_ctx_h.play(cells[i+50].secs.rs.pointps.Izhi.hObj._ref_i_est, t_h, 1)

    #for i, pop in enumerate(netParams.popParams.keys()):
    #    netParams.popParams[pop]['yRange'] = [i * 100, i * 100 + 10]


    ##### Simulate and analyze #####
    #sim.simulate()
    #sim.analyze()
    import pylab
    pylab.show()
    
    from scipy.signal import welch
    import matplotlib.pyplot as plt
    import numpy as np

    lfp = np.array(sim.allSimData['LFP'])  # shape: [num_eletrodos, num_tsteps]
    lfp = np.transpose(lfp, [1, 0])
    fs = 1.0 / simConfig.dt  # ex: 20 kHz → 0.05 ms -> fs = 20 kHz

    # Plot do sinal bruto
    plt.figure(figsize=(10, 3))
    plt.plot(lfp[0])
    plt.title("LFP - Eletrodo 0")
    plt.xlabel("Tempo (passos)")
    plt.ylabel("Potencial (uV)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('teste1.png')

    # Espectro de potência via Welch
    f, Pxx = welch(lfp[0], fs=fs, nperseg=1024)
    plt.figure(figsize=(8, 4))
    plt.semilogy(f, Pxx)
    plt.xlim([0, 100])
    plt.title('Espectro de Potência - LFP (Eletrodo 0)')
    plt.xlabel('Frequência (Hz)')
    plt.ylabel('Potência (uV^2/Hz)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('teste2.png')

    # Potência Beta (13–30 Hz)
    beta_power = np.trapz(Pxx[(f >= 13) & (f <= 30)], f[(f >= 13) & (f <= 30)])
    print(f"Potência beta: {beta_power:.4f} uV^2")

    print(sim.allSimData.keys())  # deve conter 'LFP'
    print(len(sim.allSimData['LFP']))
    


    #####################################################################################################
    ################################### 4 - Recording Data ##############################################
    #####################################################################################################


    #####################################################################################################
    ## 3.I - Extract Spikes

    class spikes:
        def __init__(self):
            self.times = []

    ## The vectors below carry the spike times of each nucleous. In each vector, there are 10 fields for
    ## the spike times of each cell of the corresponding nucleous.
    dStr_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()]
    iStr_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()]
    TH_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()]
    GPi_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()]
    GPe_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()]
    Cor_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), ##rs
               spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()] ##fs
    STN_APs = [spikes(), spikes(), spikes(),spikes(), spikes(), spikes(), spikes(),spikes(), spikes(), spikes()]

    ## Extract data from allSimData vector generated by the simulation and put the data in the vectors declared above
    for i in range(0,len(sim.allSimData.spkt)):
        if (sim.allSimData.spkid[i] >= 0 and sim.allSimData.spkid[i] <= 9):
            dStr_APs[int(sim.allSimData.spkid[i])].times = dStr_APs[int(sim.allSimData.spkid[i])].times+[sim.allSimData.spkt[i]]
        elif(sim.allSimData.spkid[i] >= 10 and sim.allSimData.spkid[i] <= 19):
            iStr_APs[int(sim.allSimData.spkid[i]-10)].times = iStr_APs[int(sim.allSimData.spkid[i]-10)].times+[sim.allSimData.spkt[i]]
        elif(sim.allSimData.spkid[i] >= 20 and sim.allSimData.spkid[i] <= 29):
            TH_APs[int(sim.allSimData.spkid[i]-20)].times = TH_APs[int(sim.allSimData.spkid[i]-20)].times+[sim.allSimData.spkt[i]]        
        elif(sim.allSimData.spkid[i] >= 30 and sim.allSimData.spkid[i] <= 39):
            GPi_APs[int(sim.allSimData.spkid[i]-30)].times = GPi_APs[int(sim.allSimData.spkid[i]-30)].times+[sim.allSimData.spkt[i]]
        elif(sim.allSimData.spkid[i] >= 40 and sim.allSimData.spkid[i] <= 49):
            GPe_APs[int(sim.allSimData.spkid[i]-40)].times = GPe_APs[int(sim.allSimData.spkid[i]-40)].times+[sim.allSimData.spkt[i]]
        elif(sim.allSimData.spkid[i] >= 50 and sim.allSimData.spkid[i] <= 69):
            Cor_APs[int(sim.allSimData.spkid[i]-50)].times = Cor_APs[int(sim.allSimData.spkid[i]-50)].times+[sim.allSimData.spkt[i]]
        elif(sim.allSimData.spkid[i] >= 70 and sim.allSimData.spkid[i] <= 79):
            STN_APs[int(sim.allSimData.spkid[i]-70)].times = STN_APs[int(sim.allSimData.spkid[i]-70)].times+[sim.allSimData.spkt[i]]            


    #####################################################################################################
    ## 3.II - Extract spike frequency of each nucleus
            
    freq_disp = [sim.allSimData.popRates['CTX_FSI'], sim.allSimData.popRates['CTX_RS'], sim.allSimData.popRates['GPe'],
                 sim.allSimData.popRates['GPi'], sim.allSimData.popRates['STN'], sim.allSimData.popRates['StrD1'],
                 sim.allSimData.popRates['StrD2'], sim.allSimData.popRates['TH']]
    for i in range(0,8):
        freq_disp[i] = round(freq_disp[i],2)

    # Após sim.runSim()

    # Plotar traces
    #fig = sim.analysis.plotTraces(include=[('STN', [0])])
#
    ## Salvar figura manualmente
    #if fig is not None:
    #    fig.savefig('./results/traces.png', dpi=300)

    # Plotar raster
    #fig_raster = sim.analysis.plotRaster(include=['CTX_RS', 'CTX_FSI', 'TH', 'GPi', 'GPe', 'STN', 'StrD2', 'StrD1'])
#
    #if fig_raster is not None:
    #    fig_raster.savefig('./results/raster.png', dpi=300)
#
    ## Plotar LFP manualmente
    #lfp_data = sim.allSimData['LFP']
#
    #plt.figure(figsize=(10, 4))
    #plt.plot(lfp_data[0])
    #plt.title('LFP signal at Electrode 0')
    #plt.xlabel('Time step')
    #plt.ylabel('LFP (uV)')
    #plt.grid(True)
    #plt.savefig('./results/LFP_signal.png', dpi=300)
    #plt.close()

    #####################################################################################################
    ## 3.III - Save the data in a .mat file
        
    if (pd == 0):
        name = 's' + str(it_num)
    elif(pd == 1 and stim_type == 'None'):
        name = 'pd' + str(it_num)
    elif(stim_type != 'None'):
        name = 'dbs' + str(it_num)
    if(freq != 0 and stim_type != 'None'):
        name = name + '_det'
    elif(freq == 0 and stim_type != 'None'):
        name = name + '_est_' + f_range
    if(str_stim == True):
        name = 'strstim_' + name
        
    
    sio.savemat(name, {
        'gsngi' : gsngi,
        'Striat_APs_dr': dStr_APs,
        'Striat_APs_indr': iStr_APs,
        'TH_APs': TH_APs,
        'GPi_APs': GPi_APs,
        'GPe_APs': GPe_APs,
        'Cor_APs': Cor_APs,
        'STN_APs': STN_APs,
        'freq_disp': freq_disp
    })
