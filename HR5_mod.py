import numpy as np
import numpy as np
from struct import unpack
import pandas as pd 
import h5py
from scipy.io import FortranFile
import matplotlib.pyplot as plt 

import configparser

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('params.ini')

Fofd = parser.get('Paths','Fofdir')

red = float(parser.get('Setting','redshift'))


# Get snapshot corresponding to redshift
hr5outs = pd.read_csv('./Time_data.dat',dtype={'Snapshot':str,'Redshift':float,\
                    'LBT':float,'dx':float})
snapno = hr5outs.loc[hr5outs['Redshift']==red,'Snapshot'].values[0]




class Galaxy:
    """Class containg galaxy parameters"""
    def __init__(self):
        self.galtsub = {}
        self.galthalo = {}
        self.galgas = {}
        self.galsink = {}
        self.galdm = {}
        self.galstar = {}




class ICL:
    """Class containg galaxy parameters"""
    def __init__(self):
        self.icldict = {}
        self.iclgas = {}
        self.iclsink = {}
        self.icldm = {}
        self.iclstar = {}


# Function to get galaxy data
def Get_galaxy(ID,HostHaloID,snapno):

    Mygal=Galaxy()

    kkk=0
    # Open Galaxy find data files
    with open(f'{Fofd}/FoF.{snapno:0>5}/GALFIND.DATA.{snapno:0>5}', mode='rb') as file: # b -> binary

        hline=0
        sline=0
        
        found=False
        while True:
            
            fof = file.read(112)
            if (not fof):
                print("Done!!")
                break
            else:
                    # Read fof halo
                    data1 = unpack('@6i11d',fof)
                    
                    thalo={'nsub':data1[0],'ndm': data1[1],'nstar': data1[2],'nsink':data1[3],'ngas':data1[4],'npall':data1[5],'mtot':data1[6],\
                    'mdm':data1[7],'mgas':data1[8],'msink':data1[9],'mstar':data1[10],'pos':data1[11:14],'vel':data1[14:17]}
                    
                    # Read subhaloes
                    for i in range(thalo['nsub']):
                            
                            data2 = unpack('@6i11d',file.read(112))
                            tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                    'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}

                            if  (hline == HostHaloID)&(sline==ID):

                                    print(hline,sline)      

                                    
                                    dm=[]
                                    # Read dm particles
                                    for j in range(tsub['ndm']):

                                        data3 = unpack('@13d1q1d1i1f',file.read(128))
                                        tdm={'pos':data3[0:3],'vel':data3[3:6],'mass':data3[6],'dum0':data3[7],'tp':data3[8],'zp':data3[9],'mass0':data3[10],\
                                        'tpp':data3[11],'indtab':data3[12],'id':data3[13],'potential':data3[14],'level':data3[15],'dum1':data3[16]}      

                                        dm.append(tdm)

                                    
                                    if tsub['ngas']>0:
                                        gas=[]   
                                        # Read gas particles
                                        for j in range(tsub['ngas']):
                                                #print('gas',j,tsub)
                                                

                                            data4 = unpack('@4d4f1d5f1i2f1q4d',file.read(128))
                                            tgas={'pos':data4[0:3],'dx':data4[3],'vel':data4[4:7],'dum0':data4[7],'density':data4[8],'temp':data4[9],'metal':data4[10],'fe':data4[11],\
                                                    'h':data4[12],'o':data4[13],'level':data4[14],'mass':data4[15],'dum1':data4[16],'id':data4[17],'potential':data4[18],'f':data4[19:22]}
                                            
                                            gas.append(tgas)


                                    if tsub['nsink']>0:

                                        sink=[]
                                        # Read sink particles
                                        for j in range(tsub['nsink']):
                                            data5 = unpack('@20d2i',file.read(168))
                                            tsink={'pos':data5[0:3],'vel':data5[3:6],'mass':data5[6],'tbirth':data5[7],'angm':data5[8:11],'ang':data5[11:14],'dmsmbh':data5[14:17],\
                                                    'esave':data5[17],'smag':data5[18],'eps':data5[19],'id':data5[20],'dum0':data5[21]}

                                            sink.append(tsink)

                                    
                                    if tsub['nstar']>0:

                                        star=[]   
                                        # Read star particles
                                        for j in range(tsub['nstar']):
                                            data6 = unpack('@13d1q1d1i1f',file.read(128))
                                            tstar={'pos':data6[0:3],'vel':data6[3:6],'mass':data6[6],'dum0':data6[7],'tp':data6[8],'zp':data6[9],'mass0':data6[10],\
                                            'tpp':data6[11],'indtab':data6[12],'id':data6[13],'potential':data6[14],'level':data6[15],'dum1':data6[16]}

                                            star.append(tstar)



                                    if tsub['nsink']>0:
                                        Mygal.galsink = pd.DataFrame(sink)
                                    if tsub['ngas']>0:
                                        Mygal.galgas = pd.DataFrame(gas)
                                    if tsub['nstar']>0:
                                        Mygal.galstar = pd.DataFrame(star)
                                    
                                    Mygal.galdm = pd.DataFrame(dm)
                                    
                                    Mygal.galthalo = thalo
                                    Mygal.galtsub = tsub 

                                    print(f"Galaxy Found!!!{sline}")
                                    return Mygal
                            
                            else:
                                    # Read dm particles
                                    for j in range(tsub['ndm']):               
                                        file.read(128)
                                        
                                    # Read gas particles
                                    for j in range(tsub['ngas']):
                                        file.read(128)

                                    if tsub['nsink']>0:
                                        # Read sink particles
                                        for j in range(tsub['nsink']):
                                            file.read(168)
                                            
                                    if tsub['nstar']>0:
                                        for j in range(tsub['nstar']):
                                            file.read(128)
                                
                                    sline = sline + 1

                    hline=hline+1

    


def Get_Unbound(HostHaloID):

    Myunbound = ICL()

    # Open Galaxy find data files
    with open(f'{Fofd}/FoF.{snapno:0>5}/background_ptl.{snapno:0>5}', mode='rb') as file: # b -> binary

        hline=0
        
        while True:
            
            fof = file.read(112)
            if not fof :
                print("Done!! Smell the success.. :)")
                break
            else:
                
                # Check if the cluster is of interest
                if (hline == HostHaloID):
                                
                    # Read fof halo
                    data2 = unpack('@6i11d',fof)
                                    
                    tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                    'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}


                    dm=[]

                    for j in range(tsub['ndm']):

                                                
                        data3 = unpack('@13d1q1d1i1f',file.read(128))
                        tdm={'pos':data3[0:3],'vel':data3[3:6],'mass':data3[6],'dum0':data3[7],'tp':data3[8],'zp':data3[9],'mass0':data3[10],\
                        'tpp':data3[11],'indtab':data3[12],'id':data3[13],'potential':data3[14],'level':data3[15],'dum1':data3[16]} 
                        
                        dm.append(tdm)
                                            


                    gas=[]
                    if tsub['ngas']>0:
                        # Read gas particles
                        for j in range(tsub['ngas']):
                            #print('gas',j,tsub)
                            
                            data4 = unpack('@4d4f1d5f1i2f1q4d',file.read(128))
                            tgas={'pos':data4[0:3],'dx':data4[3],'vel':data4[4:7],'dum0':data4[7],'density':data4[8],'temp':data4[9],'metal':data4[10],'fe':data4[11],\
                                'h':data4[12],'o':data4[13],'level':data4[14],'mass':data4[15],'dum1':data4[16],'id':data4[17],'potential':data4[18],'f':data4[19:22]}
                            
                            gas.append(tgas)
                                                
                    sink=[]
                    if tsub['nsink']>0:
                            # Read sink particles
                            for j in range(tsub['nsink']):
                                data5 = unpack('@20d2i',file.read(168))
                                tsink={'pos':data5[0:3],'vel':data5[3:6],'mass':data5[6],'tbirth':data5[7],'angm':data5[8:11],'ang':data5[11:14],'dmsmbh':data5[14:17],\
                                        'esave':data5[17],'smag':data5[18],'eps':data5[19],'id':data5[20],'dum0':data5[21]}

                                sink.append(tsink)
                    
                    star=[]
                    if tsub['nstar']>0:
                                                
                                                
                            # Read star particles
                            for j in range(tsub['nstar']):

                                data6 = unpack('@13d1q1d1i1f',file.read(128))
                                tstar={'pos':data6[0:3],'vel':data6[3:6],'mass':data6[6],'dum0':data6[7],'tp':data6[8],'zp':data6[9],'mass0':data6[10],\
                                'tpp':data6[11],'indtab':data6[12],'id':data6[13],'potential':data6[14],'level':data6[15],'dum1':data6[16]}
                                
                                star.append(tstar)

                    

                    if tsub['nsink']>0:
                        Myunbound.iclsink = pd.DataFrame(sink)
                    if tsub['ngas']>0:
                        Myunbound.iclgas = pd.DataFrame(gas)
                    if tsub['nstar']>0:
                        Myunbound.iclstar = pd.DataFrame(star)
                                    
                    Myunbound.icldm = pd.DataFrame(dm)
                                    
                    Myunbound.icldict = tsub 

                    print(f"FOF found : {hline}")
                    return Myunbound



                else:
                    # Read fof halo
                    data2 = unpack('@6i11d',fof)
                                    
                    tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                    'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}
                
                    #print(tsub)


                    # Read dm particles
                    for j in range(tsub['ndm']):               
                        file.read(128)
                                
                    # Read gas particles
                    for j in range(tsub['ngas']):
                        file.read(128)

                    if tsub['nsink']>0:
                        # Read sink particles
                        for j in range(tsub['nsink']):
                            file.read(168)
                                    
                    if tsub['nstar']>0:
                        for j in range(tsub['nstar']):
                            file.read(128)
                        

                    hline=hline+1
                    
                            

