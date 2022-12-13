from telnetlib import GA
import numpy as np
from scipy import interpolate,linalg
import numpy as np
from struct import unpack
import pandas as pd 
import h5py
from scipy.io import FortranFile
import matplotlib.pyplot as plt 
# plt.style.use('/home/ankitsingh/HR5_AGN/paper_style.mplstyle')

import configparser

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('params.ini')

Fofd = parser.get('Paths','Fofdir')

galaxyfile = parser.get('Paths','galaxyfile')
output = parser.get('Paths','outdir')


red = float(parser.get('Setting','redshift'))


# Get snapshot corresponding to redshift
hr5outs = pd.read_csv('./Time_data.dat',dtype={'Snapshot':str,'Redshift':float,\
                    'LBT':float,'dx':float})
snapno = hr5outs.loc[hr5outs['Redshift']==red,'Snapshot'].values[0]



galaxies = pd.read_csv(galaxyfile,usecols=['ID'])


with h5py.File(f"{output}HR5_galaxies.hdf5", "w") as fout:

    # Open Galaxy find data files
    with open(f'{Fofd}/FoF.{snapno:0>5}/GALFIND.DATA.{snapno:0>5}', mode='rb') as file: # b -> binary

        hline=0
        sline=0
        kkk=0

        while True:
            
            fof = file.read(112)
            if not fof or (kkk>galaxies.shape[0]-1):
                print("Done!!",kkk,sline)
                break
            else:
                    print(f'HaloID={hline},SubHaloID={sline},Galaxy={kkk}')
                
                                
                    # Read fof halo
                    data1 = unpack('@6i11d',fof)
                    
                    thalo={'nsub':data1[0],'ndm': data1[1],'nstar': data1[2],'nsink':data1[3],'ngas':data1[4],'npall':data1[5],'mtot':data1[6],\
                    'mdm':data1[7],'mgas':data1[8],'msink':data1[9],'mstar':data1[10],'pos':data1[11:14],'vel':data1[14:17]}
                    
                    
                    #print(thalo)

                    # Check if the cluster is of interest
                    if (kkk<=galaxies.shape[0]-1) & (sline == galaxies['ID'].iloc[kkk]):
                                            
                            # Read subhaloes
                            for i in range(thalo['nsub']):
                                    data2 = unpack('@6i11d',file.read(112))
                                    tsub={'ndm':data2[0],'ngas':data2[1],'nsink': data2[2],'nstar':data2[3],'npall': data2[4],'dum':data2[5],\
                                    'mtot':data2[6],'mdm':data2[7],'mgas':data2[8],'msink':data2[9],'mstar':data2[10],'pos':data2[11:14],'vel':data2[14:17]}

                                    #print(hline,thalo['pos'],tsub['pos'])
                                    
                                    posdm = np.empty((tsub['ndm'],3))
                                    massdm = np.empty(tsub['ndm'])
                                    veldm = np.empty((tsub['ndm'],3))
                                    # Read dm particles
                                    for j in range(tsub['ndm']):

                                        data3 = unpack('@13d1q1d1i1f',file.read(128))
                                        tdm={'pos':data3[0:3],'vel':data3[3:6],'mass':data3[6],'dum0':data3[7],'tp':data3[8],'zp':data3[9],'mass0':data3[10],\
                                        'tpp':data3[11],'indtab':data3[12],'id':data3[13],'potential':data3[14],'level':data3[15],'dum1':data3[16]}      

                                        posdm[j,:] = np.array(tdm['pos']) 
                                        veldm[j,:] = np.array(tdm['vel']) 
                                        massdm[j] = tdm['mass']
                                        

                                    # Data you are intrested in
                                    posg = np.empty((tsub['ngas'],3))
                                    massg = np.empty(tsub['ngas'])
                                    velg = np.empty((tsub['ngas'],3))
                                    tempg = np.empty(tsub['ngas'])
                    
                                        
                                    # Read gas particles
                                    for j in range(tsub['ngas']):
                                            #print('gas',j,tsub)
                                            

                                        data4 = unpack('@4d4f1d5f1i2f1q4d',file.read(128))
                                        tgas={'pos':data4[0:3],'dx':data4[3],'vel':data4[4:7],'dum0':data4[7],'density':data4[8],'temp':data4[9],'metal':data4[10],'fe':data4[11],\
                                                'h':data4[12],'o':data4[13],'level':data4[14],'mass':data4[15],'dum1':data4[16],'id':data4[17],'potential':data4[18],'f':data4[19:22]}

                                        
                                        # Fill the empty arrays
                                        posg[j,:] = np.array(tgas['pos']) 
                                        velg[j,:] = np.array(tgas['vel']) 
                                        massg[j] = tgas['mass']
                                        tempg[j] = tgas['temp']/tgas['density']

                                        

                                    if tsub['nsink']>0:
                                        # Read sink particles
                                        for j in range(tsub['nsink']):
                                            data5 = unpack('@20d2i',file.read(168))
                                            tsink={'pos':data5[0:3],'vel':data5[3:6],'mass':data5[6],'tbirth':data5[7],'angm':data5[8:11],'ang':data5[11:14],'dmsmbh':data5[14:17],\
                                                    'esave':data5[17],'smag':data5[18],'eps':data5[19],'id':data5[20],'dum0':data5[21]}


                                    posstar = np.empty((tsub['nstar'],3))
                                    massstar = np.empty(tsub['nstar'])
                                    velstar = np.empty((tsub['nstar'],3))

                                    if tsub['nstar']>0:
                                            
                                        # Read star particles
                                        for j in range(tsub['nstar']):
                                            data6 = unpack('@13d1q1d1i1f',file.read(128))
                                            tstar={'pos':data6[0:3],'vel':data6[3:6],'mass':data6[6],'dum0':data6[7],'tp':data6[8],'zp':data6[9],'mass0':data6[10],\
                                            'tpp':data6[11],'indtab':data6[12],'id':data6[13],'potential':data6[14],'level':data6[15],'dum1':data6[16]}

                                                
                                            posstar[j,:] = np.array(tstar['pos']) 
                                            velstar[j,:] = np.array(tstar['vel']) 
                                            massstar[j] = tstar['mass']
                            
                                        
                                    # Check if the dataset exists, delete it
                                    for dat in ['posstar','posdm','posgas','velstar',\
                                            'veldm','velgas','massstar','massdm','massgas','tgas']:

                                        if f'{hline}/{sline}/{dat}' in fout:
                                                del fout[f'{hline}/{sline}/{dat}']
                                    
                                    # Write the data from arrays to hdf5
                                    # write positions of stars in subhalo
                                    fout.create_dataset(f'{hline}/{sline}/posstar',data=posstar)
                                    fout.create_dataset(f'{hline}/{sline}/posdm',data=posdm)
                                    fout.create_dataset(f'{hline}/{sline}/posgas',data=posg)

                                    # write positions of stars in subhalo
                                    fout.create_dataset(f'{hline}/{sline}/velstar',data=velstar)
                                    fout.create_dataset(f'{hline}/{sline}/veldm',data=veldm)
                                    fout.create_dataset(f'{hline}/{sline}/velgas',data=velg)
                                        
                                    # Write mass
                                    fout.create_dataset(f'{hline}/{sline}/massstar',data=massstar)
                                    fout.create_dataset(f'{hline}/{sline}/massdm',data=massdm)
                                    fout.create_dataset(f'{hline}/{sline}/massgas',data=massg)


                                    # write temperature
                                    fout.create_dataset(f'{hline}/{sline}/tgas',data=tempg)
                                    
                                    clus = fout[f'/{hline}/']
                                    for par in ['nsub','nstar','nsink','ngas','mtot','mdm',\
                                            'mgas','msink','mstar','pos','vel']:

                                        clus.attrs[par] = thalo[par]

                                    gal = fout[f'/{hline}/{sline}']
                                    for par in ['nstar','nsink','ngas','mtot','mdm',\
                                            'mgas','msink','mstar','pos','vel']:
                                            gal.attrs[par] = tsub[par]
                                    
                                    
                                
                                    sline = sline + 1
                                
                            print(f"galaxy found!! {kkk}")

                                        
                            kkk=kkk+1


                                        
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

                
                    hline = hline + 1
