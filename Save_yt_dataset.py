from telnetlib import GA
import numpy as np
import numpy as np
import pandas as pd 
import h5py
import yt 
from yt.units import Msun, Mpc

import matplotlib.pyplot as plt 
plt.style.use('/home/ankitsingh/hr5_agn/paper_style.mplstyle')


import configparser

# load up the parameter file
parser = configparser.ConfigParser()
parser.read('params.ini')


outdir = parser.get('Paths','outdir')
plotdir = parser.get('Paths','projplot')
galaxyfile = parser.get('Paths','galaxyfile')


def galaxy_projection(HostID,galID):
    
    # Function for making projection plots of the clusters
    # It reads the data from the hdf5 files and saves the plots
    # in the plotdir path defined in params.ini

    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'Y (cMpc)')
    ax.set_ylabel(r'Z (cMpc)')
    
    # Get the position of unbound stars
    pos = fout[f'{HostID}/{galID}/posstar']

    ax.scatter(pos[:,1],pos[:,2],s=0.2,color='r',alpha=0.4)

    # stars in galaxies are to be plotted in blue
    pos =  fout[f'{HostID}/{gal}/posstar']
    ax.scatter(pos[:,1],pos[:,2],s=0.2,color='b',alpha=0.1)
    
    #fig.savefig(f'{plotdir}Proj{clusID}.png')




def analysis_yt(HostID):

    # This function saves yt dataset from the hdf5 data.
    # It uses the yt for loading a generic particle dataset of stars, gas and dm


    for gal in fout[f'/{HostID}/'].keys():

        galatt = dict(fout[f'/{HostID}/{gal}/'].attrs.items())

        # Stars
        poss =  fout[f'/{HostID}/{gal}/posstar']
        mstar =  fout[f'/{HostID}/{gal}/massstar']

        # DM
        posdm =  fout[f'/{HostID}/{gal}/posdm']
        mdm =  fout[f'/{HostID}/{gal}/massdm']


        # Gas
        posg =  fout[f'/{HostID}/{gal}/posgas']
        mgas =  fout[f'/{HostID}/{gal}/massgas']
        


        data_all = {
                        ("star","particle_position_x"): poss[:,0]- galatt['pos'][0],
                        ("star","particle_position_y"): poss[:,1]- galatt['pos'][1],
                        ("star","particle_position_z"): poss[:,2]- galatt['pos'][2],
                        ("star","particle_mass"): mstar,
                        ("dm","particle_position_x"): posdm[:,0]- galatt['pos'][0],
                        ("dm","particle_position_y"): posdm[:,1]- galatt['pos'][1],
                        ("dm","particle_position_z"): posdm[:,2]- galatt['pos'][2],
                        ("dm","particle_mass"): mdm,
                        ("gas","particle_position_x"): posg[:,0]- galatt['pos'][0],
                        ("gas","particle_position_y"): posg[:,1]- galatt['pos'][1],
                        ("gas","particle_position_z"): posg[:,2]- galatt['pos'][2],
                        ("gas","particle_mass"): mgas
                    }
        
        # get width that can encompass all particles
        width=0.2
        result = None
        while result is None:
            try:
                # connect
                bbox = np.array([[-width,width], [-width, width], [-width, width]])

                ds_all = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
                result = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
                
            except:
                width=width+0.1
                pass
        
        ds = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
        
        print(HostID,gal)
        
        ######### do something with this data#########


fout = h5py.File(f"{outdir}HR5_galaxies.hdf5", "r")

galaxies = pd.read_csv(galaxyfile,usecols=['ID','HostHaloID'])


# YT dataset analysis
for i in range(galaxies.shape[0]):    
    analysis_yt(galaxies['HostHaloID'].iloc[i])
    
for i in range(galaxies.shape[0]):    
    galaxy_projection(galaxies['HostHaloID'].iloc[i],galaxies['ID'].iloc[i])



