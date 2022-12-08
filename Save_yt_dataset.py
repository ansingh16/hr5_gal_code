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
galaxyfilefile = parser.get('Paths','galaxyfile')


def galaxy_projection(galID):
    
    # Function for making projection plots of the clusters
    # It reads the data from the hdf5 files and saves the plots
    # in the plotdir path defined in params.ini

    fig,ax = plt.subplots(1,1)
    ax.set_xlabel(r'Y (cMpc)')
    ax.set_ylabel(r'Z (cMpc)')
    
    # Get the position of unbound stars
    pos = fout[f'/{galID}/posstar']

    ax.scatter(pos[:,1],pos[:,2],s=0.2,color='r',alpha=0.4)

    # stars in galaxies are to be plotted in blue
    pos =  f[f'/{clusID}/{gal}/posstar']
    ax.scatter(pos[:,1],pos[:,2],s=0.2,color='b',alpha=0.1)
    
    fig.savefig(f'{plotdir}Proj{clusID}.png')




def analysis_yt(galID):

    # This function saves yt dataset from the hdf5 data.
    # It uses the yt for loading a generic particle dataset of stars, gas and dm


    galatt = dict(f[f'/{galID}/'].attrs.items())

    # Initialising global arrays for position and mass
    # of DM and stars
    
    # Stars
    poss =  f[f'/{galID}/posstar']
    mstar =  f[f'/{galID}/massstar']
    posstarx = np.array([])
    posstary = np.array([])
    posstarz = np.array([])
    
    # DM
    posdm =  f[f'/{galID}/posdm']
    mdm =  f[f'/{galID}/massdm']
    posdmy = np.array([])
    posdmz = np.array([])
    posdmx = np.array([])
    

    # gas
    posg =  f[f'/{galID}/posgas']
    mgas =  f[f'/{galID}/massgas']
    posgasy = np.array([])
    posgasz = np.array([])
    posgasx = np.array([])
    


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
    res=1024
    width=0.5
    result = None
    while result is None:
        try:
            # connect
            bbox = np.array([[-width,width], [-width, width], [-width, width]])

            ds_all = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
            result = yt.ParticleProjectionPlot(ds_all,'x',("star","particle_mass"))
            
        except:
            width=width+0.5
            pass
    
    ds = yt.load_particles(data_all, length_unit='Mpc', mass_unit='Msun', bbox=bbox)
    
    ######### do something with this data#########


fout = h5py.File(f"{outdir}HR5_galaxies.hdf5", "r")

galaxies = pd.read_csv(galaxyfile)


for i in range(galaxies.shape[0]):    
    galaxy_projection(galaxies.iloc[i].values[0])
    analysis_yt(galaxies.iloc[i].values[0])
    

