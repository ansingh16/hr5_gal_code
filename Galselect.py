import numpy as np
from struct import unpack
import pandas as pd 
from scipy.io import FortranFile
import configparser
from periodic_kdtree import PeriodicCKDTree


# load up the parameter file
parser = configparser.ConfigParser()
parser.read('params.ini')

# Read the parameters
galdir = parser.get('Paths','galaxycats')
galaxyfile = parser.get('Paths','galaxyfile')

red = float(parser.get('Setting','redshift'))
critmass = float(parser.get('Setting','galmass'))
boxlen = float(parser.get('Setting','boxlen'))

bounds = np.array([boxlen,-1,-1])

# Get snapshot corresponding to redshift
hr5outs = pd.read_csv('./Time_data.dat',dtype={'Snapshot':str,'Redshift':float,\
                     'LBT':float,'dx':float})
snapno = hr5outs.loc[hr5outs['Redshift']==red,'Snapshot'].values[0]

# Read the galaxy catalogue
file = f"{galdir}galaxy_catalogue_{snapno}.csv"

allgal = pd.read_csv(file,usecols=['Mstar(Msun)','ID','HostHaloID','pure','HostMtot(Msun)','x(cMpc)','y(cMpc)','z(cMpc)'])



data = allgal[['x(cMpc)','y(cMpc)','z(cMpc)']].to_numpy()
zoomgal = allgal.loc[allgal['pure']==1]

galpos = zoomgal[['x(cMpc)','y(cMpc)','z(cMpc)']].to_numpy()

T = PeriodicCKDTree(bounds,data)
d, i = T.query(galpos, k=2)
zoomgal['NND'] = d[:,1]


# Get the galaxies of interest
# Only isolated galaxies
isolated_Galaxies = zoomgal.loc[zoomgal['NND']>0.2] # no galaxies within 200kpc

print(isolated_Galaxies.shape)
galaxies = isolated_Galaxies.loc[isolated_Galaxies['Mstar(Msun)']>critmass,['ID','HostHaloID','HostMtot(Msun)','Mstar(Msun)']]

# # Move the ClusterIDs to the file
print(f"saving to : {galaxyfile}")
galaxies.to_csv(galaxyfile,index=False)

