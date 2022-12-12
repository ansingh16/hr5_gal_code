# Selecting galaxies in HR5

This repository contains the code necessary for the analysis of Intra Cluster Light (ICL) from the Horizon-Run 5 (HR5) simulation. The codes extract the clusters of interest from the galaxy catalogues of the HR5 depending on a critical mass set by the user. The code then uses the GALFIND output data to extract the galaxies and unbound particles within the clusters. The code contains following files in the order of execution:

- **params.ini** : parameter file containg all the important paths and parameters necessary for running the codes.
- **Galselect.py** : Code for selecting galaxies of interest. Currently , we select galaxies which have no neighbour within 200ckpc.
- **Get_gal_data** : Code for extracting the galaxy data. Currently, position, velocity, mass and temperature (in case of gas particles) are extracted for each galaxy. 
- **HR5_galaxy.ipynb** : Notebook containing example with working with HR5_mod module


## Parameters

- **[Paths]**
    Fofdir : Path to FOF directory of HR5 
    galaxyfile: Path to generate galaxies.csv
    galaxycats: Path to galaxy catalogs 
    outdir: Path to output directory

- **[Setting]**
    redshift: Redshift of the output study
    galmass: Critical mass above which the clusters are considered
    boxlen: boxlength of the simulation

The output file named `HR5_galaxies.hdf5` will be generated at the path `outdir` set in the params.ini file. The ouput file has the following structure:

```
HR5_galaxies
|
|___GalaxyID
    |
    |__ posstar: position of stellar particles in galaxies
    |   posdm: position of dark matter particles in galaxies
    |   posgas: position of gas particles in galaxies
    |   velstar: velocities of stellar particles in galaxies
    |   veldm: velocities of dark matter particles in galaxies
    |   velgas: velocities of gas particles in galaxies
    |   massstar: mass of stellar particles in galaxies
    |   massdm: mass of dark matter particles in galaxies
    |   massgas: mass of gas particles in galaxies
    |   tempg: temperature of gas particles in galaxies
    |
    |__ Halo_nsub: no. of subhalo
        Halo_nstar: no. of stars particles
        Halo_nsink: no. of sink particles
        Halo_ngas: no. of gas particles
        Halo_mtot: total mass
        Halo_mdm: total dark matter mass
        Halo_mgas: total gas mass
        Halo_msink: total sink mass
        Halo_mstar: total stellar mass
        Halo_pos: position of fof halo
        Halo_vel: velocity of fof halo
        nstar: no. of star particle in subhalo
        nsink: no. of sink particle in subhalo
        ngas: no. of gas particle in subhalo
        mtot: total mass of subhalo
        mdm: dark matter mass of of subhalo
        mgas: gas mass of of subhalo
        msink: sink mass of of subhalo
        mstar: stellar mass of of subhalo
        pos: of subhalo
        vel: of subhalo

```


