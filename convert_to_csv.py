# Converts the HR5 catalog written in IDL format to csv
# Puts the files with _xxx.csv extension in the folder 
# outfolder 

from os import listdir
from os.path import isfile, join
from multiprocessing import Pool
from pathlib import Path

import tarfile

def untar_file(file):
        
    # open file
    file = tarfile.open(file)
    
    # extracting file
    file.extractall('/scratch/ankitsingh/Galaxy_catalogs/')
    
    file.close()

def convert_files(snapfile):

    if snapfile.endswith('.txt'):

        snapnum = snapfile.split('_')[3].split('.txt')[0]

        with open(snapfile, 'r') as fin:
            data = fin.read().splitlines(True)

        file = snapfile.split('/')[-1].replace('txt','csv')
        

        with open('/scratch/ankitsingh/Galaxy_catalogs/'+file, 'w') as fout:
            
            data = data[2:]

            for i in range(15,1,-1):
                
                data = [item.replace(i*' ',(i-1)*' ') for item in data]
            
            data = [item.replace(' ',',') for item in data]

            header = 'ID,HostHaloID,HostMtot(Msun),HostM200c(Msun),'\
                'x(cMpc),y(cMpc),z(cMpc),Vx(km/s),Vy(km/s),Vz(km/s),'\
                'Mstar(Msun),Mstar(r<5R1/2),Mgas(Msun),Mgas(r<5R1/2),'\
                'Mcold(T<10000K),Mcold(r<5R1/2),SFR_100Myr,SFR_100Myr(r<5R1/2),'\
                'SFR_20Myr,SFR_20Myr(r<5R1/2),R1/2(M*3D)(kpc),R1/2(M*2D)(kpc),'\
                'R1/2(M*2DFace-on),R1/2(r2D),R1/2(r2DFace-on),R50(Mgas3D)(kpc),'\
                'Z(star),Z(gas),Age_m(Gyr),u(Lsun),g(Lsun),r(Lsun),i(Lsun),z(Lsun),'\
                'Age_u(Gyr),Age_g(Gyr),Age_r(Gyr),Age_i(Gyr),Age_z(Gyr),Msmbh(Msun),'\
                'Mdm(Msun),Vrot(km/s),Vsig(km/s),Vrot(r<R1/2),Vsig(r<R1/2),Lx_star,'\
                'Ly_star,Lz_star,Lx_gas,Ly_gas,Lz_gas,Theta(Lstar-Lgas),'\
                'Compactness(u),Compactness(g),Compactness(r),Compactness(i),'\
                'Compactness(z),C(u),C(g),C(r),C(i),C(z),pure,ID_prog,flag_prog,'\
                'ID_prog_gen,step_prog_gen,f_inheritance,f_new,f_accrete,ID_descen,'\
                'flag_descen,ID_descen_gen,step_descen_gen,f_succession,host_flag\n'

            data.insert(0,header)

            print(f'writing file galaxy_catalog_{snapnum}.csv')
            fout.writelines(data)




try:

        hr5cat = '/scratch/ankitsingh/Galaxy_catalogs/'
        outfolder = hr5cat

        hr5files = sorted([join(hr5cat,f) for f in listdir(hr5cat) if (isfile(join(hr5cat,f))& f.endswith('.gz')) ])

        pool = Pool(4)
        
        pool.map(untar_file, hr5files)

        
            
        hr5files = sorted([join(hr5cat,f) for f in listdir(hr5cat) if (f.startswith('sn')) ])
        txtfile = sorted([join(f,i) for f in hr5files for i in listdir(f) ])

        
        pool.map(convert_files, txtfile)

        

finally: # To make sure processes are closed in the end, even if errors happen
        pool.close()
        pool.join()
        