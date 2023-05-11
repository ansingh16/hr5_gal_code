from scipy.io import FortranFile
import pandas as pd 



import configparser
config = configparser.ConfigParser()


config.read('start_setup.ini')

simulation_dir = config['sim_dir']['simulation_dir']
snapnum = config['sim_dir']['snapnum']


def get_dt():
    info={}
    with open(f'{simulation_dir}output_00{snapnum}/info_00{snapnum}.txt') as f:
        for l,txt in enumerate(f):
            if l<=5:
                lab = txt.split('=')[0].rstrip()
                #print(int(txt.split('=')[1].split('\n')[0]))
                info[lab] = int(txt.split('=')[1].split('\n')[0])
                
            if (l>6)&(l<18):
                lab = txt.split('=')[0].rstrip()
                #print(float(txt.split('=')[1].split('\n')[0]))
                info[lab] = float(txt.split('=')[1].split('\n')[0])

    info['scale_m']    = info['unit_d']*(info['unit_l'])**3


    f = FortranFile(f'{simulation_dir}output_00{snapnum}/amr_00{snapnum}.out00001', 'r' )

    ncpu = f.read_ints( dtype='i4')
    ndim = f.read_ints( dtype='i4')
    nx,ny,nz = f.read_ints( dtype='i4')
    nlevelmax = f.read_ints( dtype='i4')
    ngridmax = f.read_ints( dtype='i4')
    nboundary = f.read_ints( dtype='i4')
    ngrid_current = f.read_ints( dtype='i4')
    boxlen = f.read_reals( dtype='float')

    noutput,iout,ifout = f.read_ints( dtype='i4')
    tout = f.read_reals( dtype='float')
    aout = f.read_reals( dtype='float')
    t = f.read_reals( dtype='float')
    dtold = f.read_reals( dtype='float')
    dtnew = f.read_reals( dtype='float')
    nstep,nstep_coarse = f.read_ints( dtype='i4')
    print(nstep,nstep_coarse)

    f.close()
    
    return info,dtnew 






def Read_sink(file):
        f = FortranFile( file, 'r' )

        nsink = f.read_ints( dtype='i4')
        nindsink = f.read_ints( dtype='i4') 
        idsink = f.read_ints( dtype='i4')
        msink = f.read_reals( dtype='float')
        xsink = f.read_reals( dtype='float')
        ysink = f.read_reals( dtype='float')
        zsink = f.read_reals( dtype='float')
        vsinkx = f.read_reals( dtype='float')
        vsinky = f.read_reals( dtype='float')
        vsinkz = f.read_reals( dtype='float')
        tsink = f.read_reals( dtype='float')
        dMsmbh = f.read_reals( dtype='float')
        dMBondi = f.read_reals( dtype='float')
        dMEdd = f.read_reals( dtype='float')
        esave = f.read_reals( dtype='float')
        Lx_gas = f.read_reals( dtype='float')
        Ly_gas = f.read_reals( dtype='float')
        Lz_gas = f.read_reals( dtype='float')
        Lx_bh = f.read_reals( dtype='float')
        Ly_bh = f.read_reals( dtype='float')
        Lz_bh = f.read_reals( dtype='float')
        Lbh = f.read_reals( dtype='float')
        Eff_bh = f.read_reals( dtype='float')

        d = {'xsink':xsink,'ysink':ysink,'zsink':zsink,\
            'idsink':idsink,'dMsmbh':dMsmbh,'dMBondi':dMBondi,\
                'dMEdd':dMEdd,'Eff_bh':Eff_bh,'msink':msink}

        Sinkout = pd.DataFrame(d)
        Sinkout.astype({'idsink': 'int32'})

        return Sinkout
