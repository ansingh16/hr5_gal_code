
# -*- coding: utf-8 -*-
import numpy as np

class domain(object):
    
    def __init__(self, shape):
        self.shape = shape
        
    @classmethod
    def read(cls, filename):
        
        f = open(filename,'r')
        line   = f.readline()
        shape  = line.strip()
        if shape=='sphere':
            line   = f.readline()
            center = np.array(line.split(),dtype=float)
            line   = f.readline()
            radius = np.float(line)
            f.close()
            return sphere(center,radius)

        elif shape=='shell':
            line   = f.readline()
            center = np.array(line.split(),dtype=float)
            line   = f.readline()
            rin    = np.float(line)
            line   = f.readline()
            rout   = np.float(line)
            f.close()
            return shell(center,rin,rout)

        elif shape=='cube':
            line   = f.readline()
            center = np.array(line.split(),dtype=float)
            line   = f.readline()
            size   = np.float(line)
            f.close()
            return cube(center,size)

        elif shape=='slab':
            line      = f.readline()
            zc        = np.float(line)
            line      = f.readline()
            thickness = np.float(line)
            f.close()
            return slab(zc,thickness)
        
        else:
            return IOError('bouh')
        

class sphere(domain):
    
    def __init__(self, center, radius):
        domain.__init__(self, 'sphere')
        self.radius = radius
        self.center = center
        
    def area(self):
        return 4*np.pi*self.radius*self.radius

    def volume(self):
        return 4./3.*np.pi*(self.radius)**3.

    # def info(self):
    #     ##print("domain INFO:")
    #     ##print("|_shape  =",self.shape)
    #     print("|_radius =",self.radius)
    #     print("|_center =",self.center)


class cube(domain):
    
    def __init__(self, center, size):
        domain.__init__(self, 'cube')
        self.center = center
        self.size   = size
        
    def area(self):
        return self.size*self.size

    def volume(self):
        return (self.size)**3.

    # def info(self):
    #     #print("domain INFO:")
    #     #print("|_shape  =",self.shape)
    #     print("|_size   =",self.size)
    #     print("|_center =",self.center)

class shell(domain):

    def __init__(self, center, rin, rout):
        domain.__init__(self, 'shell')
        self.center = center
        self.rin    = rin
        self.rout   = rout

    def area(self):
        return 4*np.pi*(self.rout**2-self.rin**2)

    def volume(self):
        return 4./3.*np.pi*(self.rout**3.-self.rin**3.)

    # def info(self):
    #     #print("domain INFO:")
    #     #print("|_shape  =",self.shape)
    #     #print("|_center =",self.center)
    #     print("|_rin    =",self.rin)
    #     print("|_rout   =",self.rout)

class slab(domain):

    def __init__(self, zc, thickness):
        domain.__init__(self, 'slab')
        self.zc = zc
        self.thickness = thickness

    # def info(self):
    #     #print("domain INFO:")
    #     #print("|_shape     =",self.shape)
    #     print("|_zc        =",self.zc)
    #     print("|_thickness =",self.thickness)


def overplot_limits(domain, color=None, linestyle=None, linewidth=None, projection=None):
    # OVERPLOT DOMAIN LIMITS
    from matplotlib.patches import Circle, Wedge, Polygon,Rectangle
    from matplotlib.collections import PatchCollection
    patches=[]
    if projection is None:
        projection = 'yx'
    
    if (domain.shape == 'sphere'):
        if (projection == 'yx'):
            xc = domain.center[0]
            yc = domain.center[1]
        elif (projection == 'xy'):
            xc = domain.center[1]
            yc = domain.center[0]
        elif (projection == 'zx'):
            xc = domain.center[0]
            yc = domain.center[2]
        elif (projection == 'xz'):
            xc = domain.center[2]
            yc = domain.center[0]
        elif (projection == 'yz'):
            xc = domain.center[2]
            yc = domain.center[1]
        elif (projection == 'zy'):
            xc = domain.center[1]
            yc = domain.center[2]
        wedge = Wedge((xc, yc), domain.radius, 0., 360.)

    elif (domain.shape == 'shell'):
        if (projection == 'yx'):
            xc = domain.center[0]
            yc = domain.center[1]
        elif (projection == 'xy'):
            xc = domain.center[1]
            yc = domain.center[0]
        elif (projection == 'zx'):
            xc = domain.center[0]
            yc = domain.center[2]
        elif (projection == 'xz'):
            xc = domain.center[2]
            yc = domain.center[0]
        elif (projection == 'yz'):
            xc = domain.center[2]
            yc = domain.center[1]
        elif (projection == 'zy'):
            xc = domain.center[1]
            yc = domain.center[2]
        wedge = Wedge((xc,yc), domain.rout, 0., 360., width=domain.rout-domain.rin)

    elif (domain.shape == 'cube'):
        if (projection == 'yx'):
            xc = domain.center[0]-domain.size/2.
            yc = domain.center[1]-domain.size/2.
        elif (projection == 'xy'):
            xc = domain.center[1]-domain.size/2.
            yc = domain.center[0]-domain.size/2.
        elif (projection == 'zx'):
            xc = domain.center[0]-domain.size/2.
            yc = domain.center[2]-domain.size/2.
        elif (projection == 'xz'):
            xc = domain.center[2]-domain.size/2.
            yc = domain.center[0]-domain.size/2.
        elif (projection == 'zy'):
            xc = domain.center[1]-domain.size/2.
            yc = domain.center[2]-domain.size/2.
        elif (projection == 'yz'):
            xc = domain.center[2]-domain.size/2.
            yc = domain.center[1]-domain.size/2.
        
        wedge = Rectangle((xc,yc),domain.size,domain.size)

    elif (domain.shape == 'slab'):
        if (projection == 'yx'):
                wedge = Rectangle((0.,0.), 1.,1.)
        elif (projection == 'xy'):
                wedge = Rectangle((0.,0.), 1.,1.)
        elif (projection == 'zx'):
            wedge = Rectangle((0., domain.zc-domain.thickness/2.),1.,domain.thickness)
        elif (projection == 'xz'):
            wedge = Rectangle((domain.zc-domain.thickness/2., 0.),domain.thickness, 1.)
        elif (projection == 'zy'):
            wedge = Rectangle((0., domain.zc-domain.thickness/2.),1.,domain.thickness)
        elif (projection == 'yz'):
            wedge = Rectangle((domain.zc-domain.thickness/2., 0.),domain.thickness, 1.)

    patches.append(wedge)
    if color is None:
        color='k'
    if linestyle is None:
        linestyle='solid'
    if linewidth is None:
        linewidth=2
    p = PatchCollection(patches, edgecolors=color, facecolors='none', linewidths=linewidth, linestyles=[linestyle],zorder=1)
    return p
