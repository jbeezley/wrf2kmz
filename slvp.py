import numpy as np
import slvpf
compute_seaprs_from_derived=slvpf.slvp.compute_seaprs_from_derived

def compute_t(t,p):
    return (t+300.)*(p/100000.)**0.287

def compute_height(ph):
    z=ph/9.81
    return (z[:-1,:,:] + z[1:,:,:])/2.

def compute_pressure(ph,phb):
    return ph+phb

def compute_seaprs(p,pb,ph,phb,t,qv):
    ph=compute_pressure(ph,phb)
    p=compute_pressure(p,pb)
    z=compute_height(ph)
    t=compute_t(t,p)

    slvp,mask=compute_seaprs_from_derived(z.T,t.T,p.T,qv.T)
    slvp=np.ma.array(slvp.T)
    slvp.mask=(mask.T != 0)
    slvp.set_fill_value(np.nan)
    return slvp

def compute_seaprs_from_derived_py(z,t,p,q):

    R=287.04
    G=9.81
    GAMMA=0.0065
    TC=273.16+17.5
    PCONST=10000

    nz,ny,nx=z.shape
    bad=np.zeros((ny,nx))
    level=-1*np.ones((ny,nx))
    t_surf=np.ones((ny,nx)) # avoid divide by zero on masked elements
    t_sea_level=np.zeros((ny,nx))
    sea_level_pressure=np.zeros((ny,nx))

    for j in xrange(ny):
        for i in xrange(nx):
            for k in xrange(nz):
                if p[k,j,i] < p[0,j,i] - PCONST:
                    level[j,i] = k
                    break

    bad[level < 0]=1

    if np.any(bad):
        print 'Trouble finding level above ground... masking output array'
    
    CEXP=GAMMA*R/G
    for j in xrange(ny):
        for i in xrange(nx):
            if not bad[j,i]:
                klo=max(level[j,i]-2,0)
                khi=min(klo+1,nz-2)

                if klo == khi:
                    raise Exception("Trapping levels are weird")
                
                plo=p[klo,j,i]
                phi=p[khi,j,i]
                tlo=t[klo,j,i] * (1.+0.608 * q[klo,j,i])
                thi=t[khi,j,i] * (1.+0.608 * q[khi,j,i])

                zlo=z[klo,j,i]
                zhi=z[khi,j,i]

                p_at_pconst=p[0,j,i] - PCONST
                t_at_pconst=thi-(thi-tlo)*np.log(p_at_pconst/phi)*np.log(plo/phi)
                z_at_pconst=zhi-(zhi-zlo)*np.log(p_at_pconst/phi)*np.log(plo/phi)

                t_surf[j,i]=t_at_pconst*(p[0,j,i]/p_at_pconst)**CEXP
                t_sea_level[j,i]=t_at_pconst+GAMMA*z_at_pconst

                sea_level_pressure[j,i]=p[0,j,i] * np.exp(2.*G*z[0,j,i]/ (R*(t_sea_level[j,i]+t_surf[j,i])))
    sea_level_pressure=np.ma.array(sea_level_pressure)
    sea_level_pressure.mask=(bad != 0)
    sea_level_pressure.set_fill_value(np.nan)

    return sea_level_pressure,bad

if __name__ == "__main__":
    import sys
    from netCDF4 import Dataset
    fname=sys.argv[1]
    f=Dataset(fname,'r')

    p=f.variables['P'][0,:]
    pb=f.variables['PB'][0,:]
    ph=f.variables['PH'][0,:]
    phb=f.variables['PHB'][0,:]
    t=f.variables['T'][0,:]
    qv=f.variables['QVAPOR'][0,:]

    slvp=compute_seaprs(p,pb,ph,phb,t,qv)

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.pylab import *
    
    print slvp.min(),slvp.max()
    imshow(np.flipud(slvp),interpolation='nearest')
    colorbar()
    savefig('slvp.png')

