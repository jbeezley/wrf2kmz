#!/usr/bin/env python

import os
import subprocess as sp
import numpy as np
from netCDF4 import Dataset

verbose=True

# set to the path where gdal binaries are located
# ex: gdalPath='/usr/loca/bin'
# None means they are in the PATH environment variable
gdalPath='/usr/local/bin'
gdalEnv={}

# from WPS/geogrid/src/misc_definitions_module.F
wrfProjectionsCode={
    0:'PROJ_LATLONG',
    1:'PROJ_LC',
    2:'PROJ_PS',
    3:'PROJ_MERC',
    4:'PROJ_GAUSS',
    5:'PROJ_CYL',
    6:'PROJ_CASSINI',
    102:'PROJ_PS_WGS84',
    105:'PROJ_ALBERS_NAD83',
    203:'PROJ_ROTLL'
}

_lambert='EPSG:9802'
_polars='EPSG:9810'
_mercator='EPSG:9805'

knownProjectionsProj4={
    'PROJ_LATLONG':None,
    'PROJ_LC':'+proj=lcc +lat_1=%(truelat1)f +lat_2=%(truelat2)f +lat_0=%(truelat1)f +lon_0=%(std_lon)f +R=6370000',
    'PROJ_PS':None,  #need to check encoding... '+proj=stere +lat_ts=%(truelat1)f +lat_0=%(truelat1)f +lon_0=%(std_lon)f -k_0=1.0 +x_0=0. +y_0=0.'
    'PROJ_MERC':None #proj.4 doesn't support lat of nat. origin/wrf only gives truelat...
                     #'+proj=merc +lat_ts=%(truelat1)f +lon_0=%(std_lon)f +x_0=0. +y_0=0.'
}

_unproj='+proj=longlat +datum=WGS84 +nodefs'
_gdal_translate='gdal_translate'
_gdalwarp='gdalwarp'

_vrtFMT='''<VRTDataset rasterXSize="%(xsize)i" rasterYSize="%(ysize)i">
  <VRTRasterBand dataType="%(dtype)s" band="1" subClass="VRTRawRasterBand">
    <SourceFilename relativetoVRT="1">%(fname)s</SourceFilename>
    <PixelOffset>%(dsize)i</PixelOffset>
    <LineOffset>%(lsize)i</LineOffset>
    <ByteOrder>%(byteorder)s</ByteOrder>
  </VRTRasterBand>
</VRTDataset>
'''

def message(s):
    if verbose:
        print s

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class ProjectionError(Exception):
    pass

class UnknownProjection(ProjectionError):
    pass

class UnsupportedProjection(ProjectionError):
    pass

class NoGDALSupport(ProjectionError):
    pass

def getGDALProg(prog):
    if gdalPath is not None:
        prog=os.path.join(gdalPath,prog)
    prog=which(prog)
    if prog is None:
        raise NoGDALSupport("Could not find %s."%prog)
    return prog

def getEPSGProjectionDef(map_proj,cen_lat,cen_lon,truelat1,truelat2,
                         std_lon,pole_lat,pole_lon):

    proj=wrfProjectionsCode.get(map_proj,None)
    if proj is None:
        raise UnknownProjection
    elif proj == 'PROJ_LATLONG':
        return None
    
    proj=knownProjectionsProj4.get(proj,None)
    if proj is None:
        raise UnsupportedProjection

    d={'cen_lat':cen_lat,'cen_lon':cen_lon,'truelat1':truelat1,'truelat2':truelat2,\
       'std_lon':std_lon,'pole_lat':pole_lat,'pole_lon':pole_lon}

    return proj % d

def georeferenceImage(img,proj,gcps,out):
    gt=getGDALProg(_gdal_translate)
    args=[gt,'-of','VRT','-a_srs',proj]
    for g in gcps:
        args.append('-gcp')
        args.append('%(ix)i' % g)
        args.append('%(iy)i' % g)
        args.append('%(lon)f' % g)
        args.append('%(lat)f' % g)
    args.append(img)
    args.append(out)
    message('running: %s' % (' '.join(args)))
    p=sp.Popen(args,shell=False,env=gdalEnv,stdout=sp.PIPE)
    p.communicate()
    if p.returncode != 0:
        print p.returncode
        raise NoGDALSupport('Failed to execute %s' % _gdal_translate)

def warpImage(img,out):
    if os.path.exists(out):
        os.remove(out)
    gw=getGDALProg(_gdalwarp)
    args=[gw,'-of','netCDF','-t_srs',_unproj,img,out]
    message('running: %s' % (' '.join(args)))
    p=sp.Popen(args,shell=False,env=gdalEnv,stdout=sp.PIPE)
    p.communicate()
    if p.returncode != 0:
        raise NoGDALSupport('Failed to execute %s' % _gdalwarp)

def createGCP(ix,iy,lon,lat):
    return {'ix':ix,'iy':iy,'lon':lon,'lat':lat}

def unProjectImage(img,proj,gcps):
    vrt=img+'.vrt'
    tmp='tmp_'+vrt

    georeferenceImage(img,proj,gcps,vrt)
    warpImage(vrt,tmp)

def vrtFromArray(fname,a):
    assert a.ndim == 2
    
    a=np.ascontiguousarray(a)
    
    bin=fname+'.bin'
    b=a.astype('<f4')
    b.tofile(bin)

    ny,nx=a.shape

    s=_vrtFMT % {'xsize':nx,'ysize':ny,'dtype':'Float32',
                 'fname':bin,'dsize':4,'lsize':4*nx,'byteorder':'LSB'}
    open(fname,'w').write(s)

def readNC(fname):
    f=Dataset(fname,'r')
    bds={'west':f.variables['lon'][0],
         'east':f.variables['lon'][-1],
         'south':f.variables['lat'][0],
         'north':f.variables['lat'][-1]}
    return np.flipud(f.variables['Band1'][:]),bds

# find gdal programs or raise error on import
getGDALProg(_gdal_translate)
getGDALProg(_gdalwarp)

if __name__ == '__main__':
    from netCDF4 import Dataset
    import sys
    (nc,img)=sys.argv[1:]

    f=Dataset(nc,'r')
    srx,sry=10,10
    arr=f.variables['NFUEL_CAT'][-1:,:-sry,:-srx].squeeze()
    lat=f.variables['FXLAT'][0,:-sry,:-srx]
    lon=f.variables['FXLONG'][0,:-sry,:-srx]
    
    ny,nx=arr.shape
    ny=ny-1
    nx=nx-1
    gcp=[createGCP(0,0,lon[0,0],lat[0,0]),
         createGCP(0,ny,lon[ny,0],lat[ny,0]),
         createGCP(nx,0,lon[0,nx],lat[0,nx]),
         createGCP(nx,ny,lon[ny,nx],lat[ny,nx])]
    proj=getEPSGProjectionDef(f.MAP_PROJ,f.CEN_LAT,f.CEN_LON,
                              f.TRUELAT1,f.TRUELAT2,f.STAND_LON,
                              f.POLE_LAT,f.POLE_LON)
    
    vrtFromArray('test.vrt',np.flipud(arr))
    georeferenceImage('test.vrt',proj,gcp,'tmp_test.vrt')
    warpImage('tmp_test.vrt','warp_test.nc')
    a=readNC('warp_test.nc')
