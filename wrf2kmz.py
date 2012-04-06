#!/usr/bin/env python

'''
wrf2kmz
Jonathan Beezley <jon.beezley@gmail.com>

This is a module/commandline script that will create kmz visualizations
from wrf netCDF files.  The module contains three main classes:

    BaseNetCDF2Raster:
        Generates georeferenced rasters from netcdf files.  It is designed
        to be overloaded for different model output conventions.  Most
        of the code should be generic enough to work with any netCDF
        output file, with the exception of the georeferencing and 
        time coding.  These functions need to be overloaded as in the
        subclass FireNetcdf2Raster that is designed to work with WRF-Fire.
        This class also contains variables that control how the rasters
        will be displayed in the kmz.

    FireRasterFile:
        Purely a convenience class used for storing defaults for various
        variables that are present in WRF-Fire output files.

    ncKML:
        A subclass of simplekml.Kml that has additional functions taking
        subclasses of BaseNetCDF2Raster to construct kml elements such
        as ground overlays.

When used as a commandline script as `python wrf2kmz.py wrfout [var1 [var2 [...]]]`
a kmz is produced containing a visualization of the WRF-Fire output file wrfout.
By default, this kmz will contain an animation of the fire perimeter.  With additional
variable names specified on the commandline, it will add ground overlays to this file.
At this time, only 2D surface variables are supported.

Additional caveat:  kmz files expect that ground overlays are unprojected, but
WRF models can output on a variety of projections.  At this time, this file will
not unproject data for ground overlays.  So, if your WRF simulation occurs using 
anything other than regular_ll, the georeferencing will be incorrect.  For the small
domains (order 10 km) typically used in fire simulations, this error is typically
quite small (less than the mesh resolution).  In the future, I plan to add an optional
image reprojecting component based off of gdal python bindings.

Dependencies:
    
    simplekml      : http://code.google.com/p/simplekml/
    matplotlib     : http://matplotlib.sourceforge.net/
    netcdf4-python : http://code.google.com/p/netcdf4-python/
'''

# standard library imports
from dateutil import parser
from cStringIO import StringIO

# dependency imports
import matplotlib

# set backend to Agg when using as a commandline script
# so it will work without a X11
if __name__ == '__main__':
    matplotlib.use('Agg')

from matplotlib import pylab
import numpy as np
from netCDF4 import Dataset

try:
    from simplekml import *
except ImportError:
    import sys
    print 'Could not find simplekml module.'
    print 'If you have setuptools installed, try `easy_install simplekml`.'
    print 'Otherwise install from http://code.google.com/p/simplekml/'
    sys.exit(1)

# verbose=False to be quiet
verbose=True
def message(s):
    '''
    Print a message to stdout if verbose is True.
    '''
    if verbose:
        print s

class MaskedArrayException(Exception):
    '''
    An exception that is raised when a variable read contains only
    masked (or invalid) data.  This allows methods that loop over all
    time steps in a file to skip those that would otherwise contain
    empty or completely transparent images.
    '''
    pass

class BaseNetCDF2Raster(object):
    '''
    Base class for creating georeferenced raster images from netCDF data.
    Subclasses must at lest override getStepTime.
    '''

    # default matplotlib norm object (maps variable range to [0,1] for coloring)
    defaultNorm=matplotlib.colors.Normalize

    # default matplotlib colormap
    defaultcmap=pylab.cm.jet

    # default colorbar label formatter
    defaultFormatter=None

    # when true use global variable range for color axis range
    # rather than different color scales for each time step
    # unimplemented!
    defaultminmaxglobal=False

    # A list containing values that are ignored in the input data.
    # These values will display as transparent in the output raster
    # and will not be used for computation of the color axis.
    defaultmaskedValues=None

    # mask all values > this value
    defaultmaskedAbove=None

    # mask all value < this value
    defaultmaskedBelow=None

    def __init__(self,file,var,norm=None,cmap=None,formatter=None,\
                 minmaxglobal=None,maskedValues=None,maskedAbove=None, \
                 maskedBelow=None,static=False,displayName=None,
                 displayDescription=None,displayColorbar=True,
                 displayAlpha=180,name=None):
        '''
        Initialize a raster class object.

        Required arguments:

            file:       opened netCDF4 Dataset or MFDataset object
            var:        netCDF4 Variable object from file
                        (the variable to visualize)

        Optional (keyword) arguments:

            norm:           matplotlib color norm (class default: defaultNorm)
            cmap:           matplotlib color map (class default: defaultcmap)
            formatter:      matplotlib label formatter 
                            (class default: defaultFormatter)
            maskedValues:   list of values to be masked in output 
                            (class default: defaultMaskedValues)
            maskedBelow:    Values less than this will be masked.
                            (class default: defaultmaskedBelow)
            maskedAbove:    Values greater than this will be masked.
                            (class default: defaultmaskedAbove)
            static:         If true, assume the variable doesn't change in each 
                            timestep.  Only one image will be created with a
                            time range from the begining of the file to the end.
                            (default False)
            displayName:    This is how the variable will be called in the kmz file. 
                            (default var._name)
            displayDescription: 
                            This is a description that will be added to the variable in 
                            kmz.  (default var.description if present)
            displayColorbar:If true, a colorbar will be created with the ground overlay.
                            (default True)
            displayAlpha:   An integer in the range 0-255 giving the transparency of the
                            ground overlay.  0 means transparent, 255 means opaque.
                            (default 180)
            minmaxglobal: unimplemented

        '''

        # make sure the variable is 2D
        ndim=var.ndim
        if file.dimensions[var.dimensions[0]].isunlimited():
            ndim=ndim-1
        if ndim != 2:
            raise Exception("Only 2D variables are supported.")

        # store arguments
        self._file=file
        self._var=var
        self._minmax=None
        self._norm=self.setToDefaultifNone(norm,self.defaultNorm)
        self._cmap=self.setToDefaultifNone(cmap,self.defaultcmap)
        self._formatter=self.setToDefaultifNone(formatter,self.defaultFormatter)
        self._minmaxglobal=self.setToDefaultifNone(minmaxglobal,self.defaultminmaxglobal)
        self._maskedValues=self.setToDefaultifNone(maskedValues,self.defaultmaskedValues)
        self._maskedAbove=self.setToDefaultifNone(maskedAbove,self.defaultmaskedAbove)
        self._maskedBelow=self.setToDefaultifNone(maskedBelow,self.defaultmaskedBelow)
        self._static=static
        self._tdim=self.getTimeDim()
        self.displayName=displayName
        self.displayDescription=displayDescription
        self.displayColorbar=displayColorbar
        self.displayAlpha=displayAlpha
        self._name=name

        if self._minmaxglobal:
            raise Exception("Global min-max computation not yet supported.")
        
        # get the number of time steps in this file from the unlimited dimension
        # of the netcdf file
        if self._tdim is None:
            self._nstep=1
        else:
            self._nstep=len(self._file.dimensions[self._tdim])

    # The following are some property getters, in case of later abstractions.
    @property
    def norm(self):
        return str(self._norm).split('.')[-1:]

    @property
    def cmap(self):
        return self._cmap.name

    @property
    def minmaxgloba(self):
        return self._minmaxglobal

    @property
    def maskedValues(self):
        return self._maskedValues

    @property
    def maskedAbove(self):
        return self._maskedAbove

    @property
    def maskedBelow(self):
        return self._maskedBelow

    @property
    def static(self):
        return self._static
    
    @staticmethod
    def setToDefaultifNone(value,default):
        '''
        Return value or default if value is None.
        '''
        if value is None:
            target=default
        else:
            target=value
        return target

    def _readArray(self,istep=None):
        '''
        Read data from the netCDF variable.  If istep is not None, then
        read a single time step.
        '''

        # read data
        if istep is None:
            a=self._var[:]
        else:
            a=self._var[istep,...].squeeze()

        # convert to a masked array
        a=np.ma.MaskedArray(a,copy=False)

        # mask according to instance properties
        if self._maskedValues is not None:
            for m in self.maskedValues:
                a=np.ma.masked_equal(a,m,copy=False)
        if self._maskedAbove is not None:
            a=np.ma.masked_greater(a,self.maskedAbove,copy=False)
        if self._maskedBelow is not None:
            a=np.ma.masked_less(a,self.maskedBelow,copy=False)

        # fill masked values with NaN's to display as transparent in matplotlib
        a.fill_value=np.nan

        # raise an exception if all elements of the array are masked
        if a.mask.all() or (a.min() == a.max()):
            raise MaskedArrayException
        return a

    def readCoordinates(self,istep=0,idx=None):
        '''
        Read coordinate arrays (as lon/lat) for this variable.  The class assumes
        that the coordinate arrays have the same shape as the variable array.

        istep: The time step to read from the coordinate arrays.
        idx:   If we are displaying as subarray from the variable, this
               contains the index bounds that we are displaying as returned
               from _getRestriction.
        '''

        # get coordinate array names
        c=self.getCoordinates()
        if c is None:
            raise Exception("Could not find coordinate array for %s" % self.getName())

        # read coordinate arrays from the file
        lon=self._file.variables[c[0]]
        lat=self._file.variables[c[1]]

        # if no index bounds given, then return the full arrays
        if idx is None:
            idx=(0,lon.shape[0]-1,0,lon.shape[1])

        # restrict coordinate arrays from index bounds
        if istep is None:
            lon=lon[idx[0]:idx[1]+1,idx[2]:idx[3]+1]
            lat=lat[idx[0]:idx[1]+1,idx[2]:idx[3]+1]
        else:
            lon=lon[istep,idx[0]:idx[1]+1,idx[2]:idx[3]+1].squeeze()
            lat=lat[istep,idx[0]:idx[1]+1,idx[2]:idx[3]+1].squeeze()
        return (lon,lat)
            
    @classmethod
    def _getRestriction(cls,a):
        '''
        Return indices of the smallest subarray covering the non-masked values of a.
        Assumes a is 2D.
        '''
        assert a.ndim == 2

        # get non-masked indices of a
        idx=(a.mask == False).nonzero()

        if len(idx) == 1 or len(idx[0]) == 0:
            # if no indices are masked return full array indices
            return (0,a.shape[0]-1,0,a.shape[1]-1)
        else:
            z=idx

            # get the smallest and largest nonmasked index in each coordinate
            imin=z[0].min()
            imax=z[0].max()
            jmin=z[1].min()
            jmax=z[1].max()

            # make sure the subarray is at least 2x2, otherwise
            # things start to break
            if imax == imin:
                if imax < a.shape[0]-1:
                    imax=imax+1
                else:
                    imin=imin-1
            if jmax == jmin:
                if jmax < a.shape[1]-1:
                    jmax=jmax+1
                else:
                    jmin=jmin-1
            z=(imin,imax,jmin,jmax)
            return z

    def getMinMax(self,istep=None):
        '''
        Returns the minimum and maximum of the variable at istep.
        '''
        if self._minmaxglobal:
            if self._minmax is None:
                a=self._readArray()
                self._minmax=(a.min(),a.max())
            minmax=self._minmax
        else:
            a=self._readArray(istep)
            minmax=(a.min(),a.max())
        return minmax

    def getUnits(self):
        '''
        Returns the units of the variable or an empty string if the variable
        has no "units" attribute.
        '''
        return self._var.__dict__.get('units','')

    def getName(self):
        '''
        Returns the name of the variable in the netCDF file.
        '''
        if self._name is None:
            return self._var._name
        else:
            return self._name

    def getDescription(self):
        '''
        Returns a description of the variable from the attribute "description"
        or None if no description is present.
        '''
        return self._var.__dict__.get('description',None)

    def getCoordinates(self):
        '''
        Returns a tuple containing the coordinate array names for the variable.
        Following WRF output standards, this is in the attribute coordinates
        containing a 2 element tuple (longitude,latitude).
        '''
        return self._var.__dict__.get('coordinates',None)

    def georeference(self,istep=None):
        '''
        Get a georeference for the variable at istep.  The georeference returned
        is a dictionary containing the west, east, south, and north bounds of the
        variable and assumes unprojected data.  This method should return the 
        coordinates of the subarray computed from the variable mask rather than
        the full variable bounds.
        '''
        a=self._readArray(istep)
        idx=self._getRestriction(a)
        lon,lat=self.readCoordinates(istep,idx)
        return {'west':lon.min(),'east':lon.max(),'south':lat.min(),'north':lat.max()}
   
    def getTimeDim(self):
        '''
        Returns the name of the unlimited (or time) dimension of the netcdf file or
        None if no unlimited dimension is used.
        '''
        for d in self._file.dimensions:
            if self._file.dimensions[d].isunlimited():
                return d
        return None

    def timereference(self,istep=None):
        '''
        Returns the time span for the variable at istep.  By default, this is 
        time span begins at the time of the current time step to the time of the
        next time step.  The time span of the last time step is computed from 
        delta time from the last two time steps in the file.

        Returns a dictionary with keywords 'start' and 'end' and values
        datetime objects.
        '''

        # if no istep given return the full file time span
        if istep is None:
            sstep=0
            estep=self._nstep
        else:
            sstep=istep
            estep=istep+1

        # get the time at the start
        start=self.getStepTime(sstep)

        # if the ending time step is valid get it, 
        # otherwise estimate the ending from the delta
        # time of the last two steps
        if estep >= self._nstep:
            end=self.getStepTime(self._nstep-1)
            dt=end-self.getStepTime(self._nstep-2)
            end=end+dt
        else:
            end=self.getStepTime(estep)

        return {'start':start,'end':end}

    def getStepTime(self,istep):
        '''
        Must be implemented by a subclasss.  Returns a datetime object 
        representing the time at istep in the output file.
        '''
        raise Exception('Unimplemented base class method')

    def getRasterFromArray(self,a,hsize=3,dpi=300):
        '''
        Returns a string containing a png psuedocolor image of the array a.

        Keyword arguments:

            hsize:  integer, width of the image in inches
            dpi:    integer, dots per inch of the image
        '''

        # get subarray restriction from mask
        idx=self._getRestriction(a)

        # restrict the array
        a=a[idx[0]:idx[1]+1,idx[2]:idx[3]+1]

        # generate a matplotlib figure object
        fig=pylab.figure(figsize=(hsize,hsize*float(a.shape[0])/a.shape[1]))
        ax=fig.add_axes([0,0,1,1])

        # get a color norm instance from the min/max of a
        norm=self._norm(a.min(),a.max())

        # add image to the axis
        ax.imshow(np.flipud(a),cmap=self._cmap,norm=norm,interpolation='nearest')

        # turn the axis off to get a bare raster
        ax.axis('off')

        # a file like string object to save the png to
        im=StringIO()

        # save and close the figure
        fig.savefig(im,dpi=dpi,format='png',transparent=True)
        pylab.close(fig)

        return im.getvalue()

    def getRaster(self,istep=None,**kwargs):
        '''
        Get a png image for the variable at istep.  See getRasterFromArray.
        '''
        a=self._readArray(istep)
        return self.getRasterFromArray(a,**kwargs)

    def getColorbarFromMinMax(self,min,max,hsize=2,dpi=200):
        '''
        Returns a colorbar as a string with given min/max values.

        Keyword arguments:
            hsize:  integer, horizontal size of the image in inches
            dpi:    integer, dots per inch of the image
        '''

        # get the color norm from the min/max
        norm=self._norm(min,max)
        
        # construct keyword arguments for ColorbarBase constructor
        kwargs={'norm':norm,'spacing':'proportional','orientation':'vertical',
                'cmap':self._cmap}
        if self._formatter is not None:
            kwargs['format']=self._formatter
        
        # create a matplotlib figure
        fig=pylab.figure(figsize=(hsize,hsize*5./2.))
        ax=fig.add_axes([.25,.03,.1,.9])

        # create Colorbar instance
        cb=matplotlib.colorbar.ColorbarBase(ax,**kwargs)

        # set the label as the units of the variable
        cb.set_label(self.getUnits(),color='1')

        # set the title as the name of the variable
        ax.set_title(self.getName(),color='1')
        for tl in ax.get_yticklabels():
            tl.set_color('1')
        
        # save png to a file-like string
        im=StringIO()
        fig.savefig(im,dpi=dpi,format='png',transparent=True)
        return im.getvalue()

    def getColorbar(self,istep=None,**kwargs):
        '''
        Returns a colorbar for the variable at istep.
        See getColorbarFromMinMax
        '''
        min,max=self.getMinMax(istep)
        return self.getColorbarFromMinMax(min,max,**kwargs)
    
    def perimeterFromContour(self,contour=0,istep=None):
        '''
        Returns a contour of the variable as a list of vertices
        of a polygon.  Generates the contour from a matplotlib
        contour plot. 
        '''
        
        # read the array
        a=self._readArray(istep)
        assert len(a.shape) == 2

        # get coordinates of the array
        lon,lat=self.readCoordinates(istep)

        # get vertices of the polygon
        c=pylab.contour(lon,lat,a,[contour]).collections[0]
        p=c.get_paths()
        poly=[]
        for i in p:
            tmp=[ (x,y) for x,y in i.to_polygons()[0] ]
            poly.append(tmp)
        return poly

class WRFNetcdf2Raster(BaseNetCDF2Raster):
    '''
    This is an implementation of BaseNetCDF2Raster for WRF output files.
    
    See BaseNetCDF2Raster docstrings for details.
    '''
    def __init__(self,*args,**kwargs):
        super(WRFNetcdf2Raster,self).__init__(*args,**kwargs)
    
    def getStepTime(self,istep):
        '''
        Returns a datetime object representing the time for time step "istep" from 
        the output file.
        '''
        time=parser.parse(self._file.variables['Times'][istep,:].tostring().replace('_',' '))
        return time

class FireNetcdf2Raster(WRFNetcdf2Raster):
    '''
    BaseNetcdf2Raster implementation for WRF-Fire files.  Generally this is the same as the
    super class WRFNetcdf2Raster, but contains additional code for working with subgrid 
    variables.  General WRF files/variables should work with this class as well.
    '''

    def __init__(self,*args,**kwargs):
        super(FireNetcdf2Raster,self).__init__(*args,**kwargs)

        # In addition to standard properties we add srx/sry, the subgrid refinement ratios
        # to restrict fire grid variables to their proper size.
        self.srx=0
        self.sry=0
        dims=self._file.dimensions
        if self._isfiregridvar(self._var):
            self.srx=len(dims['west_east_subgrid'])/(len(dims['west_east'])+1)
            self.sry=len(dims['south_north_subgrid'])/(len(dims['south_north'])+1)

    @classmethod
    def _isfiregridvar(cls,var):
        '''
        Returns true if the variable is a fire (or subgrid) variable.
        '''
        return var.dimensions[-1][-8:] == '_subgrid'

    def _readArray(self,istep=0):
        '''
        Same as base class _readArray, but ignores the extra space that is present
        for fire grid variables.
        '''
        a=WRFNetcdf2Raster._readArray(self,istep)
        return a[:-self.sry,:-self.srx]

    def readCoordinates(self,istep=0,idx=None):
        '''
        Same as base class readCoordinates, but ignores the extra space that is present
        for fire grid variables.
        '''
        if idx is None:
            c=self._var.shape
            idx=(0,c[1]-self.sry-1,0,c[2]-self.srx-1)
        return WRFNetcdf2Raster.readCoordinates(self,istep,idx)

    def getCoordinates(self):
        '''
        Same as base class getCoordinates, but returns the actual fire grid coordinates
        instead of XLONG/XLAT because WRF fills in the attributes incorrectly.
        '''
        if not self._isfiregridvar(self._var):
            return WRFNetcdf2Raster.getCoordinates()
        else:
            return ('FXLONG','FXLAT')

class LogScaleRaster(FireNetcdf2Raster):
    '''
    A raster class with defaults that generate log scale images.
    '''
    defaultNorm=matplotlib.colors.LogNorm
    defaultFormatter=matplotlib.ticker.LogFormatter(10,labelOnlyBase=False)
    defaultmaskedValues=[0.]
    defaultmaskedBelow=0.

    def __init__(self,*args,**kwargs):
        super(LogScaleRaster,self).__init__(*args,**kwargs)

class ZeroMaskedRaster(FireNetcdf2Raster):
    '''
    A raster class that masks out values == 0.
    '''
    defaultmaskedValues=[0.]
    
    def __init__(self,*args,**kwargs):
        super(ZeroMaskedRaster,self).__init__(*args,**kwargs)

class NegativeMaskedRaster(FireNetcdf2Raster):
    '''
    A raster class that masks out values < 0.
    '''
    defaultmaskedBelow=0.
    
    def __init__(self,*args,**kwargs):
        super(NegativeMaskedRaster,self).__init__(*args,**kwargs)

class FirePerimeter(FireNetcdf2Raster):
    '''
    A raster class designed to output fire perimeters.
    '''
    def __init__(self,*args,**kwargs):
        if len(args) == 1:
            args= args+(args[0].variables['LFN'],)
        if not kwargs.has_key('name'):
            kwargs['name']='fire perimeter'
        super(FirePerimeter,self).__init__(*args,**kwargs)

    def getName(self):
        return 'Fire perimeter'
    
class FireRasterFile(object):
    '''
    This is a convenience class for storing defaults for various variables
    that are present in WRF-Fire output files.
    '''

    # define display styles for individual variables
    _varClasses={
        'FGRNHFX':(LogScaleRaster,{}),
        'GRNHFX':(LogScaleRaster,{}),
        'FLINEINT':(NegativeMaskedRaster,{}),
        'FLINEINT2':(NegativeMaskedRaster,{}),
        'NFUEL_CAT':(ZeroMaskedRaster,{'static':True})
    }

    # default display style for variables not listed above
    _defaultClass=(ZeroMaskedRaster,{})
    
    # fire perimeter class
    _varPerimeter=(FirePerimeter,{})

    def __init__(self,filename,globalminmax=False):
        '''
        Initialize a FireRasterFile object. 

        filename:  string, the wrfout file to open
        '''
        self._file=Dataset(filename,'r')
        self._globalminmax=globalminmax

    def rasterFromVar(self,varname):
        '''
        Generate a sequence of images from a variable in the wrfout.
        This is mostly just here for testing, but could be useful
        somehow.
        '''

        # get a raster class instance
        v=self.rasterClassFromVar(varname)

        # generate a sequence of png images
        for istep in xrange(v._nstep):
            try:
                s=v.getRaster(istep)
            except MaskedArrayException:
                print 'skipping %s at step=%i' % (varname,istep)
                continue
            open('%s_%04i.png' % (varname,istep),'w').write(s)

            s=v.getColorbar(istep)
            open('%s_%04i_c.png' % (varname,istep),'w').write(s)

    def rasterClassFromVar(self,varname):
        '''
        Get a raster class instance from the given variable with 
        default style.
        '''
        vclass,vargs=self._varClasses.get(varname,self._defaultClass)
        if not vargs.has_key('name'):
            vargs['name']=varname
        return vclass(self._file,self._file.variables[varname],**vargs)
    
    def firePerimeterClass(self):
        '''
        Get a fire perimeter class.
        '''
        vclass,vargs=self._varPerimeter
        return vclass(self._file,**vargs)

class ncKML(Kml):
    '''
    A subclass of simplekml.Kml that works with subclasses of BaseNetCDF2Raster to 
    generate kml elements such as ground overlays easily.
    '''

    # number of meters per degree, used to calculate initial viewing
    # range so that the whole domain will appear in the view.
    _mperdegree=111325.

    # increase the view range by this factor
    _lookatexpand=1.5

    # where to store files that go into the kmz
    _filesdir='files'

    def __init__(self,kmlname='ncKML',description=None):
        '''
        Initialize an ncKML object with suitable defaults.
        '''

        # parent constructor
        super(ncKML,self).__init__(name=kmlname,description=None)

        # make the root element open by default
        self._root=self.document
        self._root.open=1

        # make a directory to save images to, but don't bomb
        # out if it already exists
        try:
            os.makedirs(self._filesdir)
        except Exception:
            pass

    def set_initial_view(self,georeference,time):
        '''
        Set the initial view after opening the kmz file.  This is
        useful because the time slider will be closed at the 
        initial time rather than displaying all time steps.
        The camera will start centered at the georeference at a
        range calculated to show the whole domain.

        Arguments:
            georeference:   dictionary like that returned from BaseNetCDF2Raster.georeference
            time:           datetime object that the time slider will start at 
        '''

        # time stamp of initial view
        t=GxTimeStamp()
        t.when=time.isoformat()
        
        # get height above the ground necessary to show the whole domain
        lon=(georeference['west']+georeference['east'])/2.
        lat=(georeference['south']+georeference['north'])/2.
        dx=(georeference['east']-georeference['west'])*self._mperdegree
        dy=(georeference['north']-georeference['south'])*self._mperdegree
        rng=self._lookatexpand*max(dx,dy)*(3**.5)/2.

        # create camera LookAt object and add it to the document
        l=LookAt(latitude=lat,longitude=lon,range=rng,gxtimestamp=t)
        self.document.lookat=l
        return l
    
    def setViewFromRaster(self,raster):
        '''
        Generates an initial view of the kmz from the given raster object.
        It calculates georeferencing and the start time from the raster object.
        '''

        # try and find a non-masked time step to show
        for i in xrange(raster._nstep):
            try:
                gref=raster.georeference(i)
                tref=raster.getStepTime(i)
            except MaskedArrayException:
                continue
            break
        else:
            raise Exception("Could not find a valid time step in variable %s" % raster.getName())
        return self.set_initial_view(gref,tref)

    @staticmethod
    def _normalizeIsteps(raster,isteps):
        '''
        Returns an iterable of all time steps in the file or [0]
        if the raster is marked as static.
        '''
        if raster.static:
            isteps=[0]
        if isteps is None:
            isteps=xrange(raster._nstep)

        if isinstance(isteps,int):
            isteps=[isteps]
        return isteps

    def groundOverlayFromRaster(self,raster, \
                                isteps=None,visible=False,           \
                                fnamefmt='%(dir)s/%(name)s_%(step)05i.png',       \
                                calpha=256):
        '''
        Add a sequence of ground overlays to the kmz from a raster object.

            raster: An instance of a subclass of BaseNetCDF2Raster

        Optional:
            isteps:     A list of time steps to use or None for all
            visible:    If the overlays should be visible by default.
            fnamefmt:   A format string that determines how to save
                        the png in the kmz.  Really only useful if the
                        same variable is displayed more than once with
                        different styles, otherwise the default is fine.
            calpha:     alpha value for the colorbar 0-255.
        '''

        # get style properties from the raster object
        name=raster.displayName
        description=raster.displayDescription
        alpha=raster.displayAlpha
        colorbar=raster.displayColorbar
        
        # This can take a while so tell the use what is going on.
        message('Creating ground overlay from %s' % raster.getName())

        if name is None:
            name=raster.getName()
        if description is None:
            description=raster.getDescription()

        # generate a new folder element to organize the images
        f=self._root.newfolder(name=name,description=description)
        if visible:
            f.visibility=1
        else:
            f.visibility=0

        # keep the folder closed
        f.open=False
        
        # get a list of steps to loop over (static variables only get generated once) 
        isteps=self._normalizeIsteps(raster,isteps)
        for i in isteps:
            if raster.static:
                tref=raster.timereference()
            else:
                tref=raster.timereference(i)
            
            try:
                # skipped time steps where all of the data is masked
                gref=raster.georeference(i)
                img=raster.getRaster(i)
            except MaskedArrayException:
                message('Skipping %i'%i)
                continue

            # get a filename to save the image to 
            fname=fnamefmt % {'dir':self._filesdir,'name':raster.getName(),'step':i}
            open(fname,'w').write(img)

            # add a ground overlay element and populate its properties
            g=f.newgroundoverlay(name='%s_%05i' % (raster.getName(),i) )
            g.latlonbox.north=gref['north']
            g.latlonbox.south=gref['south']
            g.latlonbox.east=gref['east']
            g.latlonbox.west=gref['west']
            g.color=Color.rgb(255,255,255,a=alpha)
            g.altitudemode=AltitudeMode.clamptoground
            g.timespan=TimeSpan()
            g.timespan.begin=tref['start'].isoformat()
            g.timespan.end=tref['end'].isoformat()
            g.icon.href=fname
            g.visibility=f.visibility
            
            if colorbar:
                # generate a colorbar
                img=raster.getColorbar(i)

                # get file name to save the colorbar to
                fname=fnamefmt % {'dir':self._filesdir,'name':raster.getName()+'_c','step':i}
                open(fname,'w').write(img)

                # add a new screen overlay element and populate its properties
                g=f.newscreenoverlay(name='%s_c_%05i' % (raster.getName(),i) )
                g.overlayxy=OverlayXY(x=.15,y=.5,xunits=Units.fraction,yunits=Units.fraction)
                g.screenxy=ScreenXY(x=0,y=.5,xunits=Units.fraction,yunits=Units.fraction)
                g.size.x=0
                g.size.y=.75
                g.size.xunits=Units.fraction
                g.size.yunits=Units.fraction
                g.timespan=TimeSpan()
                g.timespan.begin=tref['start'].isoformat()
                g.timespan.end=tref['end'].isoformat()
                g.color=Color.rgb(255,255,255,a=calpha)
                g.visibility=f.visibility
                g.icon.href=fname
        return f

    def polygonFromContour(self,raster,name=None,description=None, \
                           contour=0,isteps=None,     \
                           visible=False,polystyle=None,linestyle=None):
        '''
        Add a sequence of polygons to the kmz from contours of a raster object.

        raster:     An instance of a subclass of BaseNetCDF2Rasterff

        Optional:
            isteps:     A list of time steps to use or None for all
            visible:    If the overlays should be visible by default.
            polystyle:  A PolyStyle instance defining how the polygon is drawn.
            linestyle:  A LineStyle instance defining how the polygon boundary is drawn.

        The default style is a red outline of width 3.  See simplekml documentation for
        defining new styles.
        '''

        # tell the user what is going on
        message('Creating polygon from %f contour of %s' % (contour,raster.getName()))

        # get display names from raster properties
        if name is None:
            name=raster.getName()# + " contour at %i" % int(contour)

        # get a line of time steps to generate
        isteps=self._normalizeIsteps(raster,isteps)

        # create a new folder to store the perimeters in 
        f=self._root.newfolder(name=name,description=description)
        if visible:
            f.visibility=1
        else:
            f.visibility=0

        # make the folder closed
        f.open=False

        # generate default styles
        if polystyle is None:
            polystyle=PolyStyle(fill=0,outline=1,color=Color.rgb(255,0,0,255))
        if linestyle is None:
            linestyle=LineStyle(color=Color.rgb(255,0,0,255),width=3)

        for i in isteps:
            if raster.static:
                tref=raster.timereference()
            else:
                tref=raster.timereference(i)
            
            try:
                poly=raster.perimeterFromContour(contour,i)
            except MaskedArrayException:
                message('Skipping %i'%i)
                continue
            
            j=0
            for pl in poly:
                # add new polygon element to the kml class object
                p=f.newpolygon(name='%s_%05i_%05i' % (raster.getName()+'_contour',i,j))
                p.outerboundaryis=pl
                p.polystyle=polystyle
                p.linestyle=linestyle
                p.tessellate=1
                p.timespan=TimeSpan()
                p.timespan.begin=tref['start'].isoformat()
                p.timespan.end=tref['end'].isoformat()
                p.visibility=f.visibility
                j=j+1
        return f

def test(wrfout):
    '''
    Testing function/usage examples.
    '''
    if False:
        f=Dataset(wrfout,'r')
        r1=FireNetcdf2Raster(f,f.variables['UF'],name='UF')
        r2=FireNetcdf2Raster(f,f.variables['FGRNHFX'],name='FGRNHFX')
        r3=FireNetcdf2Raster(f,f.variables['ZSF'],name='ZSF')
        r4=FireNetcdf2Raster(f,f.variables['F_LINEINT2'],name='F_LINEINT2')
    
        open('UF.png','w').write(r1.getRaster(1))
        open('FGRNHFX.png','w').write(r2.getRaster(5))
        open('ZSF.png','w').write(r3.getRaster(0))
        open('F_LINEINT2.png','w').write(r4.getRaster(5))
    
        open('UF_c.png','w').write(r1.getColorbar(1))
        open('FGRNHFX_c.png','w').write(r2.getColorbar(5))
        open('ZSF_c.png','w').write(r3.getColorbar(0))
        open('F_LINEINT2_c.png','w').write(r4.getColorbar(5))

    f=FireRasterFile(wrfout)
    uf=f.rasterClassFromVar('UF')
    fgrnhfx=f.rasterClassFromVar('FGRNHFX')
    zsf=f.rasterClassFromVar('ZSF')
    flineint=f.rasterClassFromVar('F_LINEINT2')
    nfuelcat=f.rasterClassFromVar('NFUEL_CAT')
    lfn=f.rasterClassFromVar('LFN')

    n=ncKML()
    n.setViewFromRaster(nfuelcat)
    n.groundOverlayFromRaster(uf)
    n.groundOverlayFromRaster(fgrnhfx)
    n.groundOverlayFromRaster(zsf)
    n.groundOverlayFromRaster(flineint)
    n.groundOverlayFromRaster(nfuelcat)
    n.groundOverlayFromRaster(lfn)
    n.polygonFromContour(lfn,contour=0)

    n.savekmz('fire.kmz')

def main(wrfout,vars):
    '''
    Main function called when the module is used a commandline script.
    '''

    f=FireRasterFile(wrfout)
    n=ncKML()
    r=f.firePerimeterClass()
    try:
        # to make this work for non-fire output files as well
        # don't bomb out if LFN doesn't exist
        n.setViewFromRaster(r)
        n.polygonFromContour(r,contour=0)
    except Exception:
        pass
    for v in vars:
        try:
            r=f.rasterClassFromVar(v)
            n.groundOverlayFromRaster(r)
        except Exception:
            print 'Error processing %s... skipping.' % v

    n.savekmz('wrf.kmz')

if __name__ == '__main__':
    import sys
    if len(sys.argv) <= 1:
        print 'usage: %s wrfout [var1 [var2 [ ... ] ] ]' % sys.argv[0]
        print 'Outputs wrf.kmz containing the ground overlays'
        print 'of the variables specified on the commandline.'
        sys.exit(1)

    main(sys.argv[1],sys.argv[2:])
