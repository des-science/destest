import numpy as np
import fitsio as fio
import h5py
import cPickle as pickle
import os
import sys
import yaml
import cProfile, pstats
import time
# and maybe a bit optimistic...
from multiprocessing import Pool
# from mpi4py import MPI
from sharedNumpyMemManager import SharedNumpyMemManager as snmm 

if sys.version_info[0] == 3:
    string_types = str,
else:
    string_types = basestring,


def save_obj( obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj( name ):
    with open(name, 'rb') as f:
        return pickle.load(f)


def convert_mcal_to_h5( catdir='/global/project/projectdirs/des/wl/desdata/users/esheldon/matched-catalogs/', file_ext='fits', outfile='./y3v02-mcal-001-combined-blind-v1.h5' ):
    """
    Converts metacal fits files into a single h5 file with separate tables for each of the unsheared and four sheared values.
    """
    import glob
    from numpy.lib.recfunctions import append_fields

    def append_sheared_table( f, cat, i, sheared, iter_end ):
        cols = [name for name in cat.dtype.names if (name[-3:] == '_'+sheared)]
        tmp  = cat[cols]
        tmp.dtype.names = [name[:-3] for name in tmp.dtype.names]
        tmp = append_fields(tmp,'size_ratio',data=tmp['size']/cat['psf_size'])
        if i==0:
            f.create_dataset('catalog/sheared_'+sheared, maxshape=(None), shape=(total_length,), dtype=tmp.dtype)
        f['catalog/sheared_'+sheared][iter_end:iter_end+len(cat)] = tmp

    total_length = 0
    for i,name in enumerate(glob.glob(catdir+'*')):
        if file_ext not in name:
            continue
        total_length += fio.FITS(name)[1].read_header()['NAXIS2']

    f = h5py.File(outfile, 'w')
    iter_end = 0
    for i,name in enumerate(glob.glob(catdir+'*')):
        if file_ext not in name:
            continue
        print name
        cat=fio.FITS(name)[1].read()
        cols = [name for name in cat.dtype.names if (not name[-3:] in ['1p','1m','2p','2m'])]
        tmp = cat[cols]
        tmp = append_fields(tmp,'size_ratio',data=tmp['size']/tmp['psf_size'])
        if i==0:
            f.create_dataset('catalog/unsheared', maxshape=(None), shape=(total_length,), dtype=tmp.dtype)
        f['catalog/unsheared'][iter_end:iter_end+len(tmp)] = tmp
        for sheared in ['1p','1m','2p','2m']:
            append_sheared_table(f,cat,i,sheared,iter_end)
        iter_end += len(cat)

    f.close()


def convert_mcal_to_h5_v2( catdir='/global/project/projectdirs/des/wl/desdata/users/esheldon/matched-catalogs/', file_ext='fits', outfile='./y3v02-mcal-001-combined-blind-v1.h5' ):
    """
    Converts metacal fits files into a single h5 file with separate tables for each of the unsheared and four sheared values.
    This form is much faster to access, but doesn't preserve the complete recarray table structure.
    """
    import glob
    from numpy.lib.recfunctions import append_fields

    def append_sheared_table( f, cat, i, sheared, iter_end, total_length ):
        cols = [name for name in cat.dtype.names if (name[-2:] == sheared)]
        for name in cols:
            if i==0:
                f.create_dataset('catalog/sheared_'+sheared+'/'+name[:-3], maxshape=(total_length,), shape=(total_length,), dtype=cat.dtype[name], chunks=(1000000,))
            f['catalog/sheared_'+sheared+'/'+name[:-3]][iter_end:iter_end+len(cat)] = cat[name]
        if i==0:
            f.create_dataset('catalog/sheared_'+sheared+'/size_ratio', maxshape=(total_length,), shape=(total_length,), dtype=cat.dtype['size'], chunks=(1000000,))
        f['catalog/sheared_'+sheared+'/size_ratio'][iter_end:iter_end+len(cat)] = cat['size_'+sheared]/cat['psf_size']

    total_length = 0
    for i,name in enumerate(glob.glob(catdir+'*')):
        if file_ext not in name:
            continue
        total_length += fio.FITS(name)[1].read_header()['NAXIS2']

    f = h5py.File(outfile, 'w')
    iter_end = 0
    for i,fname in enumerate(glob.glob(catdir+'*')):
        if file_ext not in fname:
            continue
        print fname
        cat=fio.FITS(fname)[1].read()
        cols = [name for name in cat.dtype.names if (not name[-2:] in ['1p','1m','2p','2m'])]
        if i!=0:
            for name in cols:
                if i==0:
                    f.create_dataset('catalog/unsheared/'+name, maxshape=(total_length,), shape=(total_length,), dtype=cat.dtype[name], chunks=(1000000,))
                f['catalog/unsheared/'+name][iter_end:iter_end+len(cat)] = cat[name]
            if i==0:
                f.create_dataset('catalog/unsheared/size_ratio', maxshape=(total_length,), shape=(total_length,), dtype=cat.dtype['size'], chunks=(1000000,))
            f['catalog/unsheared/size_ratio'][iter_end:iter_end+len(cat)] = cat['size']/cat['psf_size']
            for sheared in ['1p','1m','2p','2m']:
                append_sheared_table(f,cat,i,sheared,iter_end,total_length)
        iter_end += len(cat)

    f.close()

def file_path( params, subdir, name, var=None, var2=None, ftype='txt' ):
    """
    Set up a file path, and create the path if it doesn't exist.
    """

    if var is not None:
        name += '_' + var
    if var2 is not None:
        name += '_' + var2
    name += '.' + ftype

    fpath = os.path.join(params['output'],params['param_file'][:params['param_file'].index('.')],subdir)

    if os.path.exists(fpath):
        if not params['output_exists']:
            raise IOError('Output directory already exists. Set output_exists to True to use existing output directory at your own peril.')
    else:
        os.mkdir(fpath)
        params['output_exists'] = True

    return os.path.join(fpath,name)


def write_table( params, table, subdir, name, var=None, var2=None):
    """
    Save a text table to file. Table must be a numpy-compatible array.
    """

    fpath = file_path(params,subdir,name,var=var,var2=var2)
    np.savetxt(fpath,table)


def child_testsuite( calc ):

    params, selector, calibrator = calc

    Testsuite( params, selector=selector, calibrator=calibrator, child=True )


def scalar_sum(x,length):
    # catches scalar weights, responses and multiplies by the vector length for the mean
    if np.isscalar(x):
        return x*length
    return np.sum(x)


class Testsuite(object):
    """
    Testsuite manager class. Initiated with a yaml file path or dictionary that contains necessary settings and instructions.
    """

    def __init__( self, param_file, selector = None, calibrator = None, child = False, **kwargs ):

        # Read in yaml information
        if isinstance(param_file, string_types):

            self.params     = yaml.load(open(param_file))
            self.params['param_file'] = param_file
            # Archive yaml file used for results
            self.save_input_yaml()

        else:

            self.params     = param_file

        # Set up classes that manage data
        # Source is an HDF5 file.
        if self.params['source'] == 'hdf5':
            self.source  = H5Source(self.params)

        # Source is a desdm table
        elif self.params['source'] == 'desdm':
            self.source  = DESDMSource(self.params)

        # Source is an LSST thing
        elif self.params['source'] == 'lsstdm':
            self.source  = LSSTDMSource(self.params)

        else:
            raise NameError('Something went wrong with parsing your source.')

        if selector is None:
            self.selector    = Selector(self.params,self.source)
        else:
            self.selector    = selector

        if calibrator is None:
            if self.params['cal_type'] == None:
                self.calibrator = NoCalib(self.params,self.selector)
            elif self.params['cal_type'] == 'mcal':
                self.calibrator = MetaCalib(self.params,self.selector)
            elif self.params['cal_type'] == 'classic':
                self.calibrator = ClassCalib(self.params,self.selector)
            else:
                raise NameError('Something went wrong with parsing your calibration type.')
        else:
            self.calibrator = calibrator

        # Run tests
        if 'split_mean' in self.params:

            if self.params['use_mpi'] and (not child):
                procs = comm.Get_size()
                iter_list = [self.params['split_x'][i::procs] for i in xrange(procs)]
                calcs = []
                for proc in range(procs):
                    if iter_list[proc] == []:
                        continue
                    params = self.params.copy()
                    params['split_x'] = iter_list[proc]
                    calcs.append((params,self.selector,self.calibrator))
                pool.map(child_testsuite, calcs)
            else:
                LinearSplit(self.params,self.selector,self.calibrator,self.source,self.params['split_x'],self.params['split_mean'])


    def save_input_yaml( self ):
        """
        Arxiv input yaml settings.
        """

        fpath = file_path(self.params,'',self.params['param_file'][:self.params['param_file'].index('.')],ftype='yaml')
        print 'saving input yaml to: '+fpath
        with open(fpath, 'w') as outfile:
            yaml.dump(self.params, outfile, default_flow_style=False)


class SourceParser(object):
    """
    A class to manage the actual reading or downloading of data from external sources. 
    Initiate with a testsuite param dictionary.
    To use later: source_parser.read(...). All the messy details are hidden.
    """

    def __init__( self, params ):
        self.params = params
        self.open()

    def open( self ):
        raise NotImplementedError('Subclass '+self.__class__.__name__+' should have method open().')

    def read( self ):
        raise NotImplementedError('Subclass '+self.__class__.__name__+' should have method read().')

    def close( self ):
        raise NotImplementedError('Subclass '+self.__class__.__name__+' should have method close().')


class H5Source(SourceParser):
    """
    A class to manage the actual reading or downloading of data from HDF5 sources. 
    """

    def __init__( self, params ):

        super(H5Source,self).__init__(params)

        if 'filename' not in self.params.keys():
            raise NameError('Must provide a filename for hdf5 source.')
        if 'table' not in self.params.keys():
            raise NameError('Must specify table name for hdf5 file.')
        if type(self.params['table']) is not list:
            raise TypeError('Table must be provided as a list of names (even a list of one).')

        if 'group' in self.params.keys():

            self.hdf = h5py.File(self.params['filename'], mode = 'r')
            # save all column names
            self.cols = self.hdf[self.params['group']][self.params['table'][0]].keys()
            # save length of tables
            self.size = self.hdf[self.params['group']][self.params['table'][0]][self.cols[0]].shape[0]

            # Loop over tables and save convenience information
            for t in self.params['table']:
                keys = self.hdf[self.params['group']][t].keys()
                if self.hdf[self.params['group']][t][keys[0]].shape[0] != self.size:
                    raise TypeError('Length of sheared tables in hdf5 file must match length of unsheared table.')

            if len(self.params['table'])>1:
                # save sheared column names
                self.sheared_cols = self.hdf[self.params['group']][self.params['table'][1]].keys()

        else:

            raise NameError('Need group name for hdf5 file.')

        self.close()

    def open( self ):

            self.hdf = h5py.File(self.params['filename'], mode = 'r+')

    def read( self, col=None, rows=None, nosheared=False ):

        self.open()

        def add_out( table, rows, col ):
            """
            Extract a portion of a column from the file.
            """

            if rows is not None:
                if hasattr(rows,'__len__'):
                    if len(rows==2):
                        out = self.hdf[self.params['group']][table][col][rows[0]:rows[1]] 
                else:
                    out = self.hdf[self.params['group']][table][col][rows] 
            else:
                out = self.hdf[self.params['group']][table][col][:] 

            return out

        if col is None:
            raise NameError('Must specify column.')

        out = []
        # For metacal file, loop over tables and return list of 5 unsheared+sheared values for column (or just unsheraed if 'nosheared' is true or there doesn't exist sheared values for this column) 
        # For classic file, get single column.
        # For both metacal and classic files, output is a list of columns (possible of length 1)
        for i,t in enumerate(self.params['table']):
            if i==0:
                if col not in self.hdf[self.params['group']][t].keys():
                    raise NameError('Col '+col+' not found in hdf5 file.')
            else:
                if nosheared:
                    continue
                if col not in self.sheared_cols:
                    continue
                if col not in self.hdf[self.params['group']][t].keys():
                    raise NameError('Col '+col+' not found in sheared table '+t+' of hdf5 file.')

            if rows is not None:
                out.append( add_out(t,rows,col) )
            else:
                out.append( add_out(t,None,col) )

        self.close()

        return out

    def close( self ):

        if hasattr(self,'hdf'):
            self.hdf.close()

class DESDMSource(SourceParser):
    """
    A class to manage the actual reading or downloading of data from DESDM sources. 
    """

    def __init__( self ):
        raise NotImplementedError('You should write this.')


class LSSTDMSource(SourceParser):
    """
    A class to manage the actual reading or downloading of data from LSSTDM sources. 
    """

    def __init__( self ):
        raise NotImplementedError('You should write this.')


class Selector(object):
    """
    A class to manage masking and selections of the data.
    Initiate with a testsuite object.
    Initiation will parse the 'select_cols' conditions in the yaml file and create a limiting mask 'mask_', ie, an 'or' of the individual unsheared and sheared metacal masks. The individual masks (list of 1 or 5 masks) are 'mask'.
    """

    def __init__( self, params, source ):
        self.params    = params
        self.source    = source
        self.build_limiting_mask()

    def build_limiting_mask( self ):
        """
        Build the limiting mask for use in discarding any data that will never be used.
        """

        mask = None

        # Setup mask file cache path.
        mask_file = file_path(self.params,'cache','mask',ftype='pickle')
        if self.params['load_cache']:
            # if mask cache exists, read mask from pickle and skip parsing yaml selection conditions.

            if os.path.exists(mask_file):
                mask, mask_ = load_obj(mask_file)

        if mask is None:
            # mask cache doesn't exist, or you chose to ignore it, so masks are built from yaml selection conditions
            # set up 'empty' mask
            mask = [np.ones(self.source.size, dtype=bool)]
            if self.params['cal_type']=='mcal':
                mask = mask * 5

            # For each of 'select_cols' in yaml file, read in the data and iteratively apply the appropriate mask
            for i,select_col in enumerate(self.params['select_cols']):
                cols = self.source.read(col=select_col)
                for j,col in enumerate(cols):
                    mask[j] = mask[j] & eval(self.params['select_exp'][i])

            # Loop over unsheared and sheared mask arrays and build limiting mask
            mask_ = np.zeros(self.source.size, dtype=bool)
            for imask in mask:
                mask_ = mask_ | imask

            print len(mask_),np.sum(mask_)

            # Cut down masks to the limiting mask
            # Its important to note that all operations will assume that data has been trimmed to satisfy selector.mask_ from now on
            for i in range(len(mask)):
                mask[i] = mask[i][mask_]
            mask_ = np.where(mask_)[0]
            print len(mask_),len(mask[0]),np.sum(mask[0])

            # save cache of masks to speed up reruns
            save_obj( [mask, mask_], mask_file )

        self.mask_ = snmm.createArray((len(mask_),), dtype=np.int64)
        snmm.getArray(self.mask_)[:] = mask_[:]
        mask_ = None

        self.mask = []
        for i in range(len(mask)):
            self.mask.append( snmm.createArray((len(mask[i]),), dtype=np.bool) )
            snmm.getArray(self.mask[i])[:] = mask[i][:]
            mask[i] = None

    def get_col( self, col, nosheared=False, uncut=False ):
        """
        Wrapper to retrieve a column of data from the source and trim to the limiting mask (mask_)
        """

        # x at this point is the full column
        x = self.source.read(col=col, nosheared=nosheared)

        # trim and return
        for i in range(len(x)):
            x[i] = x[i][snmm.getArray(self.mask_)]
            print 'uncut',len(x[i])
        if uncut:
            return x

        for i in range(len(x)):
            x[i] = x[i][snmm.getArray(self.mask[i])]
            print 'cut',len(x[i])
        return x

    def get_masked( self, x, mask ):
        """
        Accept a mask and column(s), apply the mask jointly with selector.mask (mask from yaml selection) and return masked array.
        """

        if mask is None:
            mask = [np.s_[:]]*5

        if type(mask) is not list:
            mask = [ mask ]

        if type(x) is not list:
            if np.isscalar(x):
                return x
            else:
                return x[snmm.getArray(self.mask[0])][mask[0]]

        if np.isscalar(x[0]):
            return x

        return [ x_[snmm.getArray(self.mask[i])][mask[i]] for i,x_ in enumerate(x) ]

    def get_mask( self, mask=None ):
        """
        Same as get_masked, but only return the mask.
        """

        if mask is None:
            return [ np.where(snmm.getArray(self.mask[i]))[0] for i in range(len(self.mask)) ]

        return [ np.where(snmm.getArray(self.mask[i]))[0][mask_] for i,mask_ in enumerate(mask) ]


class Calibrator(object):
    """
    A class to manage calculating and returning calibration factors and weights for the catalog.
    Initiate with a testsuite params object.
    When initiated, will read in the shear response (or m), additive corrections (or c), and weights as requested in the yaml file. These are a necessary overhead that will be stored in memory, but truncated to the limiting mask (selector.mask_), so not that bad.
    """

    def __init__( self, params, selector ):

        self.params = params
        self.selector = selector

    def get_w(self,mask):
        """
        Get the weights and the sum of the weights.
        """

        w  = self.selector.get_masked(self.w,mask)
        mask_ = self.selector.get_mask(mask)
        ws = [ scalar_sum(w_,len(mask_[i])) for i,w_ in enumerate(w)]
        return w,ws

    def calibrate(self,col,mask=None,return_full_w=False,weight_only=False):
        """
        Return the calibration factor and weights, given potentially an ellipticity and selection.
        """

        # Get the weights
        w,ws = self.get_w(mask)
        if return_full_w:
            w_ = w
        else:
            w_ = w[0]
        if weight_only:
            return w_

        # Get a selection response
        Rs = self.select_resp(col,mask,w,ws)

        # Check if an ellipticity - if so, return real calibration factors
        if col == self.params['e'][0]:
            Rg1 = self.selector.get_masked(snmm.getArray(self.Rg1),mask)[0]
            R = np.sum(Rg1*w[0])/ws[0]
            print 'R',R,np.sum(Rg1*w[0]),ws[0],np.mean(Rg1)
            R += Rs
            c = self.selector.get_masked(self.c1,mask)
            return R,c,w_
        elif col == self.params['e'][1]:
            Rg2 = self.selector.get_masked(snmm.getArray(self.Rg2),mask)[0]
            R = np.sum(Rg2*w[0])/ws[0]
            R += Rs
            c = self.selector.get_masked(self.c2,mask)
            return R,c,w_
        else:
            return None,None,w_

    def select_resp(self,col,mask,w,ws):
        """
        Return a zero selection response (default).
        """
        return 0.


class NoCalib(Calibrator):
    """
    A class to manage calculating and returning calibration factors and weights for a general catalog without shear corrections.
    """

    def __init__( self, params, selector ):

        super(NoCalib,self).__init__(params,selector)

        self.Rg1 = self.Rg2 = 1.
        self.c1 = self.c2 = 0.
        self.w = [1]
        if 'w' in self.params:
            self.w = self.selector.get_col(self.params['w'])


class MetaCalib(Calibrator):
    """
    A class to manage calculating and returning calibration factors and weights for a metacal catalog.
    """

    def __init__( self, params, selector ):

        super(MetaCalib,self).__init__(params,selector)

        self.Rg1 = self.Rg2 = 1.
        if 'Rg' in self.params:
            Rg1 = self.selector.get_col(self.params['Rg'][0],uncut=True)[0]
            Rg2 = self.selector.get_col(self.params['Rg'][1],uncut=True)[0]
            e1  = self.selector.get_col(self.params['e'][0],nosheared=True,uncut=True)[0]
            e2  = self.selector.get_col(self.params['e'][1],nosheared=True,uncut=True)[0]
            self.Rg1 = snmm.createArray((len(Rg1),), dtype=np.float64)
            snmm.getArray(self.Rg1)[:] = Rg1[:]
            Rg1 = None
            self.Rg2 = snmm.createArray((len(Rg2),), dtype=np.float64)
            snmm.getArray(self.Rg2)[:] = Rg2[:]
            Rg2 = None
            self.e1 = snmm.createArray((len(e1),), dtype=np.float64)
            snmm.getArray(self.e1)[:] = e1[:]
            e1 = None
            self.e2 = snmm.createArray((len(e2),), dtype=np.float64)
            snmm.getArray(self.e2)[:] = e2[:]
            e2 = None
        self.c1 = self.c2 = 0.
        if 'c' in self.params:
            self.c1 = self.selector.get_col(self.params['c'][0],uncut=True)
            self.c2 = self.selector.get_col(self.params['c'][1],uncut=True)
        self.w = [1] * 5
        if 'w' in self.params:
            self.w = self.selector.get_col(self.params['w'])

    def select_resp(self,col,mask,w,ws):
        """
        Get the selection response.
        """

        # if an ellipticity column, calculate and return the selection response and weight
        if col in self.params['e']:
            if len(mask)==1: # exit for non-sheared column selections
                return 0.
            mask_ = [ snmm.getArray(imask) for imask in self.selector.get_mask(mask) ]

        if col == self.params['e'][0]:
            Rs = np.sum(snmm.getArray(self.e1)[mask_[1]]*w[1])/ws[1] - np.sum(snmm.getArray(self.e1)[mask_[2]]*w[2])/ws[2]
            print 'Rs',Rs,np.sum(snmm.getArray(self.e1)[mask_[1]]*w[1])/ws[1],np.sum(snmm.getArray(self.e1)[mask_[2]]*w[2])/ws[2],np.sum(snmm.getArray(self.e1)[mask_[1]]*w[1]),ws[1],np.mean(snmm.getArray(self.e1)[mask_[1]])
        elif col == self.params['e'][1]:
            Rs = np.sum(snmm.getArray(self.e2)[mask_[3]]*w[3])/ws[3] - np.sum(snmm.getArray(self.e2)[mask_[4]]*w[4])/ws[4]
        else:
            return 0.

        Rs /= 2.*self.params['dg']

        return Rs


class ClassicCalib(Calibrator):
    """
    A class to manage calculating and returning calibration factors and weights for a metacal catalog.
    """

    def __init__( self, params, selector ):

        super(ClassCalib,self).__init__(params,selector)

        self.Rg1 = self.Rg2 = 1.
        if 'Rg' in self.params:
            self.Rg1 = self.selector.get_col(self.params['Rg'][0])
            self.Rg2 = self.selector.get_col(self.params['Rg'][1])

        self.c1 = self.c2 = 0.
        if 'c' in self.params:
            self.c1 = self.selector.get_col(self.params['c'][0])
            self.c2 = self.selector.get_col(self.params['c'][1])

        self.w  = [1]
        if 'w' in self.params:
            self.w = self.selector.get_col(self.params['w'])


class Splitter(object):
    """
    A class for managing splitting the data set into bins and accessing the binned data.
    Initiate with a testsuite object.
    """

    def __init__( self, params, selector, calibrator, source, nbins = None ):

        self.params     = params
        self.selector   = selector
        self.calibrator = calibrator
        self.source     = source
        self.bins       = self.params['linear_bins']
        self.x          = None
        self.y          = None
        self.xcol       = None
        self.ycol       = None
        self.order      = None

        if 'split_x' in self.params:
            for col in self.params['split_x']:
                if col not in self.source.cols:
                    raise NameError(col + ' not in source.')
        else:
            self.params['split_x'] = self.source.cols

        if nbins is not None:
            self.bins = nbins

        return

    def get_x( self, col, xbin=None, return_mask=False ):
        """
        Get the 'x' column - the column you're binning the data by. 
        If you haven't already called splitter with this x col, the data will be read from the source and the binning edges will be set up.
        Optionally give a bin number, it will return the portion of the x array that falls in that bin. Can also optionally return the mask for that bin.
        """

        # If column doesn't already exist in splitter, read the data and define bins in self.split().
        if col != self.xcol:
            self.xcol  = col
            self.order = None
            self.split(col)

        # If not asking for a bin selection, return
        if xbin is None:
            return

        # If asking for a bin selection, find the appropriate mask and return that part of the x array.
        start,end = self.get_bin_edges(xbin)
        # print 'returning x bin',start,end
        mask      = [np.s_[start_:end_] for start_,end_ in tuple(zip(start,end))] # np.s_ creates an array slice 'object' that can be passed to functions
        mask      = [ order_[mask_] for order_,mask_ in tuple(zip(self.order,mask)) ]
        if return_mask:
            return self.x[start[0]:end[0]],mask
        return self.x[start[0]:end[0]]

    def get_y( self, col, xbin=None, return_mask=False ):
        """
        Get the 'y' column - the column you're doing stuff with in bins of the x col. If you haven't called splitter.get_x(), an error will be raised, since you haven't defined what you're binning against.  
        If you haven't already called splitter with this y col, the data will be read from the source.
        Optionally give a bin number, it will return the portion of the y array that falls in that bin. Can also optionally return the mask for that bin.
        """

        if self.xcol is None:
            raise NameError('There is no x column associated with this splitter.')

        # If column doesn't already exist in splitter, read the data and order it to match x ordering for efficient splitting.
        if col != self.ycol:
            self.ycol = col
            self.y = self.selector.get_col(col,nosheared=True)
            for i,y_ in enumerate(self.y):
                self.y[i] = y_[self.order[i]]
            self.y = self.y[0]

        # If not asking for a bin selection, return
        if xbin is None:
            return

        # If asking for a bin selection, find the appropriate mask and return that part of the y array.
        start,end = self.get_bin_edges(xbin)
        # print 'returning y bin',start,end
        mask      = [np.s_[start_:end_] for start_,end_ in tuple(zip(start,end))]
        mask      = [ order_[mask_] for order_,mask_ in tuple(zip(self.order,mask)) ]
        if return_mask:
            return self.y[start[0]:end[0]],mask
        return self.y[start[0]:end[0]]

    def split( self, col ):
        """
        Reads in a column (x) and sorts it. If you allowed cache reading, it will check if you've already done this and just read it in from the pickle cach. Then finds the edges of the bins you've requested.
        """

        # Check if cache file exists and use it if you've requested that.
        sort_file = file_path(self.params,'cache','sort',var=col,ftype='pickle')
        if self.params['load_cache']:
            print 'loading split sort cache',sort_file

            if os.path.exists(sort_file):
                self.order,self.x = load_obj(sort_file)

        # Cache file doesn't exist or you're remaking it
        if self.order is None:
            print 'split sort cache not found'
            # Read x
            self.x     = self.selector.get_col(col)
            # save the index order to sort the x array for more efficient binning
            self.order = []
            for i,x_ in enumerate(self.x):
                self.order.append( np.argsort(x_) )
                self.x[i] = x_[self.order[i]]
            # save cache of sorted x and its order relative to the source
            save_obj( [self.order,self.x], sort_file )

        # get bin edges
        self.get_edge_idx()
        self.x = self.x[0]

        return

    def get_edge_idx( self ):
        """
        Find the bin edges that split the data into the ranges you set in the yaml or into a number of equal-weighted bins.
        """

        self.edges = []
        # You've provided a number of bins. Get the weights and define bin edges such that there exists equal weight in each bin.
        if not self.params['split_by_w']:

            for x_ in self.x:
                xw = np.ones(len(x_))
                normcumsum = xw.cumsum() / xw.sum()
                self.edges.append( np.searchsorted(normcumsum, np.linspace(0, 1, self.bins+1, endpoint=True)) )

        else:

            w,R = self.calibrator.calibrate(self.xcol,return_full_w=True,weight_only=True,include_Rg=True)
            for x_,w_ in tuple(zip(self.x,w)):
                xw = w_*R
                normcumsum = xw.cumsum() / xw.sum()
                self.edges.append( np.searchsorted(normcumsum, np.linspace(0, 1, self.bins+1, endpoint=True)) )

        return

    def get_bin_edges( self, xbin ):
        """
        Helper function to return the lower and upper bin edges.
        """
        return [edge[xbin] for edge in self.edges],[edge[xbin+1] for edge in self.edges]


class LinearSplit(object):
    """
    Test class to do linear splitting (operations on binned data not at the 2pt level).
    Instantiate with a testsuite object and opetionally a function to operate on the bins (not fully implemented).
    """

    def __init__( self, params, selector, calibrator, source, split_x, split_y, nbins = None, func=np.mean, **kwargs ):

        self.params = params
        self.source = source
        if self.params['split_mean'] is not None:
            for col in self.params['split_mean']:
                if col not in self.source.cols:
                    raise NameError(col + ' not in source.')
        else:
            self.params['split_mean'] = self.source.cols

        self.calibrator = calibrator
        self.splitter   = Splitter(params,selector,calibrator,source,nbins)
        self.split_x    = split_x
        self.split_y    = split_y
        self.step       = 0

        # 'step' and this separate call is meant as a placeholder for potential parallelisation
        self.iter_mean()

    def iter_mean( self ):
        """
        Loop over x columns (quantities binned by) and y columns (quantities to perform operations on in bins of x), perform the operations, and save the results
        """

        for x in self.split_x:
            print 'x col',x
            for y in self.split_y:
                n     = []
                xmean = []
                xlow  = []
                xhigh = []
                ymean = []
                ystd  = []
                print 'y col',y
                for xbin in range(self.splitter.bins):
                    # get x array in bin xbin
                    xval       = self.splitter.get_x(x,xbin)
                    n.append( len(xval) )
                    xlow.append( xval[0] )
                    xhigh.append( xval[-1] )
                    # get mean values of x in this bin
                    xmean.append( self.mean(x,xval,return_std=False) )
                    # get y array in bin xbin
                    yval,mask  = self.splitter.get_y(y,xbin,return_mask=True)
                    # get mean and std (for error) in this bin
                    ymean_,ystd_ = self.mean(y,yval,mask=mask)
                    ymean.append( ymean_ )
                    ystd.append( ystd_/np.sqrt(n[xbin]) )

                # Save results
                table = np.array([n,xlow,xmean,xhigh,ymean,ystd]).T
                print 'mean',table
                write_table(self.params, table,'test_output','linear_split',var=x,var2=y)

    def mean( self, col, x, mask=None, return_std=True, return_rms=False ):
        """
        Function to do mean, std, rms calculations
        """

        # Get response and weight.
        if mask is None:
            R,c,w = self.calibrator.calibrate(col)
        else:
            R,c,w = self.calibrator.calibrate(col,mask=mask)
        print 'Rcw',col,R,c,w

        # do the calculation
        if R is not None:

            x  = np.copy(x)-c
            Rw = scalar_sum(w*R,len(x))
            if return_std:
                Rw2 = scalar_sum(w*R**2,len(x))

        else:

            Rw  = scalar_sum(w,len(x))
            if return_std:
                Rw2 = Rw

        mean = np.sum(w*x)/Rw
        if not (return_std or return_rms):
            return mean
        if return_std:
            std=np.sqrt(np.sum(w*(x-mean)**2)/Rw2)
            if not return_rms:
                return mean,std
        if return_rms:
            rms=np.sqrt(np.sum((w*x)**2)/Rw)
            if not return_std:
                return mean,rms

        return mean,std,rms

# pr = cProfile.Profile()

if __name__ == "__main__":
    """
    """
    # pr.enable()

    # from mpi_pool import MPIPool
    # comm = MPI.COMM_WORLD
    # pool = MPIPool(comm)
    # if not pool.is_master():
    #     pool.wait()
    #     sys.exit(0)

    Testsuite( sys.argv[1] )

    # pool.close()


    # pr.disable()
    # ps = pstats.Stats(pr).sort_stats('time')
    # ps.print_stats(20)

