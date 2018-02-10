# destest

To utilize catalog interface:

params      = yaml.load(open(param_file))
source      = destest.H5Source(params)
selector    = destest.Selector(params,source)

Get a column:

selector.get_col(col) # returns a column from the catalog for column name col. If you've instantiated selector with an appropriate yaml file input, it will have used the correct cuts that you specified in the yaml file. You can use nosheared = True to not get sheared columns (it will return a tuple).

source.read(col=col) # lower level function that won't manage cuts, can also pass row for a row selection or unsheared=True to not get sheared columns.

Calibrator class to get response and weights:

if params['cal_type'] == None:
    calibrator = NoCalib(params,selector)
elif params['cal_type'] == 'mcal':
    calibrator = MetaCalib(params,selector)
elif params['cal_type'] == 'classic':
    calibrator = ClassCalib(params,selector)
    
Then:

R,c,w = calibrator.calibrate(col,mask=mask) # col is the column name string...
    
There is also a Splitter class to manage fast slicing for binning. See class LinearSplit for use.

Ex.
splitter.get_x(x) # This will nicely return a sorted column for column name x for fast bin slicing, and with get_y() return appropriately sorted column that matches the sort of x.
splitter.get_x(x,xbin=n) # return the nth bin in x based on the number of bins specified in the yaml file (similar for matching get_y().)
