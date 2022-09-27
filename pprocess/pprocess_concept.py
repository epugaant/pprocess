#!/usr/bin/env python
# IDENT         pprocess_concept.py
# LANGUAGE      Python 3.X
# AUTHOR        E. Puga
# PURPOSE       Prototype script to test concept conversion to VOTABLE of two types of catalogues:
#               - fits table data (metadata rich)
#               - generic data table (structured array) + table_info
#               It makes use of astropy Table, FITS and VOTABLE packages
# 
# EXAMPLE OF USAGE pprocess_concept.py 


import os
import glob
from typing import OrderedDict
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table, vstack

from astropy.io.fits.connect import read_table_fits
from astropy.io.votable.table import from_table
from astropy import units as u

def create_directory(path):
    '''Makes a folder (and its parents) if not present'''
    try:
        os.makedirs(path)
    except:
        pass
    return


def make_cols_info_table(tbl):
    '''Construct an info table from the columns metadata hidden in the FITS table Comment Keyword mixin columns.
    
    This is not part of the table.meta, but in the table.Columns.info attribute 
    
    Parameters
    ----------
    tbl : `~astropy.table.Table`
        input table

    Returns
    -------
    `~astropy.table.Table`
        table info with numpy fields:  name | dtype | unit | description | class | n_bad
    '''
    return tbl.info(out=None)

def table_to_votable(mr_tbl):
    '''   Function to cast a `~astropy.table.Table` into a `~astropy.io.votable.tree.Table`
    The oposite is a built in function `~astropy.io.votable.tree.Table.to_table`
    
    Parameters
    ----------
    mr_tbl : `~astropy.table.Table`
        metadata rich astropy table 

    Returns
    -------
    `~astropy.io.votable.tree.Table`
        votable tree table
    '''
    
    #TODO: serialize tables in resource construction
    from astropy.io.votable.tree import VOTableFile, Resource, Table
    #from tree import VOTableFile, Resource, Table

    votable = VOTableFile()
    # ...with one resource...
    resource = Resource()
    votable.resources.append(resource)
    # ... with one table
    resource.tables.append(Table.from_table(votable, mr_tbl))

    return votable

def split_bytes_dtype_list(dtype_info_list):
    '''Function to convert the info list dtypes bytes<N> (or str<N>) to numpy types.
    The info method for `~astropy.table.Table` contains the mixin columns info, but the datatypes are human readable types.
    This function is the inverse of dtype_info_name in /Users/epuga/opt/anaconda3/envs/pprocessor/lib/python3.10/site-packages/astropy/utils/data_info.py. 
    (dtype_info_name('S1') --> bytes1)

    Parameters
    ----------
    dtype_info_list : list
        The output is shown below for ``bytes`` and ``str`` types, with <N> being
        the number of characters. This representation corresponds to the Python
        type that matches the dtype::

            Numpy   Python      
            S<N>    bytes<N>   
            U<N>    str<N>
        (dtype_info_name('S1') --> 'bytes1')
        ('S1'<--split_bytes_dtype_list('bytes1') )
        
    Returns
    -------
    list
        modified dtype_info_list with bytes or strings modified to ``numpy`` dtype.char
    '''
    import re
    dtype_info_list_out = []
    for dtype_info in dtype_info_list:
        try:
            if np.dtype(dtype_info).isbuiltin: #This throws a Typeerror with 'bytes3', needs the function split_bytes_dtype_list
                dtype_info_list_out.append(dtype_info)
        except:
            if dtype_info.startswith(('bytes', 'str')):
                (scalar_type, length) = tuple(filter(None, re.split(r'(\d+)', dtype_info)))
                dtype_info_list_out.append(np.dtype(scalar_type).char + length)
            
    return dtype_info_list_out

def table_is_compatible(ref_table, table, verbose=True):
    '''Function to check if one table is similar to a reference table (column names and dtype)

    Parameters
    ----------
    ref_table : `~astropy.table.Table`
        Reference table
    table : `~astropy.table.Table`
        generic table
    verbose : bool, optional
        print messages, by default True

    Returns
    -------
    bool
        True if the table is similar to a reference table
    '''

    flag = True
    if all(key in table.colnames for key in ref_table.keys()):
        if verbose:
            print('OK, Minimum required table keys are present')
    else:
        if verbose:
            print('Not all minimum table keys are present')
        for key in ref_table.keys():
            if key in table.colnames:
                if not np.issubsctype(table[key].dtype, key):
                    flag = False
                    if verbose:
                        print('Dtype of column {} does not match. ')
                        print('Got {}'.format(key, table[key].dtype))
            else:
                flag = False
                if verbose:
                    print('Name of column {} does not match.'.format(key))
    return flag

def append_tables(fits_table_list):
    '''Function to aggregate the tables in a list of .fits filenames taking the first 
    one as a reference

    Parameters
    ----------
    fits_table_list : list or str
        list of filenames (dirname and filename) of FITS files

    Returns
    -------
    `~astropy.table.Table`
        aggregated table of FITS files

    Raises
    ------
    ValueError
        Warning if there is only one FITS file in the list
    ValueError
        If any of the input table is not compatible with the first (reference table)
    '''
    if isinstance(fits_table_list, str):
        fits_table_list = list(fits_table_list)

    ref_table = Table.read(fits_table_list[0], character_as_bytes=False)
    tinfo = make_cols_info_table(ref_table)
    
    try:
        for fits_table_file in fits_table_list[1:]:
            add_table = Table.read(fits_table_file, character_as_bytes=False)
            if table_is_compatible(ref_table, add_table):
                ref_table = vstack([ref_table, add_table])
            else:
                print('{} is not compatible with the reference table.'.format(add_table))
                raise ValueError
    except:
        print('There is only one element to the list')
        raise ValueError
    

    return ref_table

def metadata_rich_table(descriptions, table=None, names=None, dtype=None, units=None, meta=None):
    '''Function to create a metadata-rich `~astropy.table.Table` or modify existing table metadata
    .. Note::
        Since the origin is likely a FITS, it creates the data as list-of-dicts, 
        therefore the generation of the new table makes use of rows. Attributes 
        are colnames, dtype, meta and info.

    Parameters
    ----------
    descriptions : list, dict
        List or dict of descriptions to apply to columns.
    table : list, optional
        Data (list-of-dicts used in BinTableHDU) to initialize the table, by default None if empty of data.
    names : list, optional
        Specify column names, by default None.
    dtype : list, optional
        Column data types, it is usually a list of double tuples, by default None.
    units : list, optional
        List or dict of units to apply to columns, generally needs modification with make_cols_info_table(table)['unit'].tolist(), by default None.
    meta : dict, optional
        Metadata associated with the table, by default None.

    Returns
    -------
    :class: `astropy.table.Table`
        A metadata-rich astropy table
    '''
    if names is None:
        names = table.colnames
        
    if dtype is None:
        dtype = table.dtype
    
    if units is None:
        units = make_cols_info_table(table)['unit'].tolist()
    
    if meta is None:
        meta = table.meta
    
    if table is None:
        new_table = Table(rows=None, names=names, dtype=dtype, meta=meta, descriptions=descriptions, units=units)
    else:
        new_table = Table(rows=table.as_array(), names=names, dtype=dtype, meta=meta, descriptions=descriptions, units=units)
        
    return new_table

def extend_votable_field(votable, attr_values, attr_name='ucd', attr_key=None):
    '''Function to extend attributes in all VOTable columns that maybe parsed but not explictely defined for 
    <FIELD> element of `~astropy.io.votable.tree.Table` (e.g. descriptions) and/or VOTable 
    specific metadata (ucd).

    Parameters
    ----------
    votable : `~astropy.io.votable.tree.VOTableFile`
        VOTable element that represents an entire file.
    attr_values : list, OrderedDict
        Values of attribute or sub-element in a <FIELD>.
    attr_name : str, optional
        Attribute or sub-element name in a <FIELD>, by default 'ucd'
    attr_key : str, optional
        Substring in key within the attr_values OrderedDict, if applicable, by default None.

    Returns
    -------
    `~astropy.io.votable.tree.VOTableFile`
        Modified `~astropy.io.votable.tree.VOTableFile` with extended attribute in all fields
    '''
    if attr_key is None:
        attr_key = attr_name.upper()
    if isinstance(attr_values, (OrderedDict, dict)):
        #filter the dictionary items whose key 
        #attr_values = OrderedDict({key: val for key, val in attr_values.items() if attr_key in key})
        #flag if attr_key substring is the core key in attr_values dict
        key_is_key = True if any(attr_key in k for k in attr_values.keys()) else False
    for i, column in enumerate(votable.iter_fields_and_params(), 1):
        try:
            if isinstance(attr_values, list):
                setattr(column, attr_name, attr_values[i])
            elif isinstance(attr_values, (OrderedDict, dict)):
                print(column)
                #if key_is_key, use the input substring as dict key, else the name 
                key = 'T'+attr_key+'{0}'.format(i) if key_is_key else column.name
                setattr(column, attr_name, attr_values[key])
                print(getattr(column, attr_name))
        except AttributeError as error:
            print('Exception found!: ' + str(error))
            print('Column {} has a problem with attribute {}'.format(column.name, attr_name)) #maybe use column.name?
            
    return votable
    

if __name__ == '__main__':
   
    test = False
    if test:
        path = fits.util.get_testdata_filepath('tb.fits')
    else:
        #path = '/Users/epuga/ESDC/FP7H2020/inputs/help/data/ELAIS-N1_20171016.fits' #Catalogue as FITS Table BinTable
        path = '/Volumes/EP_DISK2/EUProjects/HELP/spire_blind/dmu22_XID+SPIRE_COSMOS_BLIND_Matched_MF.fits'
        out_dirname = '/Users/epuga/ESDC/FP7H2020/votables/help'
        create_directory(out_dirname)
            
    filename = os.path.basename(path)
    basename = os.path.splitext(filename)[0]
    in_dirname = os.path.dirname(path)
    if test:
        out_dirname = in_dirname
    
    print('Input .FITS: {}'.format(os.path.join(in_dirname, filename)))
    print('Output VOTABLE: {}'.format(os.path.join(out_dirname, basename+'.xml')))
    
    #-------------

    #Test 1: Metadata-rich fits table to votable
    #Deconstruct into numpy_array and a separate table with columns for testing purposes
    table = read_table_fits(path, memmap=True, astropy_native=True)# more tunnable way to read in table 
    tinfo = make_cols_info_table(table)
    bare_data = table.as_array() # produces a structured array (records array)
    
    #Construct one metadata-rich table
    print(split_bytes_dtype_list(tinfo['dtype'].tolist()))
    mr_table = Table(rows=table.as_array(), names=tinfo['name'].tolist(), 
                     dtype=split_bytes_dtype_list(tinfo['dtype'].tolist()), meta=table.meta, 
                     descriptions=tinfo['description'].tolist(), units=tinfo['unit'].tolist())
    

    #Test 2: serialized fits tables + external descriptions to votable
    fits_list = glob.glob(os.path.join('/Volumes/EP_DISK2/EUProjects/HELP/spire_blind', "*.fits"))
    out_dirname = '/Users/epuga/ESDC/FP7H2020/votables/help'
    
    _all_fields = ['AKARI-NEP','AKARI-SEP','Bootes','CDFS-SWIRE','COSMOS','EGS','ELAIS-N1','ELAIS-N2','ELAIS-S1','GAMA-09','GAMA-12','GAMA-15','HDF-N','Herschel-Stripe-82','Lockman-SWIRE','NGP','SA13','SGP','SPIRE-NEP','SSDF','xFLS','XMM-13hr','XMM-LSS'] 
    fits_list_clean = [x for x in fits_list if "HELP_BLIND" not in x]
    fits_list_sort = sorted(fits_list_clean, key=lambda x: [i for i,e in enumerate(_all_fields) if e in x][0])
    
    #serialization fits tables into one table
    simple_table = append_tables(fits_list_sort)
    simple_table['Bkg_SPIRE_250'].unit = u.mJy/u.beam
    simple_table['Bkg_SPIRE_350'].unit = u.mJy/u.beam
    simple_table['Bkg_SPIRE_500'].unit = u.mJy/u.beam
    simple_table['Sig_conf_SPIRE_250'].unit = u.mJy/u.beam
    simple_table['Sig_conf_SPIRE_350'].unit = u.mJy/u.beam
    simple_table['Sig_conf_SPIRE_500'].unit = u.mJy/u.beam
    
    print('')
    print('Table Info before (units, dtype)')
    simple_table.info()
    
    td = ascii.read('/Users/epuga/ESDC/FP7H2020/HELP/spire_blind/column_unit_descriptors.txt', guess=False, 
                    names=['name', 'unit', 'descriptor'], format='no_header', delimiter='\t', fast_reader=False)
    
    mr_empty_table = metadata_rich_table(td['descriptor'].tolist(), names=simple_table.colnames, 
                                         dtype=simple_table.dtype, units=make_cols_info_table(simple_table)['unit'].tolist(), 
                                         meta=simple_table.meta)
    #Create another table with the metadata. This is the step that can take long
    #mr_table = metadata_rich_table(td['descriptor'].tolist(), table=simple_table) 
    
    print('')
    print('Table Info after (units, dtype, descriptions)')
    mr_empty_table.info()

    #Convert table to votable. 
    votable = table_to_votable(mr_empty_table)
    
    votable = extend_votable_field(votable, mr_empty_table.meta, attr_name='ucd')
    #votable to .xml file
    #votable.to_xml(os.path.join(out_dirname, basename + '.xml'))
    votable.to_xml(os.path.join(out_dirname, 'dmu22_XID+SPIRE_HELP_BLIND_Matched_MF_onlymeta.xml'))
    

    #Test 3: Only 1 fits table using from_table to votable
    #simple function to create an entire VOTable astropy.io.votable.tree.Table from just 1 single table
    votable = from_table(table)
    votable.to_xml(os.path.join(out_dirname, basename + '_direct.xml'))

