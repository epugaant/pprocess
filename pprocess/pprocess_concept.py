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
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from astropy.io.fits.connect import read_table_fits
from astropy.io.votable.table import from_table

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

    #Deconstruct into numpy_array and a separate table with columns for testing purposes
    table = read_table_fits(path, memmap=True, astropy_native=True)
    tinfo = make_cols_info_table(table)
    bare_data = table.as_array() # produces a structured array (records array)
    
    #Construct metadata-rich table
    print(split_bytes_dtype_list(tinfo['dtype'].tolist()))
    mr_table = Table(rows=table.as_array(), names=tinfo['name'].tolist(), 
                     dtype=split_bytes_dtype_list(tinfo['dtype'].tolist()), meta=table.meta, 
                     descriptions=tinfo['description'].tolist(), units=tinfo['unit'].tolist())
    


    fits_list = glob.glob(os.path.join('/Volumes/EP_DISK2/EUProjects/HELP/spire_blind', "*.fits"))
    _all_fields = ['AKARI-NEP','AKARI-SEP','Bootes','CDFS-SWIRE','COSMOS','EGS','ELAIS-N1','ELAIS-N2','ELAIS-S1','GAMA-09','GAMA-12','GAMA-15','HDF-N','Herschel-Stripe-82','Lockman-SWIRE','NGP','SA13','SGP','SPIRE-NEP','SSDF','xFLS','XMM-13hr','XMM-LSS'] 
    fits_list_clean = [x for x in fits_list if "HELP_BLIND" not in x]
    fits_list_sort = sorted(fits_list_clean, key=lambda x: [i for i,e in enumerate(_all_fields) if e in x][0])
    mr_table = append_tables(fits_list_sort)
        
    #Convert to votable
    votable = table_to_votable(mr_table)
    #votable to .xml file
    #votable.to_xml(os.path.join(out_dirname, basename + '.xml'))
    votable.to_xml(os.path.join(out_dirname, 'dmu22_XID+SPIRE_HELP_BLIND_Matched_MF.xml'))
    

    #Only 1 table to xml
    #there is also a function to create an entire VOTable file with just a single table
    votable = from_table(table)
    votable.to_xml(os.path.join(out_dirname, basename + '_direct.xml'))
    #TODO: serialize fit tables


    #fits read to write for xml (no intermediate steps)
    tbl = Table.read(path)
    Table.read.list_formats()
    tbl.write(os.path.join(out_dirname, basename + '_read_write.xml'), format='votable')
    tbl.info()

