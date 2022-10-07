import os
import glob   
from typing import OrderedDict

from astropy.table import Table  
from pprocess.pprocess_concept import create_directory, make_cols_info_table, table_to_votable, extend_votable_field

empty = True

#Read in fits table
fits_list_init = glob.glob(os.path.join('/Volumes/EP_DISK2/EUProjects/HELP/herschelhelp', "*.fits"))
with open('/Users/epuga/ESDC/pprocess/pprocess/herschelhelp_files_url.txt', 'rt') as f:
    fits_list = [url.strip() for url in f.readlines()]

#ucds are in a global external ecsv file
tucd = Table.read('/Volumes/EP_DISK2/EUProjects/HELP/herschelhelp/herschelhelp_allfields_ucd.ecsv')  
ucd = OrderedDict(zip(tucd.columns[0], tucd.columns[1]))

#Process each of the files to turn them into metadata-rich votables
for path in fits_list[15:]:
    filename = os.path.basename(path)
    basename = os.path.splitext(filename)[0]
    in_dirname = os.path.dirname(path)
    out_dirname = os.path.join('/Volumes/EP_DISK2/EUProjects/HELP/herschelhelp', 'votable')
    create_directory(out_dirname)

    print('Input .FITS: {}'.format(os.path.join(in_dirname, filename)))
    print('Output VOTABLE: {}'.format(os.path.join(out_dirname, basename+'.xml')))

    #-------------
    
    table = Table.read(path, memmap=empty)
    tinfo = make_cols_info_table(table)

    #----METADATA--
    print('')
    print('Table Info before (ucd, descriptions)')
    table.info()

    #Convert table to votable. 
    if empty:
        votable = table_to_votable(table[:0].copy())
    else:
        votable = table_to_votable(table) 

    #ucds in an external file
    votable = extend_votable_field(votable, tinfo['description'].tolist(), attr_name='description')
    #ucds are provided in external table
    votable = extend_votable_field(votable, ucd, attr_name='ucd')

    print('')
    print('Table Info after (ucds, descriptions)')
    t = votable.get_first_table().to_table()# variable to hold temporary table, cannot do info directly 
    print(t.info())
    
    #-------------
    #save the final votable in an xml file changing extension
    if empty:
        votable.to_xml(os.path.join(out_dirname, basename + '_nodata.xml'))
    else:
        votable.to_xml(os.path.join(out_dirname, basename + '.xml')) 

# In an empty votable .xml, it is still necessary to include the fields after the last FIELD element
# <DATA>
# <TABLEDATA>
# </TABLEDATA>
# </DATA>
