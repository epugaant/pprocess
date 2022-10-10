import os
import yaml    
from astropy.io.fits.connect import read_table_fits    
from pprocess.pprocess_concept import create_directory, table_to_votable, extend_votable_field

empty = True 

path = '/Volumes/EP_DISK2/eRosita/resources/eFEDStables/eFEDS_c001_main_V7.4.fits'
out_dirname = '/Volumes/EP_DISK2/eRosita/votables'
create_directory(out_dirname)
        
filename = os.path.basename(path)
basename = os.path.splitext(filename)[0]
in_dirname = os.path.dirname(path)

print('Input .FITS: {}'.format(os.path.join(in_dirname, filename)))
print('Output VOTABLE: {}'.format(os.path.join(out_dirname, basename+'.xml')))

#-------------
table = read_table_fits(path, memmap=True, astropy_native=True)# more tunnable way to read in table 

print('')
print('Table Info before (descriptions, ucd)')
table.info()

# read yaml with tucd
with open(r'/Volumes/EP_DISK2/eRosita/resources/eFEDStables/eFEDS_c001_main_V7.4UCD.yml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    ucd = yaml.load(file, Loader=yaml.FullLoader)

#Convert table to votable. 
    if empty:
        votable = table_to_votable(table[:0].copy())
    else:
        votable = table_to_votable(table) 

#column descriptions are in table.meta with TCOMM header keyword
votable = extend_votable_field(votable, table.meta, attr_name='description', attr_key='COMM')
#ucds are provided in .yml file
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
