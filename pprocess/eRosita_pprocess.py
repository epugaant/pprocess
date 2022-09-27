import os
import yaml    
from astropy.io.fits.connect import read_table_fits    
from pprocess.pprocess_concept import create_directory, table_to_votable, extend_votable_field

path = '/Volumes/EP_DISK2/eRosita/eFEDS_c001_main_V7.4.fits'
out_dirname = '/Users/epuga/ESDC/eRosita/votables/help'
create_directory(out_dirname)
        
filename = os.path.basename(path)
basename = os.path.splitext(filename)[0]
in_dirname = os.path.dirname(path)

print('Input .FITS: {}'.format(os.path.join(in_dirname, filename)))
print('Output VOTABLE: {}'.format(os.path.join(out_dirname, basename+'.xml')))

#-------------

#Test 1: Metadata-rich fits table to votable
#Deconstruct into numpy_array and a separate table with columns for testing purposes
table = read_table_fits(path, memmap=True, astropy_native=True)# more tunnable way to read in table 

print('')
print('Table Info before (descriptions, ucd)')
table.info()

# read yaml with tucd
with open(r'/Volumes/EP_DISK2/eRosita/package/eFEDStables/eFEDS_c001_main_V7.4UCD.yml') as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    ucd = yaml.load(file, Loader=yaml.FullLoader)


#Convert table to votable. 
votable = table_to_votable(table)

#column descriptions are in table.meta with TCOMM
votable = extend_votable_field(votable, table.meta, attr_name='description', attr_key='COMM')
#ucds are provided in .yml file
votable = extend_votable_field(votable, ucd, attr_name='ucd')

print('')
print('Table Info after (ucds, descriptions)')
t = votable.get_first_table().to_table()# variable to hold temporary table, cannot do info directly 
print(t.info())

#votable.to_xml(os.path.join(out_dirname, basename + '.xml'))
votable.to_xml(os.path.join(out_dirname, 'eFEDS_c001_main_V7.4.xml'))
