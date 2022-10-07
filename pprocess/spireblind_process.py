import os
import yaml 
import glob   
from astropy.io import fits, ascii
from astropy.io.fits.connect import read_table_fits    
from pprocess.pprocess_concept import create_directory, append_tables, metadata_rich_table, make_cols_info_table, table_to_votable, extend_votable_field

from astropy import units as u

empty = True

#Read in fits table
fits_list = glob.glob(os.path.join('/Volumes/EP_DISK2/EUProjects/HELP/spire_blind', "*.fits"))
#Location is the same, so take first file
path = fits_list[0]
filename = os.path.basename(path)
basename = os.path.splitext(filename)[0]
in_dirname = os.path.dirname(path)

out_dirname = os.path.join(in_dirname, 'votable')
create_directory(out_dirname)

print('Input .FITS: {}'.format(os.path.join(in_dirname, filename)))
print('Output VOTABLE: {}'.format(os.path.join(out_dirname, basename+'.xml')))

#-------------
_all_fields = ['AKARI-NEP','AKARI-SEP','Bootes','CDFS-SWIRE','COSMOS','EGS','ELAIS-N1','ELAIS-N2','ELAIS-S1','GAMA-09','GAMA-12','GAMA-15','HDF-N','Herschel-Stripe-82','Lockman-SWIRE','NGP','SA13','SGP','SPIRE-NEP','SSDF','xFLS','XMM-13hr','XMM-LSS'] 
fits_list_clean = [x for x in fits_list if "HELP_BLIND" not in x]
fits_list_sort = sorted(fits_list_clean, key=lambda x: [i for i,e in enumerate(_all_fields) if e in x][0])

#serialization fits tables into one table
table = append_tables(fits_list_sort)

#----METADATA--
print('')
print('Table Info before (units, dtype)')
table.info()

table['Bkg_SPIRE_250'].unit = u.mJy/u.beam
table['Bkg_SPIRE_350'].unit = u.mJy/u.beam
table['Bkg_SPIRE_500'].unit = u.mJy/u.beam
table['Sig_conf_SPIRE_250'].unit = u.mJy/u.beam
table['Sig_conf_SPIRE_350'].unit = u.mJy/u.beam
table['Sig_conf_SPIRE_500'].unit = u.mJy/u.beam

#descriptions in an external file
td = ascii.read('/Users/epuga/ESDC/FP7H2020/HELP/spire_blind/column_unit_descriptors.txt', guess=False, 
                    names=['name', 'unit', 'descriptor'], format='no_header', delimiter='\t', fast_reader=False)
    
#Alternative way to add descriptions at the astropy.table.Table level
# if empty:
#     mr_table = metadata_rich_table(td['descriptor'].tolist(), names=table.colnames, 
#                                     dtype=table.dtype, units=make_cols_info_table(table)['unit'].tolist(), 
#                                     meta=table.meta)
# else:
#     #Create another complete table with the metadata. Warning: For large tables this step that can take long
#     mr_table = metadata_rich_table(td['descriptor'].tolist(), table=table) 

# print('')
# print('Table Info after astropy table (units, dtype, descriptions)')
# mr_table.info()

#Convert table to votable. 
if empty:
    votable = table_to_votable(table[:0].copy())
else:
    votable = table_to_votable(table) 

votable = extend_votable_field(votable, td['descriptor'].tolist(), attr_name='description')
#ucds are provided in table.meta with TUCDX header keywords, but it does not work because there are less tucds than fields
#votable = extend_votable_field(votable, table.meta, attr_name='ucd')

print('')
print('Table Info after (descriptions)')
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
