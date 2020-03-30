import glob
import os
from file_type import FileType


# packages:
# pandas, xlrd, glob, pip

class OligoHandling:
    def __init__(self, barcode_len: int, file_name: str):
        oligo_input_file_name = glob.glob(file_name + '*')

        if oligo_input_file_name.__len__() == 0:
            raise NameError('The file name ' + file_name + ' does not exist')
        elif oligo_input_file_name.__len__() > 1:
            raise NameError('There are too many files with the same name as ' + file_name + ', please make sure you '
                                                                                            'have only one file with '
                                                                                            'this specific name ')

        file_name, file_extension = os.path.splitext(oligo_input_file_name[0])
        self.file_name = file_name
        self.file_extension = file_extension
        self.barcode_len = barcode_len

    def parse_oligo(self):
        try:
            file_type = FileType(self.file_name, self.file_extension)
            oligo_input_file_df = file_type.read_file()
        except IOError:
            print('Error: File does not open')

        oligo_input_file_df = oligo_input_file_df.rename(
            columns={list(oligo_input_file_df)[0]: 'oligo_id', list(oligo_input_file_df)[1]: 'oligo_dna_sequence'})

        # oligo_input_file_barcode_df = oligo_input_file_df.join(oligo_input_file_df['oligo_dna_sequence'].str[:self.barcode_len].to_frame(), rsuffix='_barcode')
        # oligo_input_file_sort_df = oligo_input_file_barcode_df.sort_values(by=['oligo_dna_sequence_barcode'])
        # oligo_input_file_barcode_gb_df = oligo_input_file_sort_df.groupby('oligo_dna_sequence_barcode')

        print("Finished parsing the Oligo file")


