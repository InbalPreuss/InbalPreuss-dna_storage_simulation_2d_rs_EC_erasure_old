import glob
import os

from Bio.SeqIO.QualityIO import FastqGeneralIterator


#################################################################
# @ Function: get_seq_id_offset
# @ Input: seq_id
# @ Description: Calculate the sequence id offset
# @ Return: offset
#################################################################
def get_seq_id_offset(a_seq_id):
    total_seq_index = a_seq_id
    multiply_offset = 9
    number_of_id_digits = 1
    offset = 0

    while (total_seq_index - multiply_offset) > 0:
        offset += (multiply_offset * number_of_id_digits)
        total_seq_index -= multiply_offset
        multiply_offset *= 10
        number_of_id_digits += 1

    offset += ((total_seq_index - 1) * number_of_id_digits)
    return offset


#################################################################
# @ class: FastqHandling
# @ Description: FastqHandling handles the .fastq file and makes
#                a .txt file with the sorted Oligos by the barcode
#################################################################
class FastqHandling:

    #################################################################
    # @ Function: __init__
    # @ Input: number_of_barcode_letters, oligo_length, file_name
    # @ Description: init class FastqHandling
    #################################################################
    def __init__(self, a_number_of_barcode_letters, a_oligo_length, a_file_name):
        fastq_input_file_name = glob.glob(a_file_name + '*')

        if fastq_input_file_name.__len__() == 0:
            raise NameError('The file name ' + a_file_name + ' does not exist')
        elif fastq_input_file_name.__len__() > 1:
            raise NameError(
                'There are too many files with the same name as ' + self.file_name + ', please make sure you '
                                                                                     'have only one file with '
                                                                                     'this specific name ')

        # TODO:
        # import pathlib
        # p = pathlib.Path(r'data/small_data_4_barcode_8_oligo.dna')
        # p.stem, p.suffix
        file_name, file_extension = os.path.splitext(fastq_input_file_name[0])
        self.file_name = file_name
        self.file_full_name = file_name + file_extension
        self.file_extension = file_extension.replace(".", "")
        self.file_full_name_set_ids_output = "Fastq_output/" + file_name + "_set_ids.txt"
        self.file_full_name_sorted_output = "Fastq_output/" + file_name + "_sorted.dna"
        self.number_of_barcode_letters = a_number_of_barcode_letters
        self.oligo_length = a_oligo_length

    #################################################################
    # @ Function: set_oligo_id
    # @ Description: Set new ids to the Oligos,
    #                from 1,2,3,4.... to the amount of oligos in file
    #                and put the ids + Oligo into a .txt file
    #################################################################
    def set_oligo_id(self):
        seq_id = 1
        try:
            os.makedirs(os.path.dirname(self.file_full_name_set_ids_output), exist_ok=True)
            file_name_set_ids = open(self.file_full_name_set_ids_output, "w")
            for (title, sequence, quality) in FastqGeneralIterator(self.file_full_name):
                seq_and_id = str(seq_id) + " " + str(sequence) + "\n"
                file_name_set_ids.write(seq_and_id)
                seq_id += 1
        except FileExistsError:
            print('There is a problem opening the file ' + file_name_set_ids)

        file_name_set_ids.close()

    #################################################################
    # @ Function: sort_oligo
    # @ Description: Sort the Oligos by barcode and insert
    #                all the Oligos into a .txt file
    #################################################################
    def sort_oligo(self):
        try:
            os.makedirs(os.path.dirname(self.file_full_name_sorted_output), exist_ok=True)
            file_name_set_ids = open(self.file_full_name_set_ids_output, "r")
            file_name_sorted = open(self.file_full_name_sorted_output, "w")
        except FileExistsError:
            print('There is a problem opening the file ' + file_name_set_ids + " or " + file_name_sorted)

        # Insert every Oligo sequence, their id and barcode to a list
        id_by_barcode_list = []
        for line in file_name_set_ids:
            # line = id, seq
            id_and_seq = line.split()
            id_by_barcode_list.append([id_and_seq[0], id_and_seq[1][:self.number_of_barcode_letters]])

        # Sort ids by Oligo barcode
        sorted_oligo_by_barcode = sorted(id_by_barcode_list, key=lambda x: x[1])

        # Sort the oligo sequence by barcode and insert to .txt file
        # We do not want to read every time the entire lines to get to the specific line that we want,
        # therefore we calculate the location of the line and jump to that location
        for seq_id, barcode in sorted_oligo_by_barcode:
            # line_offset = prev_line_number * (oligo_length + space + new_line) + seq_id_offset + prev_line_number
            line_offset = (int(seq_id) - 1) * (self.oligo_length + 1 + 1) + get_seq_id_offset(int(seq_id)) + (
                    int(seq_id) - 1)

            # Jumps to the Oligo line and writes the line to the file
            file_name_set_ids.seek(line_offset)
            line = file_name_set_ids.readline()
            file_name_sorted.write(line)

        file_name_set_ids.close()
        file_name_sorted.close()

    def parse_fastq(self):
        self.set_oligo_id()
        self.sort_oligo()

        print("Finished parsing the Fastq file")
        return self.file_full_name_sorted_output
