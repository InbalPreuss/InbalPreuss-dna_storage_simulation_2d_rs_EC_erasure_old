import OligoHandling
import FastqHandling

# User parameters: barcode, Oligo file name, Fastq file name
NUMBER_OF_BARCODE_LETTERS = 16
OLIGO_LENGTH = 151
OLIGO_FILE_NAME = 'Oligo_Input'
FASTQ_FILE_NAME = 'Bible4_sample'

if __name__ == "__main__":
    # Parsing Oligo data
    # OligoHandling.OligoHandling(NUMBER_OF_BARCODE_LETTERS, OLIGO_LENGTH, OLIGO_FILE_NAME).parse_oligo()

    # Parsing Fastq data
    FastqHandling.FastqHandling(NUMBER_OF_BARCODE_LETTERS, OLIGO_LENGTH, FASTQ_FILE_NAME).parse_fastq()
