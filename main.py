from oligo_handling import OligoHandling
from fastq_handling import FastqHandling
from oligo_retriever import OligoRetriever


# User parameters: barcode, Oligo file name, Fastq file name


def main(config):
    # Parsing Oligo data
    if config['do_oligo_handling']:
        OligoHandling(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                      file_name=config['OLIGO_FILE_NAME']).parse_oligo()

    # Parsing Fastq data
    if config['do_fastq_handling']:
        file_name_sorted = FastqHandling(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                                         oligo_length=config['OLIGO_LENGTH'],
                                         file_name=config['FASTQ_FILE_NAME']).parse_fastq()

    # Composite Algorithm
    if config['do_retrieved_oligo']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        oligo_length=config['OLIGO_LENGTH'],
                                        k_mer=config['K_MER'])
        shrink_dict = config['shrink_dict']
        oligo_retriever = OligoRetriever(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                                         oligo_length=config['OLIGO_LENGTH'],
                                         oligo_sorted_file_name=config['file_name_sorted'],
                                         algorithm=algorithm,
                                         shrink_dict=shrink_dict,
                                         k_mer=config['K_MER'],
                                         unique_oligo_results_file=config['unique_oligo_results_file'])
        oligo_retriever.run()


if __name__ == "__main__":
    from config import config
    main(config)
