from oligo_handling import OligoHandling
from fastq_handling import FastqHandling
from decoder import Decoder
from encoder import Encoder
from mock_synthesizer import Synthesizer

# User parameters: barcode, Oligo file name, Fastq file name


def main(config):
    # Parsing Oligo data
    if config['do_oligo_handling']:
        OligoHandling(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                      file_name=config['OLIGO_FILE_NAME']).parse_oligo()

    # Encode
    if config['do_encode']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        oligo_length=config['OLIGO_LENGTH'],
                                        k_mer=config['K_MER'])
        shrink_dict = config['shrink_dict']
        encoder = Encoder(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                          oligo_length=config['OLIGO_LENGTH'],
                          binary_file_name=config['binary_file_name'],
                          algorithm=algorithm,
                          shrink_dict=shrink_dict,
                          k_mer=config['K_MER'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          z_to_binary=config['algorithm_config']['z_to_binary'],
                          results_file=config['z_results_file'])
        encoder.run()

    # Synthesize
    if config['do_synthesize']:
        synthesizer = Synthesizer(input_file=config['z_results_file'],
                                  results_file=config['synthesis_results_file'],
                                  number_of_oligos_per_barcode=config['synthesis']['number_of_oligos_per_barcode'],
                                  letter_error_ratio=config['synthesis']['letter_error_ratio'],
                                  k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],)
        synthesizer.synthesize()

    # Parsing Fastq data
    if config['do_fastq_handling']:
        file_name_sorted = FastqHandling(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                                         oligo_length=config['OLIGO_LENGTH'],
                                         file_name=config['FASTQ_FILE_NAME']).parse_fastq()
    # Decode
    if config['do_decode']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        oligo_length=config['OLIGO_LENGTH'],
                                        k_mer=config['K_MER'])
        shrink_dict = config['shrink_dict']
        decoder = Decoder(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                          oligo_length=config['OLIGO_LENGTH'],
                          oligo_sorted_file_name=config['file_name_sorted'],
                          algorithm=algorithm,
                          shrink_dict=shrink_dict,
                          k_mer=config['K_MER'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          z_to_binary=config['algorithm_config']['z_to_binary'],
                          results_file=config['binary_results_file'])
        decoder.run()


if __name__ == "__main__":
    from config import config

    main(config)
