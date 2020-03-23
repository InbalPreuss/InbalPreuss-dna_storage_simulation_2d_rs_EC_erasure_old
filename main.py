from pathlib import Path

from oligo_handling import OligoHandling
from fastq_handling import FastqHandling
from text_handling import TextFileToBinaryFile, DecoderResultToBinary, BinaryResultToText
from decoder import Decoder
from encoder import Encoder
from mock_synthesizer import Synthesizer


# User parameters: barcode, Oligo file name, Fastq file name


def main(config):
    # Parsing Oligo data
    if config['do_oligo_handling']:
        OligoHandling(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                      file_name=config['OLIGO_FILE_NAME']).parse_oligo()

    if config['write_text_to_binary']:
        text_file_to_binary = TextFileToBinaryFile(input_file=config['input_text_file'],
                                                   output_file=config['binary_file_name'],
                                                   oligo_length=config['oligo_length'],
                                                   bits_per_z=config['algorithm_config']['bits_per_z'],
                                                   k_mer=config['K_MER'])
        text_file_to_binary.run()
    # Encode
    if config['do_encode']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        oligo_length=config['oligo_length'],
                                        k_mer=config['K_MER'])
        shrink_dict = config['shrink_dict']
        encoder = Encoder(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                          oligo_length=config['oligo_length'],
                          binary_file_name=config['binary_file_name'],
                          algorithm=algorithm,
                          shrink_dict=shrink_dict,
                          k_mer=config['K_MER'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          binary_to_z=config['algorithm_config']['binary_to_z'],
                          bits_per_z=config['algorithm_config']['bits_per_z'],
                          results_file=config['encoder_results_file'])
        encoder.run()

    # Synthesize
    if config['do_synthesize']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        oligo_length=config['oligo_length'],
                                        k_mer=config['K_MER'])
        synthesizer = Synthesizer(input_file=config['encoder_results_file'],
                                  results_file=config['synthesis_results_file'],
                                  synthesis_config=config['synthesis'],
                                  number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                                  subset_size=config['algorithm_config']['subset_size'],
                                  algorithm=algorithm,
                                  k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                                  k_mer_to_dna=config['algorithm_config']['k_mer_to_dna'])
        synthesizer.synthesize()

    # Parsing Fastq data
    if config['do_fastq_handling']:
        file_name_sorted = FastqHandling(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                                         oligo_length=config['oligo_length'],
                                         file_name=config['FASTQ_FILE_NAME']).parse_fastq()
    # Decode
    if config['do_decode']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        oligo_length=config['oligo_length'],
                                        k_mer=config['K_MER'])
        shrink_dict = config['shrink_dict']
        decoder = Decoder(number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                          oligo_length=config['oligo_length'],
                          input_file=config['synthesis_results_file'],
                          algorithm=algorithm,
                          shrink_dict=shrink_dict,
                          k_mer=config['K_MER'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          z_to_binary=config['algorithm_config']['z_to_binary'],
                          results_file=config['decoder_results_file'])
        decoder.run()

    if config['decoder_results_to_binary']:
        decoder_results_to_binary = DecoderResultToBinary(input_file=config['decoder_results_file'],
                                                          output_file=config['binary_results_file'],
                                                          number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'])
        decoder_results_to_binary.run()

    if config['binary_results_to_text']:
        binary_results_to_text = BinaryResultToText(input_file=config['binary_results_file'],
                                                    output_file=config['text_results_file'],
                                                    number_of_barcode_letters=config['NUMBER_OF_BARCODE_LETTERS'],
                                                    oligo_len_binary=config['oligo_len_binary'])
        binary_results_to_text.run()


if __name__ == "__main__":
    from config import config

    main(config)
