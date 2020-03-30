import time
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
        OligoHandling(barcode_len=config['barcode_len'],
                      file_name=config['oligo_file_name']).parse_oligo()

    if config['write_text_to_binary']:
        text_file_to_binary = TextFileToBinaryFile(input_file=config['input_text_file'],
                                                   output_file=config['binary_file_name'],
                                                   payload_len=config['payload_len'],
                                                   bits_per_z=config['algorithm_config']['bits_per_z'],
                                                   k_mer=config['k_mer'])
        text_file_to_binary.run()
    # Encode
    if config['do_encode']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        payload_len=config['payload_len'],
                                        k_mer=config['k_mer'])
        shrink_dict = config['shrink_dict']
        encoder = Encoder(barcode_len=config['barcode_len'],
                          barcode_rs_len=config['barcode_rs_len'],
                          payload_len=config['payload_len'],
                          payload_rs_len=config['payload_rs_len'],
                          binary_file_name=config['binary_file_name'],
                          algorithm=algorithm,
                          shrink_dict=shrink_dict,
                          k_mer=config['k_mer'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          binary_to_z=config['algorithm_config']['binary_to_z'],
                          bits_per_z=config['algorithm_config']['bits_per_z'],
                          results_file=config['encoder_results_file'])
        encoder.run()

    # Synthesize
    if config['do_synthesize']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        payload_len=config['payload_len'],
                                        k_mer=config['k_mer'])
        synthesizer = Synthesizer(input_file=config['encoder_results_file'],
                                  results_file=config['synthesis_results_file'],
                                  synthesis_config=config['synthesis'],
                                  barcode_total_len=config['barcode_total_len'],
                                  subset_size=config['algorithm_config']['subset_size'],
                                  algorithm=algorithm,
                                  k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                                  k_mer_to_dna=config['algorithm_config']['k_mer_to_dna'],
                                  mode=config['mode'])
        synthesizer.synthesize()

    # Parsing Fastq data
    if config['do_fastq_handling']:
        file_name_sorted = FastqHandling(barcode_len=config['barcode_len'],
                                         payload_len=config['payload_len'],
                                         file_name=config['fastq_file_name']).parse_fastq()
    # Decode
    if config['do_decode']:
        algorithm = config['algorithm'](algorithm_config=config['algorithm_config'],
                                        payload_len=config['payload_len'],
                                        k_mer=config['k_mer'])
        shrink_dict = config['shrink_dict']
        decoder = Decoder(barcode_len=config['barcode_len'],
                          barcode_total_len=config['barcode_total_len'],
                          payload_len=config['payload_len'],
                          payload_total_len=config['payload_total_len'],
                          input_file=config['synthesis_results_file'],
                          algorithm=algorithm,
                          shrink_dict=shrink_dict,
                          k_mer=config['k_mer'],
                          k_mer_representative_to_z=config['algorithm_config']['k_mer_representative_to_z'],
                          z_to_binary=config['algorithm_config']['z_to_binary'],
                          results_file=config['decoder_results_file'])
        decoder.run()

    if config['decoder_results_to_binary']:
        decoder_results_to_binary = DecoderResultToBinary(input_file=config['decoder_results_file'],
                                                          output_file=config['binary_results_file'],
                                                          barcode_len=config['barcode_len'])
        decoder_results_to_binary.run()

    if config['binary_results_to_text']:
        binary_results_to_text = BinaryResultToText(input_file=config['binary_results_file'],
                                                    output_file=config['text_results_file'],
                                                    barcode_len=config['barcode_len'],
                                                    oligo_len_binary=config['oligo_len_binary'])
        binary_results_to_text.run()


if __name__ == "__main__":
    from config import config
    main(config)
