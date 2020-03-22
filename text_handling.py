from textwrap import wrap

import numpy as np


class TextFileToBinaryFile:
    def __init__(self, input_file: str, output_file: str, oligo_length: int, bits_per_z: int, k_mer: int):
        self.input_file = input_file
        self.output_file = output_file
        self.oligo_length = oligo_length
        self.bits_per_z = bits_per_z
        self.k_mer = k_mer

    def run(self):
        with open(self.input_file, 'r') as input_file, open(self.output_file, 'w') as output_file:
            oligo_len_binary = int(self.oligo_length / self.k_mer * self.bits_per_z)
            accumulation = ''
            for line in input_file:
                text_data = line
                binary_data = text_to_bits(text_data)
                # binary_data_padded, oligo_len_binary = self.transform_text_to_binary_string(text_data=text_data)
                accumulation += binary_data
                while len(accumulation) >= oligo_len_binary:
                    to_write = accumulation[:oligo_len_binary]
                    accumulation = accumulation[oligo_len_binary:]
                    output_file.write(to_write + '\n')
            z_fill = 0
            if len(accumulation) > 0:
                binary_data_padded, z_fill = self.transform_text_to_binary_string(binary_data=accumulation)
                output_file.write(binary_data_padded + '\n')

            z_fill_text = "{0:b}".format(z_fill).rjust(oligo_len_binary, '0')
            output_file.write(z_fill_text + '\n')
            a = 3

    def transform_text_to_binary_string(self, binary_data: str):
        oligo_len_binary = int(self.oligo_length / self.k_mer * self.bits_per_z)
        binary_data_len = len(binary_data)

        number_of_binary_oligos = np.ceil(binary_data_len / oligo_len_binary)
        total_binary_len = int(number_of_binary_oligos * oligo_len_binary)

        binary_data_padded = binary_data.ljust(total_binary_len, '0')
        # binary_data_padded = binary_data.zfill(total_binary_len)
        z_fill = total_binary_len - binary_data_len

        # binary_data_padded = "{0:b}".format(z_fill).rjust(oligo_len_binary, '0') + binary_data_padded
        # binary_data_padded += "{0:b}".format(z_fill).rjust(oligo_len_binary, '0')
        return binary_data_padded, z_fill


class DecoderResultToBinary:
    def __init__(self, input_file: str, output_file: str, number_of_barcode_letters: int):
        self.input_file = input_file
        self.output_file = output_file
        self.number_of_barcode_letters = number_of_barcode_letters

    def run(self):
        with open(self.input_file, 'r') as input_file, open(self.output_file, 'w') as output_file:
            for idx, line in enumerate(input_file):
                barcode_and_payload = line.strip()
                barcode, payload = barcode_and_payload[:self.number_of_barcode_letters], barcode_and_payload[
                                                        self.number_of_barcode_letters:]
                output_file.write(payload + '\n')


class BinaryResultToText:
    def __init__(self, input_file: str, output_file: str, number_of_barcode_letters: int):
        self.input_file = input_file
        self.output_file = output_file
        self.number_of_barcode_letters = number_of_barcode_letters

    def run(self):
        with open(self.input_file, 'r') as input_file, open(self.output_file, 'w') as output_file:
            last_two_payloads = ['', '']
            accumulation = ''
            for idx, line in enumerate(input_file):
                payload = line.strip()

                last_two_payloads[0] = last_two_payloads[1]
                last_two_payloads[1] = payload
                accumulation += payload
                while len(accumulation) >= 32:
                    bits = accumulation[:32]
                    accumulation = accumulation[32:]
                    text = text_from_bits(bits)
                    output_file.write(text)

            z_fill = int(last_two_payloads[1], 2)
            payload = last_two_payloads[0]
            payload = payload[z_fill:]
            text = text_from_bits(payload)
            output_file.write(text)


def text_to_bits(text: str, encoding='utf-8', errors='surrogatepass') -> str:
    bits = bin(int.from_bytes(text.encode(encoding, errors), 'big'))[2:]
    return bits.zfill(8 * ((len(bits) + 7) // 8))


def text_from_bits(bits: str, encoding='utf-8', errors='surrogatepass') -> str:
    n = int(bits, 2)
    return n.to_bytes((n.bit_length() + 7) // 8, 'big').decode(encoding, errors) or '\0'