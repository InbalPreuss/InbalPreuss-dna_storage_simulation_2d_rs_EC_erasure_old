
from reedsolomon.trimer_RS import barcode_rs_encode, barcode_rs_decode
from reedsolomon.trimer_RS import rs4096_encode, rs4096_decode


def test_reed_solomon_z_encode_decode():
    z_list = ['Z1' for i in range(120)]
    z_encoded = rs4096_encode(z_list)

    z_encoded_error = z_encoded[:119]
    z_encoded_error.append('Z2')
    z_encoded_error.extend(z_encoded[120:])
    z_list_decoded = rs4096_decode(z_encoded_error, verify_only=False)


def test_reed_solomon_barcode_encode_decode():
    barcode_encoded = barcode_rs_encode('AAAAAAAAAAAT')
    print(''.join(barcode_encoded))

    barcode_encoded_error = barcode_encoded[:11]
    barcode_encoded_error.append('C')
    barcode_encoded_error.extend(barcode_encoded[12:])

    barcode_decoded = barcode_rs_decode(barcode_encoded_error, verify_only=False)
    print(''.join(barcode_decoded))


if __name__ == '__main__':
    test_reed_solomon_z_encode_decode()
    test_reed_solomon_barcode_encode_decode()
