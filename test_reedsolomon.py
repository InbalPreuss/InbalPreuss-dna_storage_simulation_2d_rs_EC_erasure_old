
from reedsolomon.trimer_RS import barcode_rs_encode, barcode_rs_decode

barcode_encoded = barcode_rs_encode('AAAAAAAAAAAT')
print(''.join(barcode_encoded))

barcode_encoded_error = barcode_encoded[:11]
barcode_encoded_error.append('C')
barcode_encoded_error.extend(barcode_encoded[12:])


barcode_decoded = barcode_rs_decode(barcode_encoded_error, verify_only=False)
print(''.join(barcode_decoded))
a = 3