class CompositeAlgorithm:
    def __init__(self, ratio: list = 1):
        self.ratio = ratio

    def encode(self):
        pass

    def decode(self, barcode_histogram):
        unique_oligo = ''
        for hist in barcode_histogram:
            # TODO: make it for more letters (ratio) now it works only for the most common letter (max)
            most_common_letter = hist.most_common()
            unique_oligo += most_common_letter[0][0]
        return unique_oligo
