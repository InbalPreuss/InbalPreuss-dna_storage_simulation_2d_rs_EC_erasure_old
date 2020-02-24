class KMerAlgorithm:
    def __init__(self, algorithm_config, oligo_length: int, k_mer: int):
        self.binary_bits_on = algorithm_config['binary_bits_on']
        self.k_mer_representative_to_binary = algorithm_config['k_mer_representative_to_binary']
        self.oligo_length = oligo_length
        self.k_mer = k_mer


    def encode(self):
        pass

    def decode(self, shrinked_oligos):
        binary_data = []
        for col_idx in range(int(self.oligo_length/self.k_mer)):
            oligo = [shrinked_oligo[col_idx] for shrinked_oligo in shrinked_oligos]
            for key, val in self.k_mer_representative_to_binary.items():
                if set(key) == set(oligo):
                    binary_data.append(val)
                    break
        return binary_data
