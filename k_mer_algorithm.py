class KMerAlgorithm:
    def __init__(self, algorithm_config, oligo_length: int, k_mer: int):
        self.subset_size = algorithm_config['subset_size']
        self.oligo_length = oligo_length
        self.k_mer = k_mer

    def encode(self):
        pass

    def decode(self, shrunk_payload):
        binary_data = []
        for col_idx in range(int(self.oligo_length / self.k_mer)):
            oligo = [shrinked_oligo[col_idx] for shrinked_oligo in shrunk_payload]
            for key, val in self.k_mer_representative_to_binary.items():
                if set(key) == set(oligo):
                    binary_data.append(val)
                    break
        return binary_data
