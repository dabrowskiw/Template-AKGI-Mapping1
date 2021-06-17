class Sequence:
    def __init__(self, lines):
        pass

    def __str__(self):
        return ""

    def __repr__(self):
        return ""


class Read(Sequence):
    def get_seed(self, seedlength):
        return None


class Reference(Sequence):
    def __init__(self, lines):
        pass

    def calculate_kmers(self, kmersize):
        pass

    def get_kmer_positions(self, kmer):
        return None

    def count_mismatches(self, read, position):
        return None


class Mapping:
    def __init__(self, reference):
        pass

    def add_read(self, read, position):
        pass

    def get_reads_at_position(self, position):
        return None

    def __str__(self):
        return ""


def read_fasta(fastafile, klassname):
    return None


def map_reads(reads, reference, kmersize, max_mismatches):
    return None


def main():
    reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
    reference = read_fasta("data/fluA.fasta", Reference.__name__)[0]
    mapping = map_reads(reads, reference, 8, 2)
    print(mapping)


if __name__ == "__main__":
    main()
