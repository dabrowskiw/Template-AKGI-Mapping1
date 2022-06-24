from unittest import TestCase
import pytest

from mapper import Read, Reference, Mapping, read_fasta, map_reads


class MapperTestCase(TestCase):
    @pytest.mark.read
    def test_read_constructor(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("AGTCGTAGTTCAGCCTCGTTAGCTAGGCAATG", read.bases)
        self.assertEqual("Read_0", read.name)

    @pytest.mark.read
    def test_read_str(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("Read_0: AGTCGTAGTTCAGCCTCGTT...", str(read))

    @pytest.mark.read
    def test_read_seed(self):
        read = Read([">Read_0", "AGTCGTAG", "TTCAGCCT", "CGTTAGCT", "AGGCAATG"])
        self.assertEqual("AGTCG", read.get_seed(5))

    @pytest.mark.reference
    def test_reference_constructor(self):
        read = Reference([">Reference", "TTTACTGTGTCCATGGTGTATCCTGTTCCT", "GTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG"])
        self.assertEqual("TTTACTGTGTCCATGGTGTATCCTGTTCCTGTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG", read.bases)
        self.assertEqual("Reference", read.name)

    @pytest.mark.reference
    def test_reference_str(self):
        read = Reference([">Reference", "TTTACTGTGTCCATGGTGTATCCTGTTCCT", "GTTCCATGGCTGTATGGAGGATCTCCAGTATAAGAGAATG"])
        self.assertEqual("Reference: TTTACTGTGTCCATGGTGTA...", str(read))

    @pytest.mark.reference
    def test_reference_get_kmer_positions(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        self.assertEqual([9, 16], ref.get_kmer_positions("TAG"))
        self.assertEqual([0], ref.get_kmer_positions("AGTC"))

    @pytest.mark.reference
    def test_reference_count_mismatches(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        self.assertEqual(4, ref.count_mismatches(read1, 0))
        self.assertEqual(0, ref.count_mismatches(read1, 3))

    @pytest.mark.toplevel
    def test_read_fasta_reference(self):
        references = read_fasta("data/fluA.fasta", Reference.__name__)
        self.assertEqual(1, len(references))
        reference = references[0]
        self.assertEqual(Reference.__name__, reference.__class__.__name__)
        self.assertEqual(len(reference.bases), 2445)

    @pytest.mark.toplevel
    def test_read_fasta_reads(self):
        reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
        self.assertEqual(100, len(reads))
        read = reads[12]
        self.assertEqual(Read.__name__, read.__class__.__name__)
        self.assertEqual(read.bases, "GGCAAAAATAATGAATTTAACTTGTCCTTCATGAAAAAATGCCTGTTTTT")
        
    @pytest.mark.toplevel
    def test_map_reads(self):
        reference = read_fasta("data/fluA.fasta", Reference.__name__)[0]
        reads = read_fasta("data/fluA_reads.fasta", Read.__name__)
        mapping = map_reads(reads, reference, 25, 2)
        mapres = mapping.get_reads_at_position(9)
        self.assertEqual(len(mapres), 1)
        self.assertEqual(mapres[0].name, "Read_95")
        mapres = mapping.get_reads_at_position(109)
        self.assertEqual(len(mapres), 2)
        self.assertIn("Read_25", [x.name for x in mapres], "Expected Read_25 to map at position 109")
        self.assertIn("Read_54", [x.name for x in mapres], "Expected Read_54 to map at position 109")
        mapres = mapping.get_reads_at_position(177)
        self.assertEqual(len(mapres), 0)

    @pytest.mark.mapping
    def test_mapping(self):
        ref = Reference([">ref", "AGTCCTGATTAGCGGTTAGCGAAT"])
        read1 = Read([">read_1", "CCTGAT"])
        read2 = Read([">read_2", "TAGCGGT"])
        mapping = Mapping(ref)
        self.assertEqual([], mapping.get_reads_at_position(0))
        mapping.add_read(read1, 5)
        self.assertEqual([read1], mapping.get_reads_at_position(5))
        self.assertEqual([], mapping.get_reads_at_position(0))
        mapping.add_read(read2, 5)
        self.assertEqual([read1, read2], mapping.get_reads_at_position(5))

