import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location(
    "prepare_peaks", Path(__file__).resolve().parents[1] / "bin" / "prepare_peaks.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PreparePeaksTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.sizes = self.root / "chrom.sizes"
        self.sizes.write_text("chr1\t1000\nchr2\t500\n")

    def tearDown(self):
        self.temporary.cleanup()

    def run_case(self, text, compressed=False):
        source = self.root / "input.bed"
        if compressed:
            # Detect gzip by content, even when the staged name has no .gz suffix.
            source.write_bytes(gzip.compress(text.encode()))
        else:
            source.write_text(text)
        output = self.root / "peaks.bed"
        report = self.root / "report.json"
        result = MODULE.prepare_peaks(source, self.sizes, output, report)
        self.assertEqual(result, json.loads(report.read_text()))
        return output.read_text(), result

    def test_plain_bed_preserves_intervals(self):
        output, stats = self.run_case("chr1\t10\t20\nchr2\t30\t40\n")
        self.assertEqual(output, "chr1\t10\t20\nchr2\t30\t40\n")
        self.assertEqual(stats["retained_intervals"], 2)

    def test_gzip_headers_and_nonreference_chromosomes(self):
        output, stats = self.run_case(
            "# comment\ntrack name=test\nchr1\t10\t20\nchrUn\t1\t5\n", True)
        self.assertEqual(output, "chr1\t10\t20\n")
        self.assertTrue(stats["compressed_input"])
        self.assertEqual(stats["skipped_nonreference"], 1)

    def test_clipping_adjusts_narrowpeak_summit(self):
        output, stats = self.run_case("chr1\t-10\t20\tp\t10\t.\t4\t5\t6\t15\n")
        fields = output.strip().split("\t")
        self.assertEqual(fields[1:3], ["0", "20"])
        self.assertEqual(fields[9], "5")
        self.assertEqual(stats["clipped_intervals"], 1)

    def test_outside_and_end_clipping(self):
        output, stats = self.run_case("chr1\t990\t1010\nchr2\t600\t700\n")
        self.assertEqual(output, "chr1\t990\t1000\n")
        self.assertEqual(stats["skipped_outside_reference"], 1)

    def test_invalid_or_empty_peaks_fail_without_partial_output(self):
        for text in ["", "chr1\t20\t10\n", "chr1\tx\t20\n", "chrUn\t1\t5\n"]:
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    self.run_case(text)
                self.assertFalse((self.root / "peaks.bed").exists())


if __name__ == "__main__":
    unittest.main()
