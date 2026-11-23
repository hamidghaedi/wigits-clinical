#!/usr/bin/env python3
"""Unit test for sage_compat_adapter.py on a synthetic, pre-normalised DeepSomatic VCF (no bcftools/Java needed)."""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import sage_compat_adapter as sca  # noqa: E402

DICT = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chrM\tLN:16569\n"
HEADER = ("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
          "##FILTER=<ID=GERMLINE,Description=\"Non somatic variants\">\n"
          "##INFO=<ID=DS_ORIGINAL_ALLELE,Number=1,Type=String,Description=\"orig\">\n"
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">\n"
          "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"\">\n##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"\">\n"
          "##FORMAT=<ID=NAD,Number=R,Type=Integer,Description=\"\">\n##FORMAT=<ID=NDP,Number=1,Type=Integer,Description=\"\">\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKQ0205-DNA\n")
FMT = "GT:AD:DP:VAF:NAD:NDP"
PASS_ROWS = [
    # plain SNV, chr2 first to test dictionary sorting
    "chr2\t50\t.\tA\tG\t40\tPASS\t.\t" + FMT + "\t0/1:30,10:40:0.25:35,0:35",
    # multiallelic 1/2 split into two called alleles
    "chr1\t100\t.\tC\tT\t30\tPASS\tDS_ORIGINAL_ALLELE=chr1|100|C|T,G|1\t" + FMT + "\t1/0:5,8:20:0.4:30,0:30",
    "chr1\t100\t.\tC\tG\t30\tPASS\tDS_ORIGINAL_ALLELE=chr1|100|C|T,G|2\t" + FMT + "\t0/1:5,7:20:0.35:30,0:30",
    # multiallelic 0/1 where the second ALT was not called -> dropped
    "chr1\t200\t.\tA\tT\t30\tPASS\tDS_ORIGINAL_ALLELE=chr1|200|A|T,C|1\t" + FMT + "\t0/1:10,9:20:0.45:40,0:40",
    "chr1\t200\t.\tA\tC\t30\tPASS\tDS_ORIGINAL_ALLELE=chr1|200|A|T,C|2\t" + FMT + "\t0/0:10,1:20:0.05:40,0:40",
    # left-aligned deletion with low ALT support
    "chr1\t300\t.\tAT\tA\t12\tPASS\tDS_ORIGINAL_ALLELE=chr1|301|TT|T|1\t" + FMT + "\t0/1:50,2:52:0.04:60,0:60",
    # symbolic ALT
    "chr1\t400\t.\tA\t<DEL>\t12\tPASS\t.\t" + FMT + "\t0/1:50,2:52:0.04:60,0:60",
]
PANEL_ROWS = [
    # DeepSomatic GERMLINE records carry the somatic GT 0/0; ALT support decides
    "chr1\t500\t.\tG\tA\t50\tGERMLINE\t.\t" + FMT + "\t0/0:20,18:38:0.47:25,12:37",
    "chr1\t600\t.\tG\tA\t50\tGERMLINE\t.\t" + FMT + "\t0/0:40,0:40:0:25,0:25",
    "chr1\t700\t.\tG\tA\t0\tRefCall\t.\t" + FMT + "\t0/0:40,3:43:0.07:25,0:25",
    "chr1\t800\t.\tG\tA\t0\tNoCall\t.\t" + FMT + "\t./.:0,0:0:0:0,0:0",
]


class AdapterBuildTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = Path(self.tmp.name)
        (d / "ref.dict").write_text(DICT)
        (d / "pass.vcf").write_text(HEADER + "\n".join(PASS_ROWS) + "\n")
        (d / "panel.vcf").write_text(HEADER + "\n".join(PANEL_ROWS) + "\n")
        self.d = d

    def tearDown(self):
        self.tmp.cleanup()

    def test_build_dispositions_and_order(self):
        d = self.d
        sca.main(["build", "--pass-vcf", str(d / "pass.vcf"), "--panel-vcf", str(d / "panel.vcf"), "--dict", str(d / "ref.dict"),
                  "--tumor-id", "KQ0205-DNA", "--out-candidates", str(d / "cand.tsv"), "--out-ledger", str(d / "led.tsv"),
                  "--out-json", str(d / "build.json")])
        summary = json.loads((d / "build.json").read_text())
        self.assertEqual(summary["candidates"], 6)  # chr2 SNV, 2 split alleles, 1 called of 200, deletion, germline 500
        self.assertEqual(summary["counts"]["candidate_pass"], 5)
        self.assertEqual(summary["counts"]["candidate_nonpass"], 1)
        self.assertEqual(summary["counts"]["dropped_nonpass_no_alt_support"], 1)
        self.assertEqual(summary["counts"]["dropped_nonpass_no_variant_call"], 2)
        self.assertEqual(summary["counts"]["dropped_alt_not_in_called_gt"], 1)
        self.assertEqual(summary["counts"]["dropped_symbolic_or_missing_alt"], 1)
        self.assertEqual(summary["counts"]["candidate_low_alt_support_lt3"], 1)
        self.assertEqual(summary["counts"]["representation:multiallelic_split"], 4)
        self.assertEqual(summary["counts"]["representation:realigned_or_trimmed"], 1)
        rows = [l.split("\t") for l in (d / "cand.tsv").read_text().splitlines()]
        cols, body = rows[0], rows[1:]
        chrom_pos = [(r[cols.index("chrom")], int(r[cols.index("pos")])) for r in body]
        self.assertEqual(chrom_pos, [("chr1", 100), ("chr1", 100), ("chr1", 200), ("chr1", 300), ("chr1", 500), ("chr2", 50)])
        split = [r for r in body if r[cols.index("pos")] == "100"]
        self.assertEqual({(r[cols.index("alt")], r[cols.index("t_ad_alt")]) for r in split}, {("T", "8"), ("G", "7")})
        self.assertEqual(split[0][cols.index("orig_allele")].split("|")[3], "T,G")
        led = (d / "led.tsv").read_text().splitlines()
        self.assertEqual(len(led) - 1, 11)
        self.assertIn("dropped_alt_not_in_called_gt", "\n".join(led))

    def test_validate_vcf_detects_problems(self):
        d = self.d
        final = ("##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"x\">\n##INFO=<ID=TIER,Number=1,Type=String,Description=\"x\">\n"
                 "##INFO=<ID=TNC,Number=1,Type=String,Description=\"x\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"x\">\n"
                 "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"x\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n"
                 "##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKQ0205-gDNA\tKQ0205-DNA\n"
                 "chr1\t10\t.\tA\tG\t30\tPASS\tTIER=PANEL;TNC=CAG\tGT:AD:DP\t./.:20,0:20\t0/1:5,5:10\n"
                 "chr2\t10\t.\tA\tG\t30\tPASS\tTIER=PANEL;TNC=CAG;REP_C=4\tGT:AD:DP\t./.:20,0:20\t0/1:5,5:10\n")
        (d / "final.vcf").write_text(final)
        v = sca.validate_vcf(str(d / "final.vcf"), str(d / "ref.dict"), "KQ0205-gDNA", "KQ0205-DNA")
        self.assertEqual(v["records"], 2)
        self.assertEqual(v["problems"], ["undeclared INFO: ['REP_C']"])
        bad = final.replace("chr2\t10", "chr1\t5")
        (d / "bad.vcf").write_text(bad)
        v = sca.validate_vcf(str(d / "bad.vcf"), str(d / "ref.dict"), "KQ0205-gDNA", "KQ0205-DNA")
        self.assertTrue(any("out of order" in p for p in v["problems"]))

    def test_tier_parity(self):
        d = self.d
        ours = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n"
                "chr1\t10\t.\tA\tG\t30\tPASS\tTIER=PANEL\tGT\t./.\t0/1\nchr1\t20\t.\tA\tG\t30\tPASS\tTIER=LOW_CONFIDENCE\tGT\t./.\t0/1\n")
        sage = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n"
                "chr1\t10\t.\tA\tG\t30\tPASS\tTIER=PANEL\tGT\t./.\t./.\nchr1\t20\t.\tA\tG\t30\tminTumorQual\tTIER=HIGH_CONFIDENCE\tGT\t./.\t./.\n"
                "chr1\t30\t.\tA\tG\t30\tPASS\tTIER=PANEL\tGT\t./.\t./.\n")
        (d / "ours.vcf").write_text(ours)
        (d / "sage.vcf").write_text(sage)
        p = sca.tier_parity(str(d / "ours.vcf"), str(d / "sage.vcf"))
        self.assertEqual((p["intersection"], p["tier_agree"], p["tier_disagree"]), (2, 1, 1))
        self.assertEqual(p["tier_changes"], {"HIGH_CONFIDENCE->LOW_CONFIDENCE": 1})
        self.assertEqual((p["sage_pass_total"], p["sage_pass_in_adapter"]), (2, 1))


if __name__ == "__main__":
    unittest.main()
