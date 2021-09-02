import unittest
import unittest.mock as mock
import sys,os,json, requests
from tempfile import NamedTemporaryFile
import gzip
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.data_access import cs
from Scripts.data_access.db import Variant, CSVariant, CS

class TestCSSummary(unittest.TestCase):
            
    def test_init(self):
        cred_header="trait	region	cs	cs_log10bf	cs_avg_r2	cs_min_r2	low_purity	cs_size	good_cs	cs_id	v	rsid	p	beta	sd	prob	cs_specific_prob	most_severe	gene_most_severe"
        snp_header="trait	region	v	cs	cs_specific_prob	chromosome	position	allele1	allele2	maf	beta	p	se	most_severe	gene_most_severe"
        with NamedTemporaryFile(mode="w+t") as f_snp:
            f_snp.write(snp_header)
            f_snp.seek(0)
            with NamedTemporaryFile(mode="w+t") as f_cred:
                f_cred.write(cred_header)
                f_cred.seek(0)
                
                test = cs.CSSummaryReader(f_snp.name,f_cred.name)

    def test_init_fail(self):
        """load a faulty file, resulting in an error
        """
        #cred header faulty
        cred_header="trait	region	cs"
        snp_header="trait	region	v	cs	cs_specific_prob	chromosome	position	allele1	allele2	maf	beta	p	se	most_severe	gene_most_severe"
        with NamedTemporaryFile(mode="w+t") as f_snp:
            f_snp.write(snp_header)
            f_snp.seek(0)
            with NamedTemporaryFile(mode="w+t") as f_cred:
                f_cred.write(cred_header)
                f_cred.seek(0)
                with self.assertRaises(Exception) as cm:
                    test = cs.CSSummaryReader(f_snp.name,f_cred.name)
        #snp header faulty
        cred_header = "trait	region	cs	cs_log10bf	cs_avg_r2	cs_min_r2	low_purity	cs_size	good_cs	cs_id	v	rsid	p	beta	sd	prob	cs_specific_prob	most_severe	gene_most_severe"
        snp_header="trait	region	v	cs"
        with NamedTemporaryFile(mode="w+t") as f_snp:
            f_snp.write(snp_header)
            f_snp.seek(0)
            with NamedTemporaryFile(mode="w+t") as f_cred:
                f_cred.write(cred_header)
                f_cred.seek(0)
                with self.assertRaises(Exception) as cm:
                    test = cs.CSSummaryReader(f_snp.name,f_cred.name)

    def test_load(self):
        cred_header="region	cs	cs_log10bf	cs_min_r2	cs_size	good_cs	v	cs_specific_prob\n"
        snp_header="region	v	cs	cs_specific_prob\n"
        test_data = CS(
                [CSVariant(
                    Variant("1",1,"A","T"),
                    1.0,
                    1.0
                )],
                Variant("1",1,"A","T"),
                "reg1",
                1,
                1.0,
                1.0,
                1,
                True
            )
        cred_as_string = "\t".join(map(str,[
            test_data.region,
            test_data.number,
            test_data.bayes,
            test_data.min_r2,
            test_data.size,
            test_data.good_cs,
            "1:1:A:T",
            1.0
        ])) + "\n"
        snp_as_string = "\t".join(map(str,[
            test_data.region,
            "1:1:A:T",
            test_data.number,            
            test_data.variants[0].prob
        ]))+ "\n"

        with NamedTemporaryFile(mode="w+t") as f_snp:
            f_snp.write(snp_header)
            f_snp.write(snp_as_string)
            f_snp.seek(0)
            with NamedTemporaryFile(mode="w+t") as f_cred:
                f_cred.write(cred_header)
                f_cred.write(cred_as_string)
                f_cred.seek(0)
                reader = cs.CSSummaryReader(f_snp.name,f_cred.name)
                read_data = reader.get_cs()
                self.assertEqual([test_data],read_data)

class TestCSFull(unittest.TestCase):
    def test_init(self):
        cred_header = "\t".join([
            "cs_log10bf",
            "cs_min_r2",
            "cs",
            "cs_size",
            "region",
            "low_purity"
        ])
        snp_header = "\t".join([
            "v",
            "cs_specific_prob",
            "region",
            "cs",
            "lead_r2"
        ])
        with NamedTemporaryFile(mode="w+b") as cred_temp:
            gzip_cred = gzip.GzipFile(mode="wb",fileobj = cred_temp)
            gzip_cred.write(bytes(cred_header,"utf-8"))
            gzip_cred.close()
            cred_temp.seek(0)
            with NamedTemporaryFile(mode="w+b") as snp_temp:
                gzip_snp = gzip.GzipFile(mode="wb",fileobj = snp_temp)
                gzip_snp.write(bytes(snp_header,"utf-8"))
                gzip_snp.close()
                snp_temp.seek(0)
                test=cs.CSFullReader(snp_temp.name,cred_temp.name)

    def test_init_fail(self):
        """load a faulty file, resulting in an error
        """
        #cred header faulty
        cred_header = "\t".join([
            "cs_log10bf",
            "cs_min_r2",
            "cs"
        ])
        snp_header = "\t".join([
            "v",
            "cs_specific_prob",
            "region",
            "cs",
            "lead_r2"
        ])
        with NamedTemporaryFile(mode="w+b") as f_snp:
            gzip_snp = gzip.GzipFile(mode="wb",fileobj = f_snp)
            gzip_snp.write(bytes(snp_header,"utf-8"))
            gzip_snp.close()
            f_snp.seek(0)
            with NamedTemporaryFile(mode="w+b") as f_cred:
                gzip_cred = gzip.GzipFile(mode="wb",fileobj = f_cred)
                gzip_cred.write(bytes(cred_header,"utf-8"))
                gzip_cred.close()
                f_cred.seek(0)
                with self.assertRaises(Exception) as cm:
                    test = cs.CSFullReader(f_snp.name,f_cred.name)
        #snp header faulty
        cred_header = "\t".join([
            "cs_log10bf",
            "cs_min_r2",
            "cs",
            "cs_size",
            "region",
            "low_purity"
        ])
        snp_header = "\t".join([
            "v",
            "cs_specific_prob"
        ])
        with NamedTemporaryFile(mode="w+b") as f_snp:
            gzip_snp = gzip.GzipFile(mode="wb",fileobj = f_snp)
            gzip_snp.write(bytes(snp_header,"utf-8"))
            gzip_snp.close()
            f_snp.seek(0)
            with NamedTemporaryFile(mode="w+b") as f_cred:
                gzip_cred = gzip.GzipFile(mode="wb",fileobj = f_cred)
                gzip_cred.write(bytes(cred_header,"utf-8"))
                gzip_cred.close()
                f_cred.seek(0)
                with self.assertRaises(Exception) as cm:
                    test = cs.CSFullReader(f_snp.name,f_cred.name)

    def test_load(self):
        cred_header = "\t".join([
            "cs_log10bf",
            "cs_min_r2",
            "cs",
            "cs_size",
            "region",
            "low_purity"
        ])+"\n"
        snp_header = "\t".join([
            "v",
            "cs_specific_prob",
            "region",
            "cs",
            "lead_r2"
        ])+"\n"
        #create data
        test_data = CS(
            [CSVariant(
                Variant("1",1,"A","T"),
                1.0,
                1.0
            )],
            Variant("1",1,"A","T"),
            "reg1",
            1,
            1.0,
            1.0,
            1,
            True
        )
        cred_as_string = "\t".join(
            map(str,[
                test_data.bayes,
                test_data.min_r2,
                test_data.number,
                test_data.size,
                test_data.region,
                not test_data.good_cs
            ])
        )+"\n"
        snp_as_string = "\t".join(
            map(str,[
                "1:1:A:T",
                test_data.variants[0].prob,
                test_data.region,
                test_data.number,
                test_data.variants[0].lead_r2
            ])
        )+"\n"
        with NamedTemporaryFile(mode="w+b") as cred_temp:
            gzip_cred = gzip.GzipFile(mode="wb",fileobj = cred_temp)
            gzip_cred.write(bytes(cred_header,"utf-8"))
            gzip_cred.write(bytes(cred_as_string,"utf-8"))
            gzip_cred.close()
            cred_temp.seek(0)
            with NamedTemporaryFile(mode="w+b") as snp_temp:
                gzip_snp = gzip.GzipFile(mode="wb",fileobj = snp_temp)
                gzip_snp.write(bytes(snp_header,"utf-8"))
                gzip_snp.write(bytes(snp_as_string,"utf-8"))
                gzip_snp.close()
                snp_temp.seek(0)
                reader=cs.CSFullReader(snp_temp.name,cred_temp.name)
                validation = reader.get_cs()
                self.assertEqual([test_data],validation)
