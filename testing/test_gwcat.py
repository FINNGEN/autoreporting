import unittest
import unittest.mock as mock
import sys,os,json, requests
import pandas as pd
from typing import List
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts.data_access import gwcatalog_api
from Scripts.data_access.db import AlleleDB, Location, VariantData, Variant, RsidVar, Rsid
from Scripts.autoreporting_utils import Region
import itertools
from io import StringIO

class TestRequest(unittest.TestCase):
    """Test request wrappers
    """

    @mock.patch("Scripts.data_access.gwcatalog_api.time.sleep")
    def test_get(self,mock_time):
        """
        Test try_request
        Test cases: 
            Test a normal request with normal results, check that it return the orrect object
            Test a normal request that returns 404, check that it return nothing
            Test a normal request that returns 400, check that it return nothing
            Test a normal request that returns error code 500 etc repeatedly, 
            check that it calls the url as many times that is necessary, and returns None
            Test a request that causes requests.get to except
        """
        #200
        response=mock.Mock()
        response.status_code=200
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response):
            retval=gwcatalog_api.try_request("GET",url, params=params)
        self.assertTrue(type(retval)!= type(None))
        #404
        response=mock.Mock()
        response.status_code=404
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("GET",url, params=params)
        #400
        response=mock.Mock()
        response.status_code=400
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("GET",url, params=params)
        #500
        response=mock.Mock()
        response.status_code=500
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response) as mock_get:
            with mock.patch("Scripts.data_access.gwcatalog_api.print"):
                with self.assertRaises(gwcatalog_api.ResponseFailure):
                    retval=gwcatalog_api.try_request("GET",url, params=params)
        self.assertEqual(5,mock_get.call_count)
        #exception during requests.get
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",side_effect=TimeoutError("Timeout")) as mock_get:
            with mock.patch("Scripts.data_access.gwcatalog_api.print") as mock_print:
                with self.assertRaises(gwcatalog_api.ResponseFailure):
                    retval=gwcatalog_api.try_request("GET",url, params=params)
        mock_print.assert_called_with("Request caused an exception:{}".format(TimeoutError("Timeout")))

    @mock.patch("Scripts.data_access.gwcatalog_api.time.sleep")
    def test_post(self,mock_time):
        """
        Test try_request POST functionality
        Test cases:
            Test a normal request with normal results
            Test a normal request that returns 404, should 
            Test a normal request that returns 400
            Test a normal request that returns 500
            Test a request that causes an exception
        """
        #200
        response=mock.Mock()
        response.status_code=200
        url="url"
        headers="header"
        data="data"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response):
            retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        self.assertTrue(type(retval)!= type(None))
        #404
        response=mock.Mock()
        response.status_code=404
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        #400
        response=mock.Mock()
        response.status_code=400
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        #500
        response=mock.Mock()
        response.status_code=500
        url="url"
        params="params"
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",return_value=response) as mock_get:
            with self.assertRaises(gwcatalog_api.ResponseFailure):
                with mock.patch("Scripts.data_access.gwcatalog_api.print"):
                    retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        #exception during requests.post
        with mock.patch("Scripts.data_access.gwcatalog_api.requests.request",side_effect=TimeoutError("Timeout")) as mock_get:
            with self.assertRaises(gwcatalog_api.ResponseFailure):
                with mock.patch("Scripts.data_access.gwcatalog_api.print") as mock_print:
                    retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        pass 


class MockAlleleDB(AlleleDB):
    """Mock class for functions that use alleledb
    """
    def get_alleles(self,positions:List[Location]):
        return [
            VariantData(
                Variant(
                "1",
                1,
                "A",
                "C"),
                [],
                12345,
            ),
            VariantData(
                Variant("2",
                10,
                "A",
                "C"),
                ["G"],
                54321,
            )
        ]

class TestGwcat(unittest.TestCase):
    """Test gwcatalog api functions not relevant to a single db
    """
    

    def test_get_trait_name(self):
        """
        Test get_trait_name
        Test cases:
            Call with a trait that 'does exist' in gwascatalog
            Call with a trait that 'does not exist' in gwascatalog
        """
        #does exist
        trait="trait1"
        trait_name="trait_name"
        tryurl="https://www.ebi.ac.uk/gwas/rest/api/efoTraits/{}".format(trait)
        response_mock=mock.Mock()
        response_mock.json=lambda:{"trait":trait_name}
        with mock.patch("Scripts.data_access.gwcatalog_api.try_request",return_value = response_mock) as mock_try:
            return_value=gwcatalog_api.get_trait_name(trait)
        self.assertEqual(return_value,trait_name)
        mock_try.assert_called_with("GET",url=tryurl)
        #does not exist, raises exception gwcatalog_api.ResourceNotFound
        with mock.patch("Scripts.data_access.gwcatalog_api.try_request",side_effect = gwcatalog_api.ResourceNotFound):
            with mock.patch("Scripts.data_access.gwcatalog_api.print"):
                return_value=gwcatalog_api.get_trait_name(trait)
        self.assertEqual("NA",return_value)
        #TODO: handle unhandled exception

    def test_parse_efo(self):
        """
        Test parse_efo
        Test cases:
            Valid efo code
            Invalid type
        """
        #valid code
        efo = "asd/EFO_CODE"
        validate="EFO_CODE"
        retval=gwcatalog_api.parse_efo(efo)
        self.assertEqual(retval,validate)
        #invalid type, here int
        efo=123456
        with mock.patch("Scripts.data_access.gwcatalog_api.print") as mock_print:
            retval=gwcatalog_api.parse_efo(efo)
        mock_print.assert_called_with("Invalid EFO code:{}".format(efo))
        self.assertEqual("NAN",retval)

    def test_split_traits(self):
        """
        Test stplit_traits
        Test cases:
            A dataframe that works
            A dataframe with no proper columns
        """
        #correct dataframe
        variants=["1","2","3","4"]
        traits=["A","B,C","D,E","F"]
        trait_uris=["1","2,3","4,5","6"]
        data=pd.DataFrame({"variant":variants,"MAPPED_TRAIT":traits,"MAPPED_TRAIT_URI":trait_uris})
        validationvar=["1","2","2","3","3","4"]
        validationtraits=["A","B","C","D","E","F"]
        validationuris=["1","2","3","4","5","6"]
        validationdata=pd.DataFrame({"variant":validationvar,"MAPPED_TRAIT":validationtraits,"MAPPED_TRAIT_URI":validationuris})
        data2=gwcatalog_api.split_traits(data)
        self.assertTrue(validationdata.equals(data2))
        #faulty dataframe: no correct columns
        data=pd.DataFrame({"variant":variants,"MT":traits,"MTU":trait_uris})
        data2=gwcatalog_api.split_traits(data)
        self.assertEqual(type(data2),type(None))

    def test_add_alleles(self):
        """Test that add_alleles function works
        Three data items: exists in alleledb (1:1:A:C), is multiallelic (2:10:A:C,G), is not in alleledb (3:100:whatever)
        """
        
        data={
            "chrom":["1","2","3"],
            "pos":[1,10,100],
            "rsid":[12345,54321,90001],
            "datacol":["a","b","c"]
        }
        data=pd.DataFrame(data)
        alleledb = MockAlleleDB()
        output = gwcatalog_api.add_alleles(data,alleledb)
        validationdata = {
            "chrom":["1"],
            "pos":[1],
            "rsid":[12345],
            "datacol":["a"],
            "ref":["A"],
            "alt":["C"]
        } 
        validationdata = pd.DataFrame(validationdata)
        self.assertTrue(output.equals(validationdata))

    def test_resolve_alleles(self):
        alldb = MockAlleleDB()
        variantdata = alldb.get_alleles([])
        rsid_data=[
            Rsid(Location("1",1),12345),
            Rsid(Location("2",10),54321),
            Rsid(Location("3",100),90001)
        ]
        out = gwcatalog_api._resolve_alleles(rsid_data,variantdata)
        validationdata = {
            "chrom":["1","2"],
            "pos":[1,10],
            "rsid":[12345,54321],
            "ref":["A","A"],
            "alt":["C","C,G"]
        }
        validationdata = [
            RsidVar(Variant("1",1,"A","C"),12345),
            RsidVar(Variant("2",10,"A","C,G"),54321)
        ]
        #validationdata = pd.DataFrame(validationdata)
        self.assertEqual(out,validationdata)

class TestLocalDB(unittest.TestCase):
    """Test LocalDB
    """

    def test_init(self):
        col_data = ["1","2","3","4","5","6"]
        cols = ["SNP_ID_CURRENT","CHR_ID","CHR_POS","P-VALUE","PVALUE_MLOG","MAPPED_TRAIT","MAPPED_TRAIT_URI","LINK","STUDY"]
        transforms = [int,str,int,float,float,str,str,str,str]
        data = {a:list(map(b,col_data)) for (a,b) in itertools.zip_longest(cols,transforms) }
        data=pd.DataFrame(data)
        buf = StringIO()
        data.to_csv(buf,index=False,sep="\t")
        buf.seek(0)
        alldb = MockAlleleDB()
        db = gwcatalog_api.LocalDB(buf,10.0,100,alldb)

    def test_associations_for_regions(self):
        data = {
            "SNP_ID_CURRENT": [12345,54321,90001],
            "CHR_ID":["1","2","3"],
            "CHR_POS":[1,10,100],
            "P-VALUE":[0.1,0.1,0.1],
            "PVALUE_MLOG":[2.,2.,2.],
            "MAPPED_TRAIT":["A1","B1","C1"],
            "MAPPED_TRAIT_URI":["1","2","3"],
            "LINK":["link","link","link"],
            "STUDY":["study1","study2","study3"]
        }
        data=pd.DataFrame(data)
        buf = StringIO()
        data.to_csv(buf,index=False,sep="\t")
        buf.seek(0)
        alldb = MockAlleleDB()
        db = gwcatalog_api.LocalDB(buf,1.0,100,alldb)
        regions = [
            Region("1",1,2),
            Region("2",10,12),
            Region("3",100,102)
        ]
        assocs = db.associations_for_regions(regions)
        validationdata = [{
            "rsid": 12345,
            "chrom":"1",
            "pos":1,
            "ref":"A",
            "alt":"C",
            "pval":0.1,
            "pval_mlog":2.,
            "trait_name":"A1",
            "trait":"1",
            "study_link":"link",
            "study":"study1"
            }
        ]
        self.assertEqual(assocs,validationdata)


class testGWASDB(unittest.TestCase):
    def test_init(self):
        mockdb = MockAlleleDB()
        db = gwcatalog_api.GwasApi(0.1,100,1,mockdb)

    def test_associations_for_regions(self):
        mockdb = MockAlleleDB()
        db = gwcatalog_api.GwasApi(1.0,100,1,mockdb)
        ranges = [
            Region("1",1,2),
            Region("2",10,12),
            Region("3",100,102)
        ]
        
        request_val = mock.Mock() 
        data = data = {
            "SNP_ID_CURRENT": [12345,54321,90001],
            "CHR_ID":["1","2","3"],
            "CHR_POS":[1,10,100],
            "P-VALUE":[0.1,0.1,0.1],
            "PVALUE_MLOG":[2.,2.,2.],
            "MAPPED_TRAIT":["A1","B1","C1"],
            "MAPPED_TRAIT_URI":["A1","B2","C3"],
            "LINK":["link","link","link"],
            "STUDY":["study1","study2","study3"]
        }
        data=pd.DataFrame(data)
        buf = StringIO()
        data.to_csv(buf,index=False,sep="\t")
        buf.seek(0)
        request_val.text = buf.getvalue()
        
        with mock.patch("Scripts.data_access.gwcatalog_api.try_request",return_value=request_val) as mock_request:
            assocs = db.associations_for_regions(ranges)
        validationdata = [{
            "rsid": 12345,
            "chrom":"1",
            "pos":1,
            "ref":"A",
            "alt":"C",
            "pval":0.1,
            "pval_mlog":2.,
            "trait_name":"A1",
            "trait":"A1",
            "study_link":"link",
            "study":"study1"
            }
        ]
        self.assertEqual(assocs,validationdata)

if __name__=="__main__":
    unittest.main()