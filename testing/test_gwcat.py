import unittest
import unittest.mock as mock
import sys,os,json, requests
import pandas as pd
sys.path.append("../")
sys.path.append("./")
sys.path.insert(0, './Scripts')
from Scripts import gwcatalog_api

class TestGwcat(unittest.TestCase):
    @mock.patch("Scripts.gwcatalog_api.time.sleep")
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
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response):
            retval=gwcatalog_api.try_request("GET",url, params=params)
        self.assertTrue(type(retval)!= type(None))
        #404
        response=mock.Mock()
        response.status_code=404
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("GET",url, params=params)
        #400
        response=mock.Mock()
        response.status_code=400
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("GET",url, params=params)
        #500
        response=mock.Mock()
        response.status_code=500
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response) as mock_get:
            with mock.patch("Scripts.gwcatalog_api.print"):
                with self.assertRaises(gwcatalog_api.ResponseFailure):
                    retval=gwcatalog_api.try_request("GET",url, params=params)
        self.assertEqual(5,mock_get.call_count)
        #exception during requests.get
        with mock.patch("Scripts.gwcatalog_api.requests.request",side_effect=TimeoutError("Timeout")) as mock_get:
            with mock.patch("Scripts.gwcatalog_api.print") as mock_print:
                with self.assertRaises(gwcatalog_api.ResponseFailure):
                    retval=gwcatalog_api.try_request("GET",url, params=params)
        mock_print.assert_called_with("Request caused an exception:{}".format(TimeoutError("Timeout")))

    @mock.patch("Scripts.gwcatalog_api.time.sleep")
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
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response):
            retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        self.assertTrue(type(retval)!= type(None))
        #404
        response=mock.Mock()
        response.status_code=404
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        #400
        response=mock.Mock()
        response.status_code=400
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response):
            with self.assertRaises(gwcatalog_api.ResourceNotFound):
                retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        #500
        response=mock.Mock()
        response.status_code=500
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.request",return_value=response) as mock_get:
            with self.assertRaises(gwcatalog_api.ResponseFailure):
                with mock.patch("Scripts.gwcatalog_api.print"):
                    retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        #exception during requests.post
        with mock.patch("Scripts.gwcatalog_api.requests.request",side_effect=TimeoutError("Timeout")) as mock_get:
            with self.assertRaises(gwcatalog_api.ResponseFailure):
                with mock.patch("Scripts.gwcatalog_api.print") as mock_print:
                    retval=gwcatalog_api.try_request("POST",url,headers=headers, data=data)
        pass 

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
        with mock.patch("Scripts.gwcatalog_api.try_request",return_value = response_mock) as mock_try:
            return_value=gwcatalog_api.get_trait_name(trait)
        self.assertEqual(return_value,trait_name)
        mock_try.assert_called_with("GET",url=tryurl)
        #does not exist, raises exception gwcatalog_api.ResourceNotFound
        with mock.patch("Scripts.gwcatalog_api.try_request",side_effect = gwcatalog_api.ResourceNotFound):
            with mock.patch("Scripts.gwcatalog_api.print"):
                return_value=gwcatalog_api.get_trait_name(trait)
        self.assertEqual("NA",return_value)
        #TODO: handle unhandled exception

    def test_parse_float(self):
        """
        Test parse_float
        Test cases:
            number: 5e-8    => '5e-8'
            number: 1.2e-8  => '12e-9'
            number: 0       => '0e0'
            number: 12 => '1'
        """
        num=5e-8
        validate='5e-8'
        retval=gwcatalog_api.parse_float(num)
        self.assertEqual(retval,validate)
        num=1.2e-8
        validate='12e-9'
        retval=gwcatalog_api.parse_float(num)
        self.assertEqual(retval,validate)
        num=0.0
        validate='0'
        retval=gwcatalog_api.parse_float(num)
        self.assertEqual(retval,validate)
        num=12.0
        validate='1'
        retval=gwcatalog_api.parse_float(num)
        self.assertEqual(retval,validate)

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
        with mock.patch("Scripts.gwcatalog_api.print") as mock_print:
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


if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()