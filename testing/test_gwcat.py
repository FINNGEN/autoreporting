import unittest
import unittest.mock as mock
import sys,os,json, requests
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
        with mock.patch("Scripts.gwcatalog_api.requests.get",return_value=response):
            retval=gwcatalog_api.try_request(url, params)
        self.assertTrue(type(retval)!= type(None))
        #404
        response=mock.Mock()
        response.status_code=404
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.get",return_value=response):
            retval=gwcatalog_api.try_request(url, params)
        self.assertTrue(type(retval)== type(None))
        #400
        response=mock.Mock()
        response.status_code=400
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.get",return_value=response):
            retval=gwcatalog_api.try_request(url, params)
        self.assertTrue(type(retval)== type(None))
        #500
        response=mock.Mock()
        response.status_code=500
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.get",return_value=response) as mock_get:
            with mock.patch("Scripts.gwcatalog_api.print"):
                retval=gwcatalog_api.try_request(url, params)
        self.assertTrue(type(retval)== type(None))
        self.assertEqual(5,mock_get.call_count)
        #exception during requests.get
        with mock.patch("Scripts.gwcatalog_api.requests.get",side_effect=TimeoutError("Timeout")) as mock_get:
            with mock.patch("Scripts.gwcatalog_api.print") as mock_print:
                retval=gwcatalog_api.try_request(url, params)
        self.assertEqual(type(retval),type(None))
        mock_print.assert_called_with("Request caused an exception:{}".format(TimeoutError("Timeout")))
    
    @mock.patch("Scripts.gwcatalog_api.time.sleep")
    def test_post(self,mock_time):
        """
        Test try_request_post
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
        with mock.patch("Scripts.gwcatalog_api.requests.post",return_value=response):
            retval=gwcatalog_api.try_request_post(url,headers, data)
        self.assertTrue(type(retval)!= type(None))
        #404
        response=mock.Mock()
        response.status_code=404
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.post",return_value=response):
            retval=gwcatalog_api.try_request_post(url,headers, data)
        self.assertTrue(type(retval)== type(None))
        #400
        response=mock.Mock()
        response.status_code=400
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.post",return_value=response):
            retval=gwcatalog_api.try_request_post(url,headers, data)
        self.assertTrue(type(retval)== type(None))
        #500
        response=mock.Mock()
        response.status_code=500
        url="url"
        params="params"
        with mock.patch("Scripts.gwcatalog_api.requests.post",return_value=response) as mock_get:
            with mock.patch("Scripts.gwcatalog_api.print"):
                retval=gwcatalog_api.try_request_post(url,headers, data)
        self.assertTrue(type(retval)== type(None))
        self.assertEqual(5,mock_get.call_count)
        #exception during requests.get
        with mock.patch("Scripts.gwcatalog_api.requests.post",side_effect=TimeoutError("Timeout")) as mock_get:
            with mock.patch("Scripts.gwcatalog_api.print") as mock_print:
                retval=gwcatalog_api.try_request_post(url,headers, data)
        self.assertEqual(type(retval),type(None))
        mock_print.assert_called_with("Request caused an exception:{}".format(TimeoutError("Timeout")))
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
        mock_try.assert_called_with(tryurl)
        #does not exist
        with mock.patch("Scripts.gwcatalog_api.try_request",return_value = None) as mock_try:
            with mock.patch("Scripts.gwcatalog_api.print"):
                return_value=gwcatalog_api.get_trait_name(trait)
        self.assertEqual("NA",return_value)

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

if __name__=="__main__":
    os.chdir("./testing")
    unittest.main()