import unittest
from rfsed import ExtractEq
from rfsed.util import raw_data_example, save_tests, catalog_example

class TestExtractEq(unittest.TestCase):
    def test_ExtractEq(self):
        datapath=raw_data_example()
        staname='NE301'
        stalat=6.183432
        stalon=53.48334
        savepath=save_tests()
        Request_window=[-50, 150]
        filename = savepath + 'eq_data.h5'
        catalog = catalog_example()
        ExtractEq.ExtractEq(datapath=datapath, catalog=catalog, filename=filename, stalat=stalat, stalon=stalon, Request_window=Request_window)
        
        
def suite():
    return unittest.makeSuite(TestExtractEq, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')