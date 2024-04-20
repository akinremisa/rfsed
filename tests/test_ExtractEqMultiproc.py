import unittest
from rfsed import ExtractEqMultiproc
from rfsed.ExtractEqMultiproc import Add_sta_time_2cat, Get_Eqtimes_multiproc, ExtractEq_Multiproc
from rfsed.util import raw_data_example, save_tests, catalog_example
from glob import glob
from obspy import read_events

class ExtractEqMultiproc(unittest.TestCase):
    def test_Add_sta_time_2cat(self):
        catalog = catalog_example()
        Request_window=[-50, 150]
        stalat=6.183432
        stalon=53.48334
        Add_sta_time_2cat(catalogfile=catalog, stalat=stalat, stalon=stalon, Request_window=Request_window)
    def test_Get_Eqtimes_multiproc(self):
        catalog = catalog_example()
        Request_window=[-50, 150]
        stalat=6.183432
        stalon=53.48334
        catalog_time = Add_sta_time_2cat(catalogfile=catalog, stalat=stalat, stalon=stalon, Request_window=Request_window)
        nproc=25
        Get_Eqtimes_multiproc(catalog=catalog_time, nproc=nproc)

    def test_ExtractEq_Multiproc(self):
        datapath=raw_data_example()
        datafiles=datafiles = sorted(glob("%s*.dat"%(datapath)))
        staname='NE301'
        catalog = catalog_example()
        Request_window=[-50, 150]
        stalat=6.183432
        stalon=53.48334
        nproc=25
        savepath=save_tests()
        filename = savepath + 'eq_data_multiproc.h5'
        catalog_time = Add_sta_time_2cat(catalogfile=catalog, stalat=stalat, stalon=stalon, Request_window=Request_window)
        timewindow=Get_Eqtimes_multiproc(catalog=catalog_time, nproc=nproc)
        ExtractEq_Multiproc(datafiles=datafiles, nproc=nproc, filename=filename,  timewindow=timewindow)
        
        
def suite():
    return unittest.makeSuite(ExtractEqMultiproc, 'test')


if __name__ == '__main__':
    unittest.main(defaultTest='suite')