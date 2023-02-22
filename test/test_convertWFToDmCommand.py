'''
Authors
  - C. Selmi: written in February 2023
'''
import os
import numpy as np
from astropy.io import fits
import shutil
import unittest
import unittest.mock as mock
from test.test_helper import testDataRootDir
from plico_dm_characterization.convertWFToDmCommand import Converter

class TestConverterWF(unittest.TestCase):

    @mock.patch('plico_dm_characterization.influenceFunctionsMaker.IFMaker._storageFolder', autospec=True)
    def setUp(self, mock_path_iff):
        mock_path_iff.return_value = os.path.join(testDataRootDir(), 'IFFunctions')
        self.cc = Converter('20230222_102118')

    def tearDown(self):
        del(self.cc)

    
    def testConverter(self):
        path = os.path.join(testDataRootDir(), 'imaAlpao88.fits')
        hduList = fits.open(path)
        wf = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))

        command = self.cc.fromWfToDmCommand(wf)
        self.cc.getInteractionMatrix()
        mask = self.cc.getMasterMask()
        self.cc.setAnalysisMask(mask)
        self.cc.setAnalysisMaskFromMasterMask()
        self.cc.getAnalysisMask()
        self.cc.getReconstructor()
        #zernike_matrix = self.cc.getCommandsForZernikeModeOnDM(3)
