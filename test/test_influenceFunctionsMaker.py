'''
Authors
  - C. Selmi: written in February 2023
'''
import os
import numpy as np
import time
from astropy.io import fits
import shutil
import unittest
import unittest.mock as mock
from test.test_helper import testDataRootDir
from plico_dm_characterization.influenceFunctionsMaker import IFMaker

class TestInfluenceFunctionsMaker(unittest.TestCase):
    
    def setUp(self):
        self.dm, self.interf = self._createDMAndInterf()
        self.iff = IFMaker(self.interf, self.dm)

    def tearDown(self):
        del(self.dm, self.interf)
    
    def _createDMAndInterf(self):
        class DM():
            def __init__(self, nActs):
                self.nActs = nActs
                self.cmd = np.zeros(nActs)
            def get_number_of_actuators(self):
                return self.nActs
            def get_shape(self):
                return self.cmd
            def set_shape(self, cmd):
                self.cmd = cmd
                return self.cmd
        
        class Interf():
            def wavefront(self, n_images):
                path = os.path.join(testDataRootDir(), 'imaAlpao88.fits')
                hduList = fits.open(path)
                image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
                return image
        
        dm = DM(3)
        interf = Interf()
        return dm, interf
    
    @mock.patch('plico_dm_characterization.influenceFunctionsMaker.IFMaker._storageFolder', autospec=True)
    @mock.patch('plico_dm_characterization.type.modalBase.ModalBase._storageFolder', autospec=True)
    @mock.patch('plico_dm_characterization.type.modalAmplitude.ModalAmplitude._storageFolder', autospec=True)
    @mock.patch('plico_dm_characterization.type.commandHistory.CmdHistory._storageFolder', autospec=True)
    def testIFMaker(self, mock_path_cmd, mock_path_ma, mock_path_mb, mock_path):
        mock_path.return_value = os.path.join(testDataRootDir(), 'IFFunctions')
        cmd_matrix_tag = 'testModalBase'
        mock_path_mb.return_value = os.path.join(testDataRootDir(), 'ModalBase')
        amplitude_tag = 'testModalAmplitude'
        mock_path_ma.return_value = os.path.join(testDataRootDir(), 'ModalAmplitude')
        mock_path_cmd.return_value = os.path.join(testDataRootDir(), 'CommandHistory')
        tt = self.iff.acquisitionAndAnalysis(cmd_matrix_tag, amplitude_tag)
        time.sleep(3)
        tt2 = self.iff.acquisitionAndAnalysis(cmd_matrix_tag, amplitude_tag,
                                              shuffle=True, template=np.array([1, -1, 1]))

        cube = self.iff.getCube()
        iff = IFMaker.loadAnalyzerFromIFMaker('20230222_102118')
        iff._storageFolder()
        
        import platform
        if platform.system() == 'Windows':
            pass
        else:
            shutil.rmtree(os.path.join(testDataRootDir(), 'IFFunctions', tt))
            shutil.rmtree(os.path.join(testDataRootDir(), 'IFFunctions', tt2))
