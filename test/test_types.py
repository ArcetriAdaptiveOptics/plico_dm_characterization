'''
Authors
  - C. Selmi: written in 2023
'''
import os
import numpy as np
import time
import shutil
import unittest
import unittest.mock as mock
from test.test_helper import testDataRootDir

class TestTypes(unittest.TestCase):

    @mock.patch('plico_dm_characterization.type.modalBase.ModalBase._storageFolder', autospec=True)
    def testModalBaseClass(self, mock_path):
        from plico_dm_characterization.type.modalBase import ModalBase
        mb = ModalBase()
        mock_path.return_value = os.path.join(testDataRootDir(), 'ModalBase')
        mb.getHadamardMatrix(11)
        modal_base = mb.getZonalMatrix(11)
        mb.getModalBase()
        mb.getTag()
        mb.getFitsFileName()
        mb.saveAsFits('testModalBase', modal_base)
        mb.saveAsH5('testModalBase', modal_base)
        mb.loadFromFits('testModalBase')
        mb.loadFromH5('testModalBase')

    @mock.patch('plico_dm_characterization.type.modalAmplitude.ModalAmplitude._storageFolder', autospec=True)
    def testModalAmplitudeClass(self, mock_path):
        from plico_dm_characterization.type.modalAmplitude import ModalAmplitude
        ma = ModalAmplitude()
        mock_path.return_value = os.path.join(testDataRootDir(), 'ModalAmplitude')
        ma.getModalAmplitude()
        ma.getTag()
        ma.getFitsFileName()
        modal_amplitude = np.zeros(11)*0.1
        ma.saveAsFits('testModalAmplitude', modal_amplitude)
        ma.saveAsH5('testModalAmplitude', modal_amplitude)
        ma.loadFromFits('testModalAmplitude')
        ma.loadFromH5('testModalAmplitude')

    @mock.patch('plico_dm_characterization.type.commandHistory.CmdHistory._storageFolder', autospec=True)
    def testCommandHistoryClass(self, mock_path):
        from plico_dm_characterization.type.commandHistory import CmdHistory
        cmdH = CmdHistory(11)
        mock_path.return_value = os.path.join(testDataRootDir(), 'CommandHistory')
        
        amp = np.zeros(11)*0.1
        modal_base = np.zeros((11, 11))
        np.fill_diagonal(modal_base, 1)
        n_rep = 1
        template = np.array([1, -1, 1])
        acts_vec = np.arange(11)

        cmd1, tt1 = cmdH.tidyCommandHistoryMaker(acts_vec, amp, modal_base, n_rep, template)
        time.sleep(2)
        cmd2, tt2 = cmdH.shuffleCommandHistoryMaker(acts_vec, amp, modal_base, n_rep, template=None)
        
        cmdH1 = CmdHistory.load(tt1, 0)
        cmd = cmdH1.getCommandHistory()
        self.assertTrue((cmd==cmd1).all())
        ind_list = cmdH1.getIndexingList()
        
        import platform
        if platform.system() == 'Windows':
            pass
        else:
            shutil.rmtree(os.path.join(testDataRootDir(), 'CommandHistory', tt1))
            shutil.rmtree(os.path.join(testDataRootDir(), 'CommandHistory', tt2))
