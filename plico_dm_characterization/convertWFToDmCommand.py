'''
Authors
  - C. Selmi:  written in 2021
'''
import numpy as np
import os
from plico_dm_characterization.influenceFunctionsMaker import IFMaker
from scipy.linalg import hadamard
from plico_dm_characterization.ground import zernike


class Converter():
    '''
    HOW TO USE IT::

        from plico_dm_characterization.convertWFToDmCommand import Converter
        tt = '...'
        cc = Converter(tn)
        cmd = cc.fromWfToDmCommand(wf)
    '''

    def __init__(self, tt_an, zonal=True):
        an = IFMaker.loadAnalyzerFromIFMaker(tt_an)
        self._cube = an.getCube()
        self.zonal = zonal
        #analisi
        self._analysisMask = None
        self._intMat = None
        self._rec = None


    def fromWfToDmCommand(self, wf):
        '''
        Parameters
        ----------
        wf: numpy masked array
            image to convert in a command for deformable mirror

        Return
        ------
        cmd: numpy array [nActs]
            command for deformable mirror
        '''
        #manca l'utilizzo delle coordinate
        new_mask = np.ma.mask_or(wf.mask, self.getMasterMask())
        self.setAnalysisMask(new_mask)
        wf_masked = np.ma.masked_array(wf.data, mask=new_mask)
        rec = self.getReconstructor()
        command = np.dot(rec, wf_masked.compressed())
        if self.zonal is not True:
            mat = hadamard(128)
            hadaMat = mat[0:self._cube.shape[2], 0:self._cube.shape[2]]
            command = np.matmul(hadaMat, command)
        return command

    def getCommandsForZernikeModeOnDM(self):
        '''
        Function for the calculation of the zernike command (in Unit Mirror).
        The first 10 Zernike are calculate as the product of the theoretical
        modes on cube master mask and the reconstructor.
        '''
        zernike_cube = self._createZernikeOnDM()
        self.setAnalysisMaskFromMasterMask()
        rec = self.getReconstructor()

        zernike_command_matrix_list = []
        for i in range(zernike_cube.shape[2]):
            comm = np.dot(rec, zernike_cube[:, :, i].compressed())
            zernike_command_matrix_list.append(comm)
        zernike_command_matrix = np.array(zernike_command_matrix_list)
        return zernike_command_matrix.T

    def _createZernikeOnDM(self):
        ima = self._cube[:, :, 0]
        mask = self.getMasterMask()
        image = np.ma.masked_array(ima.data, mask=mask)
        coef, mat = zernike.zernikeFit(image, np.arange(10) + 1)
        surf_list = []
        for i in range(coef.size):
            surf = zernike.zernikeSurface(image, 1., mat[:, i])
            surf_list.append(surf)
        zernike_cube = np.ma.dstack(surf_list)
        return zernike_cube

    def saveZernikeMat(self, zernike_command_matrix, dove):
        from astropy.io import fits
        fits.writeto(os.path.join(dove, 'ZernikeMat.fits'), zernike_command_matrix)
        return

    def getMasterMask(self):
        '''
        Returns
        -------
        master_mask: [pixels, pixels]
                    product of the masks of the cube
        '''
        aa = np.sum(self._cube.mask.astype(int), axis=2)
        master_mask = np.zeros(aa.shape, dtype=np.bool)
        master_mask[np.where(aa > 0)] = True
        return master_mask

    def setAnalysisMask(self, analysis_mask):
        ''' Set the analysis mask chosen

        Parameters
        ----------
        analysis_mask: numpy array [pixels, pixels]
        '''
        self._intMat = None
        self._rec = None
        self._analysisMask = analysis_mask

    def setAnalysisMaskFromMasterMask(self):
        ''' Set the analysis mask using the master mask of analysis cube
        '''
        self.setAnalysisMask(self.getMasterMask())

    # def setDetectorMask(self, mask_from_ima):
    #     ''' Set the detector mask chosen
    #
    #     Parameters
    #     ----------
    #     detector_mask: numpy array [pixels, pixels]
    #     '''
    #     self._analysisMask = None
    #     self._intMat = None
    #     self._rec = None
    #     self._analysisMask = mask_from_ima

    def getAnalysisMask(self):
        '''
        Returns
        -------
        analysis_mask: numpy array [pixels, pixels]
        '''
        return self._analysisMask

    def _getMaskedInfluenceFunction(self, idx_influence_function):
        return np.ma.array(self._cube[:, :, idx_influence_function],
                           mask=self.getAnalysisMask())

    def _createInteractionMatrix(self):
        if self._analysisMask is None:
            self.setAnalysisMaskFromMasterMask()
        n_acts_in_cube = self._cube.shape[2]
        n_interferometer_pixels_in_mask = \
                    self._getMaskedInfluenceFunction(0).compressed().shape[0]
        self._intMat = np.zeros((n_interferometer_pixels_in_mask,
                                 n_acts_in_cube))
        for i in range(n_acts_in_cube):
            self._intMat[:, i] = \
                            self._getMaskedInfluenceFunction(i).compressed()

    def _createSurfaceReconstructor(self, rCond=1e-15):
        self._rec = self._createRecWithPseudoInverse(rCond)

    def _createRecWithPseudoInverse(self, rCond):
        return np.linalg.pinv(self.getInteractionMatrix(), rcond=rCond)

    def getInteractionMatrix(self):
        '''
        Returns
        -------
        intMat: numpy array
                interaction matrix from cube
        '''
        if self._intMat is None:
            self._createInteractionMatrix()
        return self._intMat

    def getReconstructor(self):
        '''
        Returns
        -------
        rec = numpy array
            reconstructor calculated as pseudo inverse of the interaction matrix
        '''
        if self._rec is None:
            self._createSurfaceReconstructor()
        return self._rec
