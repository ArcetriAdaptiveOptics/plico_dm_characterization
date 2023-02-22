'''
Authors
  - C. Selmi: written in February 2023
'''
import os
import numpy as np
import unittest
from astropy.io import fits
from test.test_helper import testDataRootDir
from plico_dm_characterization.ground import zernike

class TestZernike(unittest.TestCase):
    
    def test_zernike(self):
        path = os.path.join(testDataRootDir(), 'imaAlpao88.fits')
        hduList = fits.open(path)
        image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))

        coeff, mat = zernike.zernikeFit(image, np.arange(10)+1)
        zernike.zernikeSurface(image, coeff, mat)
