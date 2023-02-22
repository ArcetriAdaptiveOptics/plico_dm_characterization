'''
Authors
  - C. Selmi: written in February 2023
'''
import os
import numpy as np
import unittest
from astropy.io import fits
from test.test_helper import testDataRootDir
from plico_dm_characterization.ground import geo

class TestGeo(unittest.TestCase):
    
    def test_geo(self):
        path = os.path.join(testDataRootDir(), 'imaAlpao88.fits')
        hduList = fits.open(path)
        image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
        
        ima = geo.draw_mask(image, 400, 400, 30)
        geo.qpupil(image.mask)
        geo.rotate(image, 45)
