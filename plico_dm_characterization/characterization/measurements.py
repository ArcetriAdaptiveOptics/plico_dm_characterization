'''
Authors
  - C. Selmi:  written in 2019 and 2022
'''
import os
import numpy as np
import time
from astropy.io import fits
import shutil
from plico_dm_characterization.ground import zernike
from plico_dm_characterization.ground.timestamp import Timestamp
from plico_dm_characterization.influenceFunctionsMaker import IFMaker
from plico_dm_characterization.ground.tracking_number_folder import TtFolder
from plico_dm_characterization.configuration import config
from plico_dm_characterization.type.modalBase import ModalBase
from plico_dm_characterization.type.modalAmplitude import ModalAmplitude
from plico_dm_characterization.convertWFToDmCommand import Converter

class MeasurementAcquisition():
    '''
        HOW TO USE IT::

        from plico_dm_characterization.characterization.measurements import MeasurementAcquisition
        ma = MeasurementAcquisition(dm, interf)
    '''

    def __init__(self, deformable_mirror, interferometer):
        self.dm = deformable_mirror
        self.interf = interferometer

    def flattening(self, tn_iff, zernike_to_remove=None):
        '''
        Parameters
        ----------
        zernike_to_remove: numpy array
            zernike to be removed from image before alignment
            example: to remove piston, tip and tilt use zernike = np.arange(3) + 1
        tn_iff: string
            tracking number of if function to use for the reconstructor calculation
        '''
        converter = Converter(tn_iff)

        fold_for_meas = config.FLAT_ROOT_FOLD
        dove, tt = TtFolder(fold_for_meas).createFolderToStoreMeasurements()
        ima = self.interf.wavefront()
        if zernike_to_remove is not None:
            coef, mat = zernike.zernikeFit(ima, zernike_to_remove)
            surf = zernike.zernikeSurface(ima, coef, mat)
            ima_ttr = ima - surf
            command = converter.fromWfToDmCommand(ima_ttr)
        else:
            command = converter.fromWfToDmCommand(ima)

        fits.writeto(os.path.join(dove, 'imgstart.fits'), ima.data)
        fits.append(os.path.join(dove, 'imgstart.fits'), ima.mask.astype(int))
        fits.writeto(os.path.join(dove, 'flatDeltaCommand.fits'), command)

        pos = self.dm.get_shape()
        cmd_to_apply = pos - command
        maxComm = max(abs(cmd_to_apply))
        print('max command= %f' % maxComm)

        if maxComm>float(config.MAX_COMMAND_TO_APPLY):
            raise OSError('Actuator command too large')
        else:
            self.dm.set_shape(cmd_to_apply)

        self._commandToApply = cmd_to_apply
        
        wf = self.interf.wavefront()
        fits.writeto(os.path.join(dove, 'imgflat.fits'), wf.data)
        fits.append(os.path.join(dove, 'imgflat.fits'), wf.mask.astype(int))

        fits.writeto(os.path.join(dove, 'flatCommand.fits'), self._commandToApply)
        return tt

    def closeLoop(self, n_meas, delay=0):
        '''
        Parameters
        ----------
        n_meas: int
            number of times to repeat the flattening measurement

        Returns
        -------
        tn_list: list
            list of flattening tracking numbers
        '''
        tt_list = []
        for i in range(0, n_meas):
            tt = self.flattening()
            tt_list.append(tt)
            time.sleep(delay)
        return tt_list

    def opticalMonitoring(self, n_images, delay):
        '''
        Parameters
        ----------
        n_images: int
            number of images to acquire
        delay: int [s]
            waiting time (in seconds) between two image acquisitions

        Returns
        ------
        tt: string
            tracking number of measurements
        '''
        fold_for_meas = config.OPD_SERIES_ROOT_FOLD
        dove, tt = TtFolder(fold_for_meas).createFolderToStoreMeasurements()

        zer_list = []
        t0 = time.time()
        for i in range(n_images):
            ti = time.time()
            dt = ti - t0
            masked_ima = self.interf.wavefront(1)
            name = Timestamp.now() + '.fits'
            fits_file_name = os.path.join(dove, name)
            fits.writeto(fits_file_name, masked_ima.data)
            fits.append(fits_file_name, masked_ima.mask.astype(int))

            coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
            vect = np.append(dt, coef)
            zer_list.append(vect)

            fits_file_name = os.path.join(dove, 'zernike.fits')
            fits.writeto(fits_file_name, np.array(zer_list), overwrite=True)

            time.sleep(delay)

    def _actuator_linearity(self, actuator_to_test, vector_of_amplitude_for_meas):
        '''
        Parameters
        ----------
        actuators_to_test: int
            vector containing the number of actuator to test for linearity
        vector_of_amplitude_for_meas: numpy array
            vector containing the amplitude values to be given to the actuators
            to make the measurement

        Returns
        -------
        tt = string
            tracking number of measurements
        '''
        cmd0 = np.zeros(self.dm.get_number_of_actuators())

        ima_amp_list = []
        for j in range(vector_of_amplitude_for_meas.shape[0]):
            cmd = cmd0
            cmd[actuator_to_test] = vector_of_amplitude_for_meas[j]
            
            template = np.array([1,-1,1])
            ima_list = []
            for temp in template:
                self.dm.set_shape(cmd*temp)
                ima = self.interf.wavefront()
                ima_list.append(ima)
            final_ima = (ima_list[0] - 2.*ima_list[1] + ima_list[2]) / (2.*template.shape[0]-1)
            ima_amp_list.append(final_ima)
        return ima_amp_list
        
    def linearity(self, vector_of_actuators_to_test, vector_of_amplitude_for_meas): 
        fold_for_meas = config.LINEARITY_ROOT_FOLD
        dove, tt = TtFolder(fold_for_meas).createFolderToStoreMeasurements()
        fits.writeto(os.path.join(dove, 'actuators.fits'),
                     vector_of_actuators_to_test, overwrite=True)
        fits.writeto(os.path.join(dove, 'amplitude.fits'),
                     vector_of_amplitude_for_meas, overwrite=True)

        for i in range(vector_of_actuators_to_test.shape[0]):
            actuator_to_test = vector_of_actuators_to_test[i]
            ima_amp_list = self._actuator_linearity(actuator_to_test, vector_of_amplitude_for_meas)
            cube = np.ma.dstack(ima_amp_list)
            fits.writeto(os.path.join(dove, 'act_%03d.fits' %actuator_to_test),
                         cube.data)
            fits.append(os.path.join(dove, 'act_%03d.fits' %actuator_to_test),
                         cube.mask.astype(int))
        self.dm.set_shape(np.zeros(self.dm.get_number_of_actuators()))
        return tt

