'''
Authors
  - C. Selmi:  written in 2019 and 2022
'''
import numpy as np
import os
#import glob
from astropy.io import fits
from matplotlib import pyplot as plt
from plico_dm_characterization.configuration import config


def linearity(tn):
    '''
    Parameters
    ----------
    tn: string
        tracking number of linearity measurements

    Returns
    -------
    linearityMatrix: numpy array [n_amp: n_acts]
        
    '''
    fold_for_meas = os.path.join(config.LINEARITY_ROOT_FOLD, tn)
    acts_list = os.listdir(fold_for_meas)
    acts_list.sort()
    #del(tt_list[0])
    del(acts_list[-2:])
    hduList = fits.open(os.path.join(fold_for_meas, 'amplitude.fits'))
    amp = hduList[0].data
    hduList = fits.open(os.path.join(fold_for_meas, 'actuators.fits'))
    acts = hduList[0].data

    cube = None
    for act in acts_list:
        hduList = fits.open(os.path.join(fold_for_meas, act))
        cube_act = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
        if cube is None:
            cube = cube_act
        else:
            cube = np.ma.dstack((cube, cube_act))

    linearityMatrix = _createLinearityMatrix(amp, acts, cube, 10)
    aVector, bVector, fit, diff, rms = _calcFitDiffAndRMSAct(amp, acts, linearityMatrix)
    _plotResults(linearityMatrix, fit, diff, acts, amp)
    return linearityMatrix, rms

def _createLinearityMatrix(amp, acts, cube, pixelCut):
    linearityMatrix = np.zeros((amp.size, acts.size))
    for i in range(amp.size):
        for j in range(acts.size):
            k = i*acts.size+j
            ima = cube[:,:,k]
            y, x = np.where(ima==np.amax(ima))
            imgCut = ima[y[0]-pixelCut:y[0]+pixelCut, x[0]-pixelCut:x[0]+pixelCut]
            m = np.mean(imgCut)
            value = m * amp[i]
            linearityMatrix[i, j] = value
    return linearityMatrix

def _calcFitDiffAndRMSAct(amp, acts, linearityMatrix):
    diff = np.zeros((amp.size, acts.size))
    aVector = np.zeros(acts.size)
    bVector = np.zeros(acts.size)
    fit = np.zeros((amp.size, acts.size))
    rms = np.zeros(acts.size)
    rmsLine = np.zeros(acts.size)

    for i in range(amp.size):
        for j in range(linearityMatrix.shape[1]):    
            a, b = np.polyfit(amp, linearityMatrix[:, j], 1)
            bVector[j] = b
            aVector[j] = a
            fit[:, j] = a*amp+b
            diff[:, j] = linearityMatrix[:, j] - fit[:, j] 
            rms[j] = diff[:, j].std()
            #rms[j]= np.sqrt(np.mean(np.square(diff[:,j])))
            rmsLine[j] = np.sqrt(np.mean(np.square(fit[:, j])))
    return aVector, bVector, fit, diff, rms

def _plotResults(linearityMatrix, fit, diff, acts, amp):
    n_figure = linearityMatrix.shape[1]
    plt.figure(figsize=(10, 16))
    
    for i in range(n_figure):
        plt.subplot(n_figure, 1, i+1)
        plt.plot(amp, linearityMatrix[:, i]*1e09, '-o', label='misure')
        plt.plot(amp, fit[:, i]*1e09, '-o', label='fit')
        plt.plot(amp, diff[:, i]*1e09, '-o', label='diff')
        plt.legend(loc='upper left')
        plt.xlabel('amp [UC]')
        plt.ylabel('displacement [nm]')
        plt.title('act #%d' %acts[i])

        
    
