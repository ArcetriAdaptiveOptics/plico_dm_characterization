# plico_dm_characterization

This is part a component of the [plico][plico] framework and it design to calibrate and characterize 
Deformable Mirror (DMs).
This application uses [plico_interferometers][plico_interferometers] and [plico_dm][plico_dm].


## Prerequisites
Before starting to use this library make sure you are using the interferometer and deformable 
mirror of interest through the applications [plico_interferometers][plico_interferometers] and 
[plico_dm][plico_dm].


## Description (and use) of the calibration algorithms
__Acquisition and analysis of Influence Functions__

The influence functions are measured by moving each actuator (zonal commands) by the same 
quantity in the positive and negative direction relative to 
the mirror's initial position. Once the interferometric images corresponding to the positive 
and negative displacement of the actuator have been acquired, they are subtracted and normalized 
with respect to the chosen amplitude command. The resulting image shows the phase value for each 
interferometer pixel of that particular actuator influence function. This operation is then repeated 
for all actuators from which the mirror is composed.

The commands to be used to perform measurement and data analysis are:

- Open a terminal

- Define the deformable mirror:
	```
	import plico_dm
	dm =  plico_dm.deformableMirror(hostServer, portServer)
	```
- Define the interferometer:
	```
	import plico_interferometer
	interf = plico_interferometer.interferometer(hostServer, portServer)

	```
- Define the class for Influence Function acquisition:
	```
	from plico_dm_characterization.influenceFunctionsMaker import IFMaker
	iff = IFMaker(interf, dm)
	```
- Define the modal base and amplitude to use:
	```
	modalBaseTag = 'YourModalBase.fits'; ampTag = 'YourAmplitude.fits'
	```
	To create your fits files containing the modal base to be used use the following instructions
	```
	from plico_dm_characterization.type.modalBase import ModalBase
	mb = ModalBase()
	your_modal_base = mb.getHadamardMatrix(nActs) for global IFs or your_modal_base = mb.getZonalMatrix(nActs) for zonal IFs
	mb.saveAsFits(‘YourModalBase’,  your_modal_base)
	```
	To create your fits files containing the vector of amplitudes with which to make measurements use the following instructions
	```
	from plico_dm_characterization.type.modalAmplitude import ModalAmplitude
	ma = ModalAmplitude()
	your_amplitude = np.array(nActs); the vector containing the amplitude of the modes (or zonal command)  in mirror units 
	ma.saveAsFits(‘YourAmplitude’,  your_amplitude)
	```
- Start the acquisition:
	```
	tn = iff.acquisitionAndAnalysis(modalBaseTag, ampTag, shuffle=False, template=None, n_rep=1)
	```

Acquisition is done by creating a matrix with the history of the commands to be applied to the mirror 
and the subsequent serial application of these: the history matrix is realized depending on the 
shuffle, template and n_rep parameters passed in the argument. If shuffle is equal to None, the modal 
base will be applied following the order of the elements in it; otherwise, shuffle=True, the elements 
will be applied in random order keeping track of them. If not specified, template= None, the algorithm 
uses the sequence 1, -1, 1 for application of each mode: by appropriately combining the three 
measurements, linear drift is eliminated. Any other definition of the template vector can be applied 
by passing it as an argument to the measurement and analysis function (e.g. the sequence 1, -1 
corresponding to the push/pull). The n_rep parameter allows the template vector to be repeated several 
times.
The process of acquiring and analyzing influence functions returns the tracking number string 
containing the data: all acquired interferograms and the cube containing the analysis of these with 
all the information characterizing the measurement (YourModalBase, YourAmplitude, shuffle, template 
and n_rep) are saved.

__From Wavefront to Deformable Mirror command__

The influence functions obtained above constitute the Interaction Matrix: calculation of the pseudo 
inverse of this matrix returns the Reconstructor. The reconstructor allows the passage from wavefront 
to command in mirror units (UC). 
Having then acquired an interferogram, it is possible to pass it as an argument to the function that 
performs the conversion: after appropriately combining the mask with which the interferogram was 
acquired and the mask resulting from the cube of measurements taken, the algorithm returns the 
command in units of the mirror by executing the product between the reconstructor and the wavefront 
to be converted.
The commands to be used to perform the conversion are
- Define the tracking number of the influence function measurements to be used
- Use the commands
	```
	from plico_dm_characterization.convertWFToDmCommand import Converter
	cc = Converter(tt, zonal=True)
	```
If zonal is not True the command is calculated by multiplying the standard command to Hadamard matrix.
	```
	wf = interf.wavefront()
	cmd = cc.fromWfToDmCommand(wf)
	```


__Get command for Zernike modes on DM__

With the class Converter it is also possible to calculate the commands to be given to the mirror to obtain the first n_modes Zernike modes: this computation is done by multiplying the theoretical Zernike surfaces, constructed on the intersection mask of the measurement cube, with the reconstructor.
The commands to be used for this is
```
zernike_command_matrix = cc.getCommandsForZernikeModeOnDM(n_modes)
```
where zernike_command_matrix has the shape [nActs, n_modes] and the commands are expressed in meters.
It is possible to save the matrix using the command
```
cc.saveZernikeMat(zernike_command_matrix, file_path)
```


## Data storage
At the end of the calibration of a deformable mirror all essential data for its reuse in a different 
system setup,  i.e. DM not in front of interferometer, are saved in gDrive at the following link. 
This folder is part of the PLICO shared drive where all laboratory instrumentation data linked to 
the framework are stored.
Within the folder for the mirror you want to use you can find
- a read_me file in which information on the influence functions produced for calibration is collected
- the interaction matrix: stored inside the folder IFFunctions at the relative tracking number
- all the information about the calibration: modal base, amplitude of modes, command history for the DM…
- the flat command:  stored inside the folder Flattening at the relative tracking number
- other information about the interferometer as his configuration file or the detector mask used
- the Zerbine command matrix for the mirror: stored inside the folder Zernike



[plico]: https://github.com/ArcetriAdaptiveOptics/plico
[plico_interferometers]: https://github.com/ArcetriAdaptiveOptics/plico_interferometer
[plico_dm]: https://github.com/ArcetriAdaptiveOptics/plico_dm
