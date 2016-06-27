# Photoluminescence-Work-up

This repository is for my Matlab functions/scripts that handle PL data work-up. These are intedned to work for all types of PL data from any spectrometer, so long as they are imported as column vectors or a matrix of multiple column vectors (ie. 2500x5).

Abs_avg serves to determine the average absorbance of a sample over the range that it was excited for PL. The .m file should explain the inputs and outputs in its comments.

QY_Calc calcualtes the quantum yield of a sample provided another reference sample. Again, the .m file is commented to explain its use.

PL_Master incorporates the functionality of QY_Calc (with some improved plotting functionality) as well as converting PL data from wavelength (nm) to energy (ev), normalizing, plotting, and fitting of spectra to 1st order Gaussian functions at the user's request.

If you use these functions and run into any issues or desire more functionality, please do not hesistate to contact me with bugs/concerns/requests at MarkCBabin@gmail.com.
