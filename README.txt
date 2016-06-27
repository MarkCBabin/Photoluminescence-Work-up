Photoluminescence-Work-up

This repository is for my Matlab functions/scripts that handle PL data work-up. These are intended to work for all types of PL data from any spectrometer, so long as they are imported as column vectors or a matrix of multiple column vectors (ie. 2500x5).

PL_master_run is a script that will add the path of PL_master and Abs_avg to Matlab so that the files do not need to be copied into the same directory as the data being processed. It also implements a graphic UI to open files, rather than from the standard folder system implemented in Matlab. 
      It should be noted that this script must be edited so that the path being added to Matlab matches the location of PL_Master, at present a dummy directory is listed. 

Abs_avg serves to determine the average absorbance of a sample over the range that it was excited for PL. The .m file should explain the inputs and outputs in its comments.

QY_Calc calcualtes the quantum yield of a sample provided another reference sample. Again, the .m file is commented to explain its use.

PL_master incorporates the functionality of QY_Calc (with some improved plotting functionality) as well as converting PL data from wavelength (nm) to energy (ev), normalizing, plotting, and fitting of spectra to 1st order Gaussian functions at the user's request.

If you use these functions and run into any issues or desire more functionality, please do not hesistate to contact me with bugs/concerns/requests at MarkCBabin@gmail.com.
