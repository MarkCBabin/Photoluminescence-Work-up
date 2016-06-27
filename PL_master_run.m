disp(sprintf('Use this function to open the data you would like to process with the PL_Master package.\n'))


addpath 'C:\Users\Mark\Documents\School\Undergraduate\Research\Roberts Lab\Analysis'
%this adds the path where PL_master, Abs_avg, and QY_Calc are stored. The
%use of quotations are only required if a folder name in the string has a
%space(as it does for me in "Roberts Lab" above). 

uiopen
%this will allow the user to select and open whichever files they so choose
%through the typical UI of their OS.

disp(sprintf('You are now able to run PL_Master, Abs_avg, and QY_Calc\n'));