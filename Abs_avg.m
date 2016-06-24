function abs_avg = Abs_avg(x,abs,wid,exc)
% This function averages the absorbance of a species around a central
% wavelength with a given width. This is to be used in conjuction with
% QY_Calc to more accurately determine the quantum yield of a samples
% The function arguements are:
% x       =   Wavelengths over which abs is measured
% abs     =   Optical density or abs of sample
% wid     =   Slit width of fluorometer
% exc     =   Excitation wavelength used in fluorometer

x_min = exc-wid./2;
x_max = exc+wid./2;
% Setting the minimum and maximum wavelengths to consider based on the
% slit width

x_bar1 = x.*(x>=x_min);
x_bar  = x_bar1.*(x_bar1<=x_max);
%Cutting the wavelength down to fit within x_min and x_max

wide = size(abs,2);
%this calcuates how wide (how many columns) abs is (when imported as a
%matrix

x_wide = repmat(x_bar,1,wide);
%creating a matrix of duplicates of the wavelength vector to match the size
%of the abs vector

abs_bi = abs.*logical(x_wide);
%we convert the x matrix to a binary (logical) matrix where every non-zero
%value becomes one, then multiply that with the abs input to isolate the
%absorbance values of interest

abs_cut = abs_bi(any(abs_bi,2),:);
%here we cut all rows full of zeros, leaving behind only the absorbance
%values we're interested in.


abs_avg = mean(abs_cut,1);
%finally, we average each column, giving us a vector of our average
%absorbanes for each value.