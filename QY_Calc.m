function y = QY_Calc(x,y_samp,y_ref,n_samp,n_ref,abs_samp,abs_ref,QY_ref)
% QY(x,y_samp,y_ref,n_samp,n_ref,abs_samp,abs_ref,QY_ref) calculates the
% quantum yield of a sample (samp) by using a reference (ref). 
% The function arguements are:
% x        =   Wavelengths PL are measured over for both ref and samp
% y_samp   =   PL counts for the sample
% y_ref    =   PL counts for the reference
% n_samp   =   Index of refraction of sample's solvent
% n_ref    =   Index of refraction of reference's solvent
% abs_samp =   UV-Vis absorbance of sample at excitation wavelength
% abs_ref  =   UV-Vis absorbance of reference at excitation wavelength
% QY_ref   =   Quantum Yield of reference

Int_ref = trapz(x,y_ref);
Int_samp = trapz(x,y_samp);
Trans_ref = 1-10^(-1.*abs_ref);
Trans_samp = 1-10^(-1.*abs_samp);

%Here, the function integratss the area under each PL plot (using the trapz
%funtion)

y = QY_ref.*Int_samp./Int_ref.*Trans_samp./Trans_ref.*(n_samp./n_ref)^2


