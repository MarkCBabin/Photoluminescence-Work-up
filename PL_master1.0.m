%This function converts a PL spectrum from wavelength to energy (nm ->
%eV). It also has the ability to fit PL data with a guassian curve and 
%calculate the QY of measurements taken, provided a spectra of a known
%reference is included in the data. Prompts will be generated to ask if
%these functionalities are desired and give instruction on their
%implementation.
%
%PL_master should be used in conjuction with Abs_avg and PL_master_run to 
%determine the absorption average over the exciation wavelength, which will 
%be used as the input abs_avg (see below). Use/edit the script 
%PL_master_run to prevent having to copy this file into every folder used. 
%
%The function input arguements are:
% x            =       wavelengths over which PL is measured (in nm)
% PL           =       PL counts, background corrected
% abs_avg      =       avg abs at exciation wavelength (from abs_avg.m)
%The function output arguements are:
% ev           =       conversion of x from wavelength to ev
% PL_ev        =       PL spectra with Jacobian performed
% PL_norm      =       normalized PL spectra
%
%If the additional functionalities are performed the following variables
%will be output to your workspace:
% coeffs       =       coefficients of the gaussian fit (see details below)
% FWHM         =       full width at half max of each gaussian fit
% selected_x   =       x_values used in each gaussian fitting
% QY           =       calculated quantum yield of each sample
%
% coeffs will be a matrix of the form
%
%      a1  a1  a1  ...  a1
%      b1  b1  b1  ...  b1
%      c1  c1  c1  ...  c1
%
%where each column of a1, b1, and c1 correspeonds to the coefficients of a 
%gaussian fit of the form a1*e^(x-b1)^2/c1 for a particular column of data. 
%
%You may want to input a semi-colon after your function call in the
%comand line to prevent a long output. (See below)
%   [a,b,c] = PL_Conv(x,PL,abs_avg);
%
%Additionally, this file assumes that you input all of your PL data as a
%matrix wherein each point aligns with the same wavelength (i.e. data was
%collected over the same range with the same resoluton).


%% PL_master function; converting PL in nm to eV, normalizing, and plotting

function [ev,PL_ev,PL_norm] = PL_master(x,PL,abs_avg)
close all
%this will close all figure windows open

ev = 1239.84193./x;
%convert x (in wavelength) to energy (in eV)

wid = size(PL,2);
%this collects the width of the PL input matrix (i.e. if the matrix is 1x4,
%this value will be 4)

ev_wid = repmat(ev,1,wid);
%This creates a new 'x' vector in eV insead of wavelength that is of the
%same dimensions as PL (allowing you to input a matrix of PL data)

abs_avg1 = repmat(abs_avg,length(x),1);
%like before, this corrects an abs_avg input to be the correct dimensions
%for multiplication (.*)

PL_cor = PL./abs_avg1;
%corrects PL for absorbance at emission wavelength - this correlates one
%photon in for one photon out.

PL_ev = PL_cor*1239.84193./(ev_wid.^2);
%convert data by multiplying by the Jacobian (hc/E^2)

maxi = max(PL_ev);
%this finds the maximum value in each PL column, but returns it as a simple
%vector, which means that ./ will not work, hence the next step that
%converts our file to the same dimensions as our PL input.

max_long = repmat(maxi,length(x),1);

PL_norm = PL_ev./max_long;
%here we normalized each column in the PL matrix to its respective maximum
%value.

resp = cell(1,wid);
%this creates a cell array (basically a matrix of cells that can hold
%vectors - this is helpful because words are considered 1x n vectors). This
%will be used to store the inputs of the for loop below in which a user can
%enter the labels for each column in the data matrix which will end up on
%the legend on the final plots.

for i = 1:wid
    str = sprintf('What is the name of column %d?\n', i);
    resp{i} = input(str, 's');
end

figure;
plot(ev,PL_ev)
title('Plot of PL Spectra')
xlabel('Energy (eV)')
ylabel('Counts')
legend(resp)
figure;
plot(ev,PL_norm)
title('Plot of Normalized PL Spectra')
xlabel('Energy (eV)')
ylabel('Normalized PL')
legend(resp)
%this will plot both your normalized spectra (Figure 2) and the raw spectra
%in eV (Figure 1) with a legend comprised of the entries to the above
%prompt, saved in the "resp" array.

prompt1 = 'Would you like to fit your data with Gaussians?\n If so, type any number.\n If not, hit enter.\n';
str1 = input(prompt1);
if isempty(str1)
    prompt2 = '\nWould you like to calculate the quantum yield of these samples?\n If so, enter which column (number) your QY standard is in.\n If not, just hit enter.\n';
    str2 = input(prompt2);
    if isempty(str2)
        fprintf('Thanks, you are done!\n')
        return
    else
        QY_Calc1(x,PL,abs_avg,wid,str2,resp);
    end
else
    Gaussfitting(ev,PL_norm,wid,resp)
    prompt2 = '\nWould you like to calculate the quantum yield of these samples?\n If so, enter which column (number) your QY standard is in.\n If not, just hit enter.\n';
    str2 = input(prompt2);
    if isempty(str2)
        fprintf('Thanks, you are done!\n')
        return
    else
        QY_Calc1(x,PL,abs_avg,wid,str2,resp);
    end
end
     
%This decision tree of if statesments allows the user whether or not they
%would like to run the next two subfuncitons. In the case of QY_Calc1,
%these prompts also generate some of the input necessary for the function
%to run. 
end 

%% Gaussian Fitting of data and determination of fwhm/fit quality

function Gaussfitting(ev,PL,wid,resp)
fprintf('\n\nUse the mouse to select the bounds over which the program will fit your data with a Gaussian.\nPlease select traces in the order they appear in the legend (top to bottom)\nand please click form LEFT to RIGHT.\nIt may help to maximize the plot on your screen.\n')

[xx,yy] = ginput(2.*wid);
%this creates two vectors (xx and yy) that correspond to the x and y
%coordinates selected by the user on Figure 2 using their mouse.

coeffs = zeros(3,wid);
FWHM = zeros(1,wid);
%here, we create the vectors that will be populated in the for loop below.
%It is important to create these beforehand to improve the speed of the
%program, particularly if you will be iterating through many data sets.

for i = 1:wid
    y_min_index = 1+(length(PL).*(i-1));
    y_max_index = 0+(length(PL).*i);
    x_min_index = 2.*i-1;
    x_max_index = 2.*i;
    domain1 = [xx(x_min_index),xx(x_max_index)];
    if wid == 1
        outliers = excludedata(ev,PL(y_min_index:y_max_index),'domain',domain1);
        test1 = fit(ev,PL(y_min_index:y_max_index),'gauss1','Exclude', outliers);
    else 
        outliers = excludedata(ev,PL(y_min_index:y_max_index)','domain',domain1);
        test1 = fit(ev,PL(y_min_index:y_max_index)','gauss1','Exclude', outliers);  
    end
    coeff = coeffvalues(test1)';
    a1 = coeff(1);
    b1 = coeff(2);
    c1 = coeff(3);
    coeffs(1,i) = a1;
    coeffs(2,i) = b1;
    coeffs(3,i) = c1; 
    FWHM(i) = 2.*sqrt(2.*log(2)).*c1;
    resp_i = resp{i};
    if wid == 1
        figure;
        subplot(2,1,1);
        plot(test1,ev,PL(y_min_index:y_max_index),'Residuals')
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit1 = sprintf('Residuals of the fit for %s',resp_i);
        title(tit1,'FontSize',15)    
        subplot(2,1,2);
        plot(test1,ev,PL(y_min_index:y_max_index))
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit2 = sprintf('Plot and fit of %s, with FWHM of %s',resp_i,c1);
        title(tit2,'FontSize',15)      
    else
        figure;
        subplot(2,1,1);
        plot(test1,ev,PL(y_min_index:y_max_index)','Residuals')
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit1 = sprintf('Residuals of the fit for %s',resp_i);
        title(tit1,'FontSize',15)    
        subplot(2,1,2);
        plot(test1,ev,PL(y_min_index:y_max_index)')
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit2 = sprintf('Plot and fit of %s, with FWHM of %s',resp_i,c1);
        title(tit2,'FontSize',15)
    end
end

%here we fit each set of data to a gaussian function using the bounds
%slected by the user
%The x/y_min/max  indicies indicate the boundries over which to fit
%They y's are required for a data matrix so that each column is treated
%individually while the x's are generated from the user selected values
%from here, we create a subplot system that plots both the PL data against
%its fit as well as a plot of the residuals for each fit. The titles are
%generated based on the index of the "resp" variable which is generated in
%the main function by user input to give each column in the matrix a name
%for legened and plotting purposes.

fprintf('The next displays will repeat until you are satisfied with your fittings.\n')

prompta = '\nWould you like to re-select the x-boundries for your Gaussian fitting?\n If so, enter the number of the FIGURE you would like to correct.\n If not, type zero.\n';
whil_index = input(prompta);

while whil_index >= 3 
    colnum = whil_index - 2;
    disp(figure(whil_index))
    [xx(2.*colnum -1),yy(2.*colnum-1)] = ginput(1);
    [xx(2.*colnum),yy(2.*colnum)] = ginput(1);
    y_min_index = 1+(length(PL).*(colnum-1));
    y_max_index = 0+(length(PL).*colnum);
    x_min_index = 2.*colnum-1;
    x_max_index = 2.*colnum;
    domain1 = [xx(x_min_index),xx(x_max_index)];
    if wid == 1
        outliers = excludedata(ev,PL(y_min_index:y_max_index),'domain',domain1);
        test1 = fit(ev,PL(y_min_index:y_max_index),'gauss1','Exclude', outliers);
    else 
        outliers = excludedata(ev,PL(y_min_index:y_max_index)','domain',domain1);
        test1 = fit(ev,PL(y_min_index:y_max_index)','gauss1','Exclude', outliers);  
    end
    aaa = coeffvalues(test1)';
    aa = aaa(1);
    bb = aaa(2);
    cc = aaa(3);
    coeffs(1,colnum) = aa;
    coeffs(2,colnum) = bb;
    coeffs(3,colnum) = cc;
    FWHM(colnum) = 2.*sqrt(2.*log(2)).*cc;
    resp_c = resp{colnum};
    if wid == 1
        figure(whil_index);
        subplot(2,1,1);
        plot(test1,ev,PL(y_min_index:y_max_index),'Residuals')
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit1 = sprintf('Residuals of the fit for %s',resp_c);
        title(tit1,'FontSize',15)    
        subplot(2,1,2);
        plot(test1,ev,PL(y_min_index:y_max_index))
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit2 = sprintf('Plot and fit of %s, with FWHM of %s',resp_c,c1);
        title(tit2,'FontSize',15)      
    else
        figure(whil_index);
        subplot(2,1,1);
        plot(test1,ev,PL(y_min_index:y_max_index)','Residuals')
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit1 = sprintf('Residuals of the fit for %s',resp_c);
        title(tit1,'FontSize',15)    
        subplot(2,1,2);
        plot(test1,ev,PL(y_min_index:y_max_index)')
        xlabel('Energy (eV)','FontSize',15)
        ylabel('Counts','FontSize',15)
        tit2 = sprintf('Plot and fit of %s, with FWHM of %s',resp_c,c1);
        title(tit2,'FontSize',15)
    end
%This while loop gives users the options to alter the bounds over which
%the fit is computed. It is done iteratively, allowing the user to select
%which plot to re-fit, and can continue until the user is happy with the
%fit generated. This will over-write the values of "coeffs" and "FWHM" as
%well as the figure that is being corrected.

prompta = '\nWould you like to re-select the x-boundries for your Gaussian fitting?\n If so, enter the number of the FIGURE you would like to correct.\n If not, type zero.\n';
whil_index = input(prompta);
end

assignin('base','coeffs',coeffs);
assignin('base','FWHM',FWHM);
assignin('base','selected_x',xx);
%here we save the values of the selected x values (from ginput above), the 
%calculated full width at half max, and the fit coefficients to the main
%workspace so that the user can see them and use them outside of this
%function.
end

%% Option to calc QY from this data
function QY_Calc1(x,PL,abs_avg,wid,str2,resp)

prompt3 = 'Please enter the index of refraction for the solvents for sample and reference, respectively\n Please use the format [n_sample,n_ref]\n';
str3 = input(prompt3);
prompt4 = 'Now please enter the value of your referene QY:\n';
str4 = input(prompt4);
%this if statement either terminates the code (if there is no input) or
%generates two new prompts whose inputs are treated as numbers rather than
%strings (of charactrs). The use of input(x) instead of input(x,'s') causes
%the change.

area = trapz(x,PL);
%This step calculates the area of each plot (numeric intergration using
%trapazoids).

trans_ref = 1-10^(-1.*abs_avg(str2));
trans_samp = zeros(1,wid);
QY = zeros(1,wid);

for i = 1:wid
    trans_samp(i) = 1-10^(-1.*abs_avg(i));
    if i == str2
        QY(i) = str4.*(area(i)./area(str2))*(trans_ref./trans_samp(i))
    else
        QY(i) = str4.*(area(i)./area(str2))*(trans_ref./trans_samp(i))*(str3(1)./str3(2))^2
    end
end
%here, we calculate the QY for each sample iteratively
%This equation is pulled from "Standards for Photoluminescence Quantum 
%Yield Measurements in Solution (IUPAC Technical Report)" by Brouwer (2001)

assignin('base','QY',QY);
%here, the QYs are assigned a variable in the main workspace (base) as the
%name "QY" as there is no output variable in the function header for QY
%(this is because QY is an optional output, so putting it as a standard
%output causes an error message.

hold on
figure
bar(QY)
ylabel('QY')
set(gca,'XTick',1:wid)
set(gca,'XTickLabel',resp)
hold off

%Now we plot a bar graph of these QYs and label the x-axis with the names
%given in resp (this is acheived by setting the number of x-ticks to
%force the bar plot to allow a manual changing of the labeling, which is
%done in the second "set" command

fprintf('PL Master has completed its functionality.\n')
end
