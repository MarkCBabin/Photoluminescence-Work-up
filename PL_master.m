%NOTE: You may want to input a semi-colon after your function call in the
%comand line to prevent a long output. (See below)
% [a,b,c] = PL_Conv(x,PL,abs_avg);

% Additionally, this file assumes that you input all of your PL data as a
% matrix wherein each point aligns with the same wavelength (i.e. data was
% collected over the same range with the same resoluton).

%Lastly, you will need to copy this file (and abs_avg.m if you plan on
%using that) into the same folder as your data for it to work properly.

function [ev,PL_ev,PL_norm] = PL_master(x,PL,abs_avg)
%This function converts a PL spectrum from wavelength the energy (nm ->
%eV). It should be used in conjuction with abs_avg to determine the
%absorption average over the exciation wavelength.

%The function input arguements are:
% x            =         wavelengths over which PL is measured
% PL           =         PL counts, background corrected
% abs_avg      =         avg abs at exciation wavelength (from abs_avg.m)
% The function output arguements are:
% ev           =         conversion of x from wavelength to ev
% PL_ev        =         PL spectra with Jacobian performed
% PL_norm      =         normalized PL spectra

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

PL_cor = PL.*abs_avg1;
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
%vectors - which is helpful because words are considered 1x"n" vectors). This
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



%% Gaussian Fitting of data and determination of fwhm/fit quality

disp(sprintf('\n\nUse the mouse to select the bounds over which the program will fit your data with a Gaussian.\nPlease select traces in the order they appear in the legend (top to bottom)\nand please click form LEFT to RIGHT.\nIt may help to maximize the plot on your screen.'));

[xx,yy] = ginput(2.*wid);
%this creates two vectors (xx and yy) that correspond to the x and y
%coordinates selected by the user on Figure 2 using their mouse.

assignin('base','selected_x',xx);
%here we save the values of the selected x values above to the main
%workspace so that the user can see them and use them outside of this
%function.

for i = 1:wid
    y_min_index = 1+length(PL).*(i-1);
    y_max_index = length(PL).*i;
    x_min_index = 2.*i-1;
    x_max_index = 2.*i;
    domain1 = [xx(x_min_index),xx(x_max_index)];
    outliers = excludedata(ev,PL(y_min_index:y_max_index)','domain',domain1);
    test1 = fit(ev,PL(y_min_index:y_max_index)','gauss1','Exclude', outliers);
    figure;
    plot(test1,ev,PL(y_min_index:y_max_index)')
    xlabel('Energy (eV)')
    ylabel('Counts')
    resp_i = resp{i};
    titl = sprintf('Plot and fit of %s',resp_i);
    title(titl)
end

%here we fit each set of data to a gaussian function using the bounds
%slected by the user
%The x/y_min/max  indicies indicate the boundries over which to fit
%They y's are required for a data matrix so that each column is treated
%individually while the x's are generated from the user selected values

















%% Option to calc QY from this data
prompt = 'Would you like to calculate the quantum yield of these samples?\n If so, enter which column (number) your QY standard is in.\n If not, just hit enter.\n';
str = input(prompt);
if isempty(str)
    disp(sprintf('Thanks, you are done!'))
    return
else
    prompt2 = 'Please enter the index of refraction for the solvents for sample and reference, respectively\n Please use the format [n_sample,n_ref]\n';
    str2 = input(prompt2);
    prompt3 = 'Now please enter the value of your referene QY:\n';
    str3 = input(prompt3);
end
%this if statement either terminates the code (if there is no input) or
%generates two new prompts whose inputs are treated as numbers rather than
%strings (of charactrs). The use of input(x) instead of input(x,'s') causes
%the change.

area = trapz(x,PL);
%This step calculates the area of each plot (numeric intergration using
%trapazoids).

ref_area = area(str);
%this calcualtes the integrated area of the reference

for i = 1:wid
    trans_ref = 1-10^(-1.*abs_avg(str));
    trans_samp(i) = 1-10^(-1.*abs_avg(i));
    QY(i) = str3.*area(i)./area(str)*trans_samp(i)./trans_ref*(str2(1)./str2(2))^2;
end
%here, we calculate the QY for each sample iteratively

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


disp(sprintf('PL Master has completed its functionality'))
%Now we plot a bar graph of these QYs and label the x-axis with the names
%given in resp (this is acheived by setting the number of x-ticks to
%force the bar plot to allow a manual changing of the labeling, which is
%done in the second "set" command
