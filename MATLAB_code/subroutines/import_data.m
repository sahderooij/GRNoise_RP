function [Data,Temperature,Power,Header] = import_data(filename)
% This is a multipurpose import function that can read any data with the
% format:
% * A Header of a variable number of lines (max 50)
% * 2 empty lines (uses this for header end recongnition)
% Any number of the following set of temperature/data
% * Subheader containing in the first line the temperature and in the
% second line the column headers
% * The actual data in up to 10 columns
% * 1 empty line
% An example with 3 data columns is given below
% 
% This function will extract from the header the read power at the KID
% (requires the string 'Power at KID:') and from each subheader the
% temperature (requires the first line of the subheader to be 'Temperature
% in K:xx.xxx').
% 
% 
% EXAMPLE FORMAT:
% C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
% Power at KID:73dBm
% resonance Frequency in GHz :4.893858
% Q=40792.348421, Qc=56509.887020, Qi=146662.340671  
%
%
% Temperature in K:0.099953
% GHz	dB	Rad
% 4.893115000	-6.722063065	2.685341760
% 4.893118870	-6.721852779	2.685013925
% 4.894663000	-6.540477276	2.794459459
%
% Temperature in K:0.105161
% GHz	dB	Rad
% 4.893118740	-6.718620300	2.685038958
% 4.893122571	-6.718233109	2.684799274
%
%
%INPUT:
%   filename = fullpath of the file.
%
%OUTPUT:
% For nT the number of temperatures (subsections) in the file
%   Data = a cell vector of length nT, containing double arrays with the
%          data of each temperature.
%   Temperature = a double vector of length nT, containing the temperatures
%                   found in the subheaders
%   Power = a double giving the Read Power of the KID as found in the
%           header.
%   Header = cell vector containing the header lines
%
%SUBROUTINES:
%
%REQUIRED FILES:
%   SRON data file in the format given above.
%
%VERSION: 1.0
%
%DATE: October 29, 2012
%AUTHOR: Reinier Janssen
%==========================================================================

%==========================================================================
% Set internal variables
%==========================================================================
format('long','e'); %set display format
nTmax = 200; %Maximum number of temperatures in the read file

%==========================================================================
% Open file and read header
%==========================================================================
fid = fopen(filename);
if fid == -1
    error(['ERROR import_data: Cannot find ' filename ])
else
    disp(['Reading: ',filename])
end

%First read through the header line by line until two lines are found
emptylines = 0; %Counter counting the number of consequtive empty lines
p=0; %Header line counter
Power = 0/0; %Step a default value for Power.
while emptylines < 2
    p=p+1;
    Header{p,1} = fgetl(fid); %gets the full next line of the header.
    if isempty(Header{p,1}) %Check if the line is empty
        emptylines = emptylines +1;
    else
        emptylines = 0; %reset emptylines\
        %Check for read power and extract
        PowerIndex = strfind(Header{p,1},'Power at KID:');
        if ~isempty(PowerIndex)
            PowerString = Header{p,1}(PowerIndex+13:end);
            Power = cell2mat(textscan(PowerString,'%f'));
        end
    end
    
    if p>50
        %Error catching for if double empty line is not found.
        error('ERROR import_data: Too many header lines.\n')
    end
end

%==========================================================================
% Loop over all subheader/data combinations and recover data
%==========================================================================
%After the header we start looping over subheader/data combinations
nT = 0; %temperature counter
tempTemperature = zeros(nTmax,1); %Temporary Temperatures vector
tempData = cell(nTmax,1); %Temporary Data Cell vector

while nT < nTmax
    temp = textscan(fid,'%*s%*s%*2c%f%*s',1);
    if isempty(temp{1})
        break;
    else
        nT = nT +1;
    end
    tempTemperature(nT,1) = temp{1};
    temp = cell2mat(textscan(fid,'%f%f%f%f%f%f%f%f%f%f','headerlines',2));
    if ~isempty(temp) 
        dataindexes = ~isnan(temp(1,:));
        tempData{nT,1} = temp(:,dataindexes);
    else
        tempData{nT,1} = [];
    end
end
fclose(fid);
%remove excess variables
clear temp

Temperature = tempTemperature(1:nT,1);
Data = tempData(1:nT,1);
%==========================================================================
%END OF ROUTINE
%==========================================================================
end