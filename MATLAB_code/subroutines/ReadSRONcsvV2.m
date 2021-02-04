function [header,data] = ReadSRONcsvV2(path,name,RamiHeaderCorruption)
%This function extracts the internal power, read-out power and temperature
%from the header of the Delft S21 datafiles. This can be done for various
%different headers based on the different devices used in the Delft setup.
%
% Does not nessecarily require a csv-file to be read.
%
% NOTE: For analysing the header of Delft noise files use the NoiseHeader function.
% NOTE: The predefined locations of this function might only be correct after
%       the importing by ReadSRONcsv.m
%
%INPUT:
% path == full path specifying the location of the file (ending with '\')
% name == name of the file including extension
% RamiHeaderCorruption == True (1) or False (0) [Optional, Default == 0]
%
% RamiHeaderCorruption Functionality:
% In the files written by the Delft Setup the labview output has a
% corruption (or at least something strange) making matlab take the first
% line of data as part of the header (Likely due to the lack of column
% names). This is compensated in this part of the routine.
% NOTE: Does not seem to be required anymore in the newest versions of
% matlab (tested in R2012a)
%
%OUTPUT:
% header == cell array containing the stringarrays of the header.
%               Typically: Cell vector with 1 cell per header line.
% data == double matrix containing the read data.
%
%SUBROUTINES:
%
%REQUIRED FILES:
%
%VERSION: 1.1
%   V1.0 original script by Reinier Janssen
%   V2.0(2012-08-07,RJ) Added comments, modified:
%                       header = 0 => header =[] in case of no header.
%
%DATE: August 7, 2012
%AUTHOR: Reinier Janssen

%==========================================================================
% Check if the optional input RamiHeaderCorruption is provided.
if nargin == 2
    RamiHeaderCorruption = 0;
end

%==========================================================================
%Simple import of the file
pathtot = strcat(path,name);
datacells = importdata(pathtot);
if strcmp(class(datacells),'struct')
    %If a header is present a struct will be created. The data and header
    %will then be extracted from different fields.
    data = getfield(datacells,'data');
    header = getfield(datacells,'textdata');
else
    %If a header is not present a regular matrix will be read. This is
    %immediately the data and an empty header is left.
    data = datacells;
    header = [];
end

%==========================================================================
%CORRECTION FOR DELFT HEADER - DATA SEPARATION ERROR
%==========================================================================
%In the files written by the Delft Setup the labview output has a
%corruption (or at least something strange) making matlab take the first
%line of data as part of (last line of) the header (Likely due to the lack 
%of column names). This is compensated in this part of the routine.
if RamiHeaderCorruption == 1
    %Correct their sizes
    datasize = size(data);
    headersize = size(header);
    datasize(1) = datasize(1)+1;
    headersize(1) = headersize(1)-1;
    
    %Copy the data into the correct positions
    CopyOfData = zeros(datasize);
    CopyOfData(2:end,:)=data;
    for p=1:datasize(2)
        CopyOfData(1,p) =  str2num(header{size(header,1),p});
    end
    data = CopyOfData;
    clear CopyOfData
    
    %Clean up the header
    headerold = header;
    clear header
    % Loop over all header components.
    for p=1:headersize(1) %corrected header size. 
        %Will not loop over last header line that is actually data.
        for q=1:headersize(2)
            % If there is something in the header, copy it to the new one.
            %Done this awkwardly because cells don't copy easily.
            if isempty(headerold{p,q}) == 0
                header{p,q} = headerold{p,q};
            end
        end
    end
    clear headerold
end
    