function [KID,KIDParameters] = S21analysisV5
% JOCHM 10-2014 Made KID circle callibration more reliable for deep KIDs.
% no version change

% This function will search in the specified directory for all (KID,Pread)
% combinations and read the S21(T) sweep. It will then process all S21(F)
% sweeps to obtain Fres,Q and other KID specific parameters, as well as
% calibrate the resonance circle. After using Kaplan to convert T into
% #qp's using various material parameter, responses and responsivities are
% determined. Finally all information is stored to figure and csv files.
%
% For good practice always: CHECK DELTA AND RESPONSE VALUES
%
% NOTE: 2D (power and temperature) or 1D power data is handled by treating 
%   each (KID,power) combination as a new KID. Only the T dependence of
%   (KID,min(P)) is shown to screen unless 'showallplots' is set to 1 on
%   line 73. Further Power Dependence can be analysed using
%   PowerDependenceS21.m
%
% NOTE: Function assumes External subroutines are located in the 
%   ..\subroutines subdirectory, while this function is in the main
%   directory.
%
% NOTE: (TODO) These modifications could be made in the future:
%   *stopoffset variable for each resonator
%   *replace subplot(3,3,4) or (3,3,5) with a graph of Fres(T)
%
%INPUT:
%   SPECIFIED INTERNALLY
%
% S21subdir == subdirectory specifying the location of the S21 files.
%               Starts with "\", ends without "\"
% 
% ChipInfo == a struct containing information about the measured chip in 
%           the following elements
%   ChipInfo.path == Root path containing the measurement data.
%   ChipInfo.ResonatorInfoFile == Name (Location) of the file containing information
%                                 about each resonator (see below for the
%                                 layout). NOTE: include the '.dat' or
%                                 '.txt'. 4 cols: KID# Fres[GHz] Q Al area
%                                 [um^2]
%   ChipInfo.FilmTc == Transition temperature measured for the film [K].
%   ChipInfo.material == Material of which the resonator active area is made.
%                        Supported Materials are [Delta/kbTc,N0]:
%                        * 'Al' [1.6,1.70e10]
%                        * 'Ta' [1.8,4.08e10]
%                        * 'Ti' [1.76,1.70e10 WRONG]
%                        * 'TiN' [1.81,8.7e9]
%                        * 'AlBCS' (Default) [1.76,1.70e10]
%                        * 'AlPascale' [2.01,1.70e10]
%                        * 'NbTiN'[1.76,3.7e10]
%                   NOTE: Case Sensitive
%   ChipInfo.Thickness == Thickness of the resonator active area [um].
%
% EXTRA: stopoffset, showallplots, TCfrac, FITF0. (See below).
%==========================================================================
% Set ChipInfo inside function.
%==========================================================================
%
%==========================================================================
%OUTPUT:
%   KID(N). == struct array of N = # KIDs, which for all 
%                     temperatures (nT = # temperatures).
%       KID(N).KIDnumber = KIDnumber (ID)
%       KID(N).filename = a single string containing the full
%                         path+filename that contains the S21datafile.
%       KID(N).Temperature(nT,1) = vector containing the chip temperature [K]
%       KID(N).ReadPower(nT,1) = vector containing the used read power [dBm]
%       KID(N).S21data{nT,1}(:,1:5) = cell array containing the S21data at each Temperature
%              Each cell contains a double array with the columns:
%              [F (GHz), |S21|, phase (rad), Re S21 (normalized), Im S21 (normalized)]
%       KID(n).Ql(nT,1) == Measured Loaded Q at each temperature [From S21 fit]
%       KID(n).Qi(nT,1) == Measured Internal Q at each temperature [From S21 fit]
%       KID(n).Qc(nT,1) == Measured Coupling (External) Q at each temperature [From S21 fit]
%       KID(n).Fres(nT,1) == Measured Resonance Frequency [From S21 fit, subsampling resolution]
%       KID(n).S21min(nT,1) == Minimum value S21 [From S21 fit, subsampling resolution]
%       KID(n).Paramerrors(nT,1:5) == errors in [Fres,S21min(dB),Ql,Qi,Qc] obtained from S21 fitting.
%       KID(n).InternalPower(nT,1) == KID internal power [dBm]
%       KID(n).Nqp(:,1) == Number of thermally created quasiparticles
%       KID(n).Delta == Assumed superconducting gap [J]
%       KID(n).Volume == Volume of the active material [um^3]
%       KID(n).Tc == Assumed transition temperature [K]
%       KID(n).Area == Designed area of the active material [um^2]
%       KID(n).Thickness == Thickness of the active material [um]
%       KID(n).Fdesign == Resonator Design Frequency [GHz]
%       KID(n).Qdesign == Resonator Design Q (Qc)
%       KID(n).ddxdNqp == 2 element vector containing the coefficients of a
%                          linear fit to (Fres(T)-Fres(0))/Fres(0) vs Nqp
%       KID(n).dBdNqp == 2 element vector containing the coefficients of a
%                          linear fit to 1/Qi(T) vs Nqp
%       KID(n).ResponsivityM1(:,1:2) == phase (1) and radius (2) responsivity
%                                       as a function of temperature as
%                                       determined by Method 1 (see above)
%       KID(n).ReImF0(:,1:2) == The value of S21(Re(1),Im(2)) in the calibrated
%                               resonance circle at Fres(T=0[K])
%       KID(n).Response(:,1:2) == The (change in) phase (2) and radius (1) at the 
%                                 base temperature resonance frequency wrt
%                                 the circle at T == 0 [K]
%       KID(n).responsivityM2(:,1:2) == phase (1) and radius (2) responsivity
%                                       as a function of temperature as 
%                                       determined by Method 2 (see above)
%   KIDParameters == Nx16 array contains the most important information 
%                  about all N=#(KID,P) combinations at base temperature
%                  It contains the parameters (Dim 2):
%               1 = KIDID
%               2 = Tbase [K]
%               3 = Fres(Tbase) [GHz]
%               4 = Fdesign [GHz]
%               5  = Q
%               6 = Qi
%               7 = Qc
%               8 = Qdesign
%               9 = Pread [dBm]
%               10 = Pint [dBm]
%               11 = ddx/dNqp
%               12 = dR/dNqp'
%               13 = Tc [K]
%               14 = Delta [J]
%               15 = Film Thickness [um]
%               16 = Active Area [um^2]
%
%   In addition all results are output as 4 types of files (Described below):
%   In the S21subdir where all measurements are:
%   One overview file with information of all KIDs at the base temperature.
%   One file per KID containing S21(F) at base temperature
%   One file per KID containing analysis results at all measurement temperatures.
%   In the measurement main directory:
%   A summary file containing 16 columns with the values of KIDParameters.
%
%
%SUBROUTINES:
% import_S21data
%   import_data(External)hybri
% findparameters
% FitS21main5 (External)
% QPnumber
% findresponse
% savethedata
% 	WriteSRONcsv (External)
% plotTdependence
%   GenerateColorsFromMap (External)
%       colormapstorage.mat (External)
%
%REQUIRED FILES:
%
%VERSION: 4.3
%   Based upon the original "process_S21_2" program
%   V1.0 original script by Stephen Yates & Jochem Baselmans
%   V1.1(2011-03-16,JB) Added ALPascale option, using delta/lTc = 2.10
%   V1.2(2012-08-03,RJ) clean-up and addition of comments.
%                       Changed code (not effect) in which unique KIDnumbers
%                       are determined.
%   V2.0(2012-08-21,RJ) Addition of more comments, switch to FitS21_2
%                       routine,
%   V2.0(2012-09-05,RJ) Fixed a 2011b compatibility issue of matlabs unique
%                       function
%   V2.1 (2012-10-18,RJ) Made is mac/unix compatible by using filesep
%                       instead of /. Output Fres now in Hz (ipv GHz) for
%                       backwards compatibility.
%   V2.2 (2012-10-26,RJ) Updated some comments
%   V4.0 (2012-10-30,RJ) Modified the import_S21data function to use the
%                        general import_data subroutine. 
%   V4.1 (2012-11-07,RJ) Renamed to from HybridsS21 to RESPONSEanalysis.
%                       Added a clear statement to remove useless variables
%                       from the workspace before saving.
%   V4.2 (2013-02-14,RJ) Now uses FitS21main3
%   V4.3 (2013-03-06,RJ) Now made robust against poor data in which no KID
%                       can be found. (It will give warning but will not
%                       crash the code anymore).
%   V5.0 (2013-11-20,RJ) Now uses FitS21main5 en stores error estimates
%DATE: September 5, 2012
%AUTHOR: Reinier Janssen
%
%==========================================================================
% Layout of user INPUT defined files
%==========================================================================
% ResonatorInfoFile: A four column - tab separated - file without header
% lines containing for each KID the following information.
% ID    Fres[GHz]   Qdesign     Active Area [um^2]
% 17	4.2         14000   	42156
%
% KID S21 measurements (SRON): A single file per KIDnr-Power combination.
% Temperatures are added sequentially in these files.
%   EXAMPLE LAYOUT SRON S21 FILE
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=40792.348421, Qc=56509.887020, Qi=146662.340671  
%
%
%Temperature in K:0.099953
%GHz	dB	Rad
%4.893115000	-6.722063065	2.685341760
%4.893118870	-6.721852779	2.685013925
%4.894663000	-6.540477276	2.794459459
%
%Temperature in K:0.105161
%GHz	dB	Rad
%4.893118740	-6.718620300	2.685038958
%4.893122571	-6.718233109	2.684799274
%   END EXAMPLE LAYOUT SRON S21 FILE
%==========================================================================
% Layout of OUTPUT files
%==========================================================================
% KID_KIDid_PdBm_Tdep.csv: A 26 column - comma separated - file with header line
% containing for a KID at a single power the following parameters:
% 1=KIDID,    2=Temperature,	3=Q,	4=Qc,	5=Qi
% 6=Fres,   7=S21@Fres,	8=power,	9=Pinternal,	10=response
% 11=dtheta/dNqp FIT,	12=dtheta/dNqp numerical,     13=Nqp,	14=dx/dN,	15=V[um3]
% 16=Delta [J],	 17=Re point,	18=Im point,	19=drdnqp=ampresp,	20=dRdB
% 21=grad 1/Qi vs nqp, 22=Tc,   23=Fres des, 24=Q des, 25=Area des,
% 26=thickness
%
% (response = )
% (dRdB = numerical variation of radius with 1/Qi, normalized to lowest T)
%
% 
% KID_KIDid_PdBm_S21lowT.csv: A 5 column - comma separated - file with
% header line containing for a KID at a single power the transmission as a
% function of frequency at the lowest available temperature
% Freq[GHz],S21[dB],theta[rad],Re Norm, Im Norm
% 5.754142630000e+00,-3.573115700000e-02,-6.639363800000e-02,4.806554635962e-01,1.857603917372e-02
%
% allTdep.csv: A 26 column - comma separated - file containing the same
% parameters as KID_KIDid_PdBm_Tdep.csv. However, this file contains information
% of all analysed KIDs but only at base temperature.
%==========================================================================

%sub-dir inside ChipPath where the S21 T dependence is located. Starts with \ ends without \.
S21subdir='/S21/2D';                                          
%Note: pwd is current matlab directory. Usually directory where this m-file
%is located.
if nargin==0
    ChipInfo.path = ...                                         %Root path containing the measurement data.
    '../data/A1A2_test';
    ChipInfo.ResonatorInfoFile = ...                            %Name of the file containing resonator information
        [ChipInfo.path '/KIDparam.txt'];
    ChipInfo.FilmTc = 1.255;                                    %Transition temperature measured for the film [K].
	ChipInfo.material = 'Al';                                %Material of which the resonator active area is made.
    ChipInfo.Thickness = 0.050;                                  %Thickness of the resonator active area [um].
    ChipInfo.Linewidth = 1.73; %linewidht in um, assumes that 4th colum of txt is length, not area
end

%==========================================================================
% Internal variables which may be user modified.
%==========================================================================
global TCfrac FITF0 InterpolatPhaseref
stopoffset = ones(30).*0;    %Removes the highest temperatures from the data analysis. The index correponds to KID number.
                             %Used to stop before the measurement has lost any KID due to dipjumping, extremely shallow dips, or similar.
showallplots = 0; %IF 1 the T dependences will be shown in plot for all (KID,P) combinations.
                  %IF 0 [default], only for the minimum read-out power the T dependence will be shown on screen.
                  %The other figures will be save and closed.
TCfrac=1/5;      %Fraction of Tc from where the responsivity fit starts. Default: 1/5
FITF0=1;          %Switch used to determine which method is used for determination of F0. 
                  %Please see the FitS21main5 routine for details.
% 0 == use minimum value of measured |S21| to determine fres and S21min.
%       Fit is only used to determine Q.
% 1 == use a logspace Lorentzian Fit to obtain Q, 
%       Fres and S21min. (FitS21_3 by Pieter de Visser). 
% 2 == use the function by Khalil et.al. to fit to |S21|(dB) (FitS21_2)
InterpolatPhaseref=1;    %1=find phaseref by interpolartion (default), best for VNA deep KIDs. 0 takes phase at measured datapoint closes to Fres.
%==========================================================================


%==========================================================================
% Setting some routine default values
%==========================================================================
global S21path
format('long','e'); %Set display format of numbers to 7 digits
%Enable subroutines by adding path in search path.
addpath([pwd,filesep,'subroutines']);
%close all; %close all open plots to remove clutter.
%==========================================================================
% Read in KID design values from the ResonatorInfoFile
%==========================================================================
fid=fopen(ChipInfo.ResonatorInfoFile);
if fid <0
    error(['Cannot find ' ChipInfo.ResonatorInfoFile] )
end
designvalues=cell2mat(textscan(fid,'%f %f %f %f'));%kidID F0 Q active area in um^2
designvalues(:,4) = designvalues(:,4)*ChipInfo.Linewidth;
fclose(fid);
%==========================================================================
% In the datapath find all files and identify KIDs.
%==========================================================================
S21path=[ChipInfo.path S21subdir filesep]; %Path containing raw S21 data
RawS21files=dir([S21path 'KID*.dat']); %Reading the names of all the raw S21 files into an array of cells
NumberOfS21Files = length(RawS21files);

KIDnumbers = zeros(1,NumberOfS21Files);
for n=1:NumberOfS21Files %Loop over all S21 files
    %Determine the KIDid (its number, XX) from each filename.
    KIDnumbers(n) = cell2mat(textscan(RawS21files(n).name,'%*3s %f %*s'));
end
%Determine all unique KIDs.
KIDnumbers = unique(KIDnumbers);

%Print to screen some of the reading information.
fprintf('Search for S21 data performed in:\n')
disp(S21path)
fprintf(['Inside the path the following KIDs are available: ',num2str(KIDnumbers),'\n'])
fprintf(['Inside this path a total number of ',num2str(NumberOfS21Files),' files were found.\n'])

%==========================================================================
% Analyse each S21 file. For each raw file a results file is
% created with a similar name to the S21data file, but with the
% ''extension.
% All results save in KID struct array.
% A graph is only given if a single power per KID is measured. Graphs are
% saved to disk in any case.
%==========================================================================

%Predefine some parameters:
transformation = cell(NumberOfS21Files,1);
%Stores for each (KID#,P) combination the calibration transformation of the 
%resonance circle in the complex plane
fitTbase = cell(NumberOfS21Files,1);
%Stores for each (KID#,P) combination the [F(GHz),|S21|(dB)] that is the
%evaluated version of the fit to the base temperature S21 transmission.
for n=1:NumberOfS21Files %Loop over all S21 files
    %Extract from the filename piece of information: KIDnumber
    info=textscan(RawS21files(n).name,'%*3s %f %*1s %*f %*s'); %Get KIDnumber
    %Note: second float is the power [dBm]. Replace %*f by %f to get Pread = info{2}
    KID(n).KIDnumber = info{1}; %Storing KIDnumber
    KID(n).filename = [S21path,RawS21files(n).name]; %Storing original data file
    %======================================================================
    
    %Import the S21 data from the specified file.
    [KID(n).S21data,KID(n).Temperature,KID(n).ReadPower] = import_S21data(KID(n).filename,stopoffset(info{1}));
    format short
    if n == 1
        disp([(1:1:length(KID(n).Temperature))' KID(n).Temperature])
    end
    format long
    % Assigns:
    %   KID(n).S21data{:,1}(:,1:3)
    %   KID(n).Temperature(:,1)
    %   KID(n).ReadPower(:,1)
    %======================================================================
    
    %Analyse the S21data to obtain:
    
    [KID,transformation{n},fitTbase{n}]=findparameters(n,KID);
    % Requires: KID(n).Temperature, KID(n).S21data{:,1}(:,1:3), KID(n).ReadPower
    % Assigns: 
    %   KID(n).S21data{:,1}(:,2) [update, normalized & converted to magnitude space]
    %   KID(n).S21data{:,1}(:,4:5)
    %   KID(n).Ql(:,1)
    %   KID(n).Qi(:,1)
    %   KID(n).Qc(:,1)
    %   KID(n).Fres(:,1) [From S21 fit, subsampling resolution]
    %   KID(n).S21min(:,1) [From S21 fit, subsampling resolution]
    %   KID(n).Paramerrors(:,1:5)
    %   KID(n).InternalPower
    %   
    %   transformation{n} = 2 element vector containing the mean (over T) phase
    %   rotationa and Real axis translation applied to the resonance circle
    %   for calibration.
    %   fitTbase{n} = FitResult of the fit performed by FitS21 subroutine.
    %======================================================================  
    
    % Calculated the number of quasi-particles in the film at each temperature.
    KID = QPnumber(n,KID,ChipInfo,designvalues);
    % Requires: KID(n).Temperature, KID(n).KIDnumber
    % Assigns: 
    %   KID(n).Nqp(:,1)
    %   KID(n).Delta
    %   KID(n).Volume
    %   KID(n).Tc
    %   KID(n).Area
    %   KID(n).Thickness
    %   KID(n).Fdesign
    %   KID(n).Qdesign
    %======================================================================
    
    if length(KID(n).Temperature) < 4 %Value 4 is arbitrary. Can be anything >= 2
        %If the number of temperatures is too low reponse measurements are useless
        
        %Must make some assignments to avoid problems later.
        KID(n).ddxdNqp = zeros(1,2);
        KID(n).dBdNqp = zeros(1,2);
        KID(n).ResponsivityM1 = zeros(length(KID(n).Temperature),2);
        KID(n).ResponsivityM2 = zeros(length(KID(n).Temperature),2);
        KID(n).Response = zeros(length(KID(n).Temperature),2);
        KID(n).ReImF0 = zeros(length(KID(n).Temperature),2);
    else
        %If the number of temperatures is high enough reponse measurements
        %are useful.
        % Calculate the temperature respons of the KID in various ways.
        KID = findresponse(n,KID,ChipInfo.FilmTc);
        % Requires: KID(n).Temperature, KID(n).Fres, KID(n).Nqp, KID(n).Ql,
        %           KID(n).Qi, KID(n).Qc,KID(n).S21data
        % Assigns:
        %   KID(n).ddxdNqp
        %   KID(n).dBdNqp
        %   KID(n).ResponsivityM1(:,1:2)
        %   KID(n).ReImF0(:,1:2)
        %   KID(n).Response(:,1:2)
        %   KID(n).ResponsivityM2(:,1:2)
    end
    
    %Plot the data to screen at least once. Such that it will be saved for
    %later use. Close it to avoid clutter.
    plotTdependence(n,KID,fitTbase{n},transformation{n},-1*KID(n).KIDnumber)
end
%==========================================================================
% Saving data externally
%==========================================================================
% Save the data to correct output files in either the S21 or Root Directory
KIDParameters = savethedata(KID,ChipInfo.path);
%Returns a summary file with most important parameters
%Returns the KIDParameters array, which contains the most important
%information about all N=#(KID,P) combinations at base temperature
%KIDParameters is a Nx16 array containing
%Dim 2:
%   1 = KIDID
%   2 = Tbase [K]
%   3 = Fres(Tbase) [GHz]
%   4 = Fdesign [GHz]
%   5 = Q
%   6 = Qi
%   7 = Qc
%   8 = Qdesign
%   9 = Pread [dBm]
%   10 = Pint [dBm]
%   11 = ddx/dNqp
%   12 = dR/dNqp'
%   13 = Tc [K]
%   14 = Delta [J]
%   15 = Film Thickness [um]
%   16 = Active Area [um^2]

%==========================================================================
% Plot the data to screen.
%==========================================================================
%If there are multiple powers per KID
if showallplots
    %User requests all plots are shown
    for n=1:NumberOfS21Files %Loop over all S21 files (all KIDs)
        plotTdependence(n,KID,fitTbase{n}(:,:),transformation{n}(:,:),n)
    end
elseif NumberOfS21Files == length(KIDnumbers)
    %No mixing with Power sweep.
    %Plot the results of all resonators to file and to screen to show user.
    for n=1:NumberOfS21Files %Loop over all S21 files (all KIDs)
        plotTdependence(n,KID,fitTbase{n}(:,:),transformation{n}(:,:),KID(n).KIDnumber)
    end
else
    %Mixing with Power sweep. A slightly different method will be used to
    %plot only the lowest power to screen for each resonator.
    for n=1:length(KIDnumbers)
        KIDid = KIDnumbers(n);
        KIDlocations = find(KIDParameters(:,1) == KIDid);
        [~,IndexPmin]=min(KIDParameters(KIDlocations,9));
        IndexPmin = KIDlocations(IndexPmin);
        plotTdependence(IndexPmin,KID,fitTbase{IndexPmin}(:,:),transformation{IndexPmin}(:,:),KID(n).KIDnumber)
    end
end

%==========================================================================
% Workspace saving and wrap-up
%==========================================================================
clear NumberOfS21Files RawS21files S21subdir ans fid info n showallplots

%Remove Path containing Subroutines from pathlist.
rmpath([pwd,filesep,'subroutines']);
%==========================================================================
%END OF HybridsS21_1RJ MAIN FUNCTION
%==========================================================================
end

function [data,Temperature,Power] = import_S21data(file,stopoffset)
%This function imports the S21 data from specified file. 
%From the main header only the power information is obtained.
%From the subheaders the temperature information is obtained.
%Always outputs the data, Temperatures and Power in order of unique
%monotonically increasing temperature
%
%INPUT:
% file == full file path
% stopoffset == number of temperatures at the end of the measurement disregarded
%OUTPUT:
% FOR N the number of temperatures, the following elements are read
% data == Nx1 cell array. Each cell contains a Mx3 array with [F (GHz), |S21| (dB), Phase (Rad)]
% Temperatures == Nx1 vector containing at the temperatures at which S21 is taken
% Power == Nx1 vector containing at the power (dBm) at which S21 is taken
%
%==========================================================================
% Example S21 data format inside the file.
%==========================================================================
%C:\KID Metingen\ADR\B6Ch7__7_11_06 16_02\S21\2D\KID41_73dBm_JB1.dat
%Power at KID:73dBm
%resonance Frequency in GHz :4.893858
%Q=40792.348421, Qc=56509.887020, Qi=146662.340671  
%
%
%Temperature in K:0.099953
%GHz	dB	Rad
%4.893115000	-6.722063065	2.685341760
%4.893118870	-6.721852779	2.685013925
%4.894663000	-6.540477276	2.794459459
%
%Temperature in K:0.105161
%GHz	dB	Rad
%4.893118740	-6.718620300	2.685038958
%4.893122571	-6.718233109	2.684799274
%==========================================================================

%Matlab textscan routine hard to understand (Reinier Janssen, 6-8-2012).

%==========================================================================
[tempdata,Temperature,Power,~] = import_data(file);

%Sort the temperatures and cut off the predefined number from the high T end.
%Making monotonic increasing Temperature in which all Temperature are
%unique. Non-unique values give problems later during interpolation steps.
[tempT,IncreasingTindex] = unique(Temperature(:,1),'first'); %'sorted' means picks the first of each encountered T duplicate.

%Number of temperatures corrected for the predefined high temperature cut-off
if stopoffset >= length(IncreasingTindex)
    %Would reduce the number of temperatures to zero or lower.
    %Give warning an ignore stopoffset
    fprintf('WARING [import_data]: Variable stopoffset is larger than the number of temperatures.\n')
    fprintf('Assuming stopoffset = 0 for this resonator.\n')
    Ntemperatures = length(IncreasingTindex);
else
    %Normal operation
    Ntemperatures = length(IncreasingTindex)-stopoffset; 
end

%Fill the output data in the correct way (monotically increasing T).
Temperature = tempT(1:Ntemperatures,1); %Temperature
data = tempdata(IncreasingTindex(1:Ntemperatures),1); %S21data
%General Matlab Note: tempdata and data are cell arrays usually to be
%addressed by {} brackets. However when copying cell to cell () brackets
%are to be used.
Power = -1*Power*ones(Ntemperatures,1); %Power

%Print information to screen
fprintf(['Selected temperatures range from ' num2str(Temperature(1,1)) ' to ' num2str(Temperature(end,1)) ' K\n']);

%==========================================================================
%END OF import_data SUBROUTINE
%==========================================================================
end

function [KID,transformation,fitTbase] = findparameters(n,KID)
%This function analyses the S21 data of all temperatures individually. For
%each temperature it does the following:
%1) Normalize |S21|(dB) to zero and convert it to magnitude space
%2) Fit the |S21|curve using one of the FitS21_x subroutine as selected by 
% FITF0 to obtain accurate values of Fres,Ql,Qc,Qi,|S21(Fres)|
%3) Calibrate the S21 resonance circle.

%INPUT:
% n == index in KIDstruct of KID considered.
% KID == large struct array containing all data.
%INPUT VIA GLOBAL:
% FITF0 == Switch index to determine the exact fitting routine used to 
% determine Fres and other parameters more accurately (To sub-frequency
% sampling resolution).
% The options for FITF0 are:
% 0 == use minimum value of measured |S21|. Q if fitted using FitS21_1
% 1 == use the original fit routine by Yates and Baselmans. (FitS21_1)
% 2 == use the function by Khalil et.al. to fit to |S21|(dB) (FitS21_2)
% 3 == [FUTURE]use the function by Khalil et.al. to fit to S21 (complex number, FitS21_3)
% 4 == [FUTURE] make a point-by-point fit in the complex plane. AND/OR make
%       use of theta(F) to determine bias point.
%
%==========================================================================
global FITF0 InterpolatPhaseref
Ntemperatures=length(KID(n).Temperature); %Find the number of temperatures.
S21data=KID(n).S21data; %Assign shorthand for S21data.

%Predefine a number of vectors with the number of temperatures
fres = zeros(Ntemperatures,1); % resonance frequency
Q = zeros(Ntemperatures,3); %[Ql,Qi,Qc]
S21min = zeros(Ntemperatures,1); %|S21(Fres)| (dB)
Paramerrors = zeros(Ntemperatures,5); %errors in [Fres,S21min(dB),Ql,Qi,Qc]

%==========================================================================
% Calibrate and convert |S21|
%==========================================================================
for p=1:Ntemperatures %Loop over all measured Temperatures
    %Smooth and normalize the |S21| data in dB space
    filtered=smooth(S21data{p}(:,2),3); %smoothing the |S21| data in dB space
    %normalise S21 in log space to the max(|S21|)
    S21data{p}(:,2)=S21data{p}(:,2)-max(filtered);
    
    %Convert to magnitude space
    S21data{p}(:,2)=10.^(S21data{p}(:,2)/20); %from log(magn) to magnitude
    %Division by 20 because dB is always defined in power. S21 is a
    %measurement in voltage. Since P \propto V^2: 
    %|S21| [dB] = 10*log10(|S21|^2) = 20*log10(|S21|)
    
    %Estimate the resonance parameters using fit and estimation routines.
    %Actual mode used is determined by FITF0 parameter.
    if S21data{p}(1,1) ~= 0 %catchong cols of 0's
        if p==1
            %At base temperature we want the evaluated fit function for future
            %plotting.
            [fres(p),Q(p,:),S21min(p),fitTbase,FitErrors] = FitS21main5(S21data{p}(:,1:3),FITF0);
        else
            [fres(p),Q(p,:),S21min(p),~,FitErrors] = FitS21main5(S21data{p}(:,1:3),FITF0);
        end
        Paramerrors(p,:) = FitErrors';
    else
        fres(p) = NaN;
        Q(p,:) = NaN;
        S21min(p) = NaN;
        disp(['problem in S21 datafile: T = ' num2str(KID(n).Temperature(p)) ]);
    end
end

%SAVE THE ESTIMATES IN THE KID STRUCT
KID(n).Ql=Q(:,1);
KID(n).Qi=Q(:,2);
KID(n).Qc=Q(:,3);
KID(n).Fres=fres;
KID(n).S21min=S21min;
KID(n).S21data = S21data;
KID(n).Paramerrors = Paramerrors;
%Calculate the Internal Power
KID(n).InternalPower = 10*log10((2/pi)*10.^(KID(n).ReadPower/10).*(KID(n).Ql.^2./KID(n).Qc));

%==========================================================================
% Rotate the phase of the resonance circle such that:
% * Fres(T) is located in the direction of (-0.5,0)
% * Circle at T0 is centered at (0,0)
% * Points far off-resonance are, for all temperatures, located at (0.5,0)
%==========================================================================
%Determine the phase at resonance for each KID.
phaseres = zeros(Ntemperatures,1); %Phase at resonance
for p=1:Ntemperatures
    if InterpolatPhaseref==1
        %Interpolate Phase data to get on resonance phase.
        phaseres(p) = interp1(S21data{p}(:,1),unwrap(S21data{p}(:,3)),fres(p),'linear');
    else
        %ALTERNATIVE (ORIGINAL)
        [~,F0index] = min(S21data{p}(:,2));
        phaseres(p) = S21data{p}(F0index,3);
    end
end
% Use the mean rotation to avoid problems for deepest KIDs
rotation = mean(phaseres);

%Use the fact the data is sorted that the lowest temperature is at index 1.
%Here the maximum value of Fres will be located. This is the resonance
%circle used as reference.

%rotating the circle
newphase=cell(Ntemperatures,1);
S21Real=cell(Ntemperatures,1);
S21Imag=cell(Ntemperatures,1);
OffsetReal=zeros(Ntemperatures,1);
for p=1:Ntemperatures
    newphase{p} = unwrap(S21data{p}(:,3)) - rotation; %Correct the phases of the S21 data
    S21Real{p} = S21data{p}(:,2).*cos(newphase{p}(:)); %Calculate Re{S21} with the calibrated phases
    S21Imag{p} = S21data{p}(:,2).*sin(newphase{p}(:)); %Calculate Im{S21} with the calibrated phases
    
     %Shift the S21 circle to center of the complex plane.
     OffsetReal(p) = (max(S21Real{p}(:))+10.^(S21min(1)/20))/2; 
     S21Real{p}(:)=S21Real{p}(:)-OffsetReal(p);
    
    %Update the S21data with calibrated S21 circle.
    S21data{p}(:,4)=S21Real{p}(:);
    S21data{p}(:,5)=S21Imag{p}(:);
end

%SAVE THE UPDATED S21data IN THE KID STRUCT
KID(n).S21data = S21data;

%Save the transformation for the use in plotting.
transformation(1)=rotation; %Mean phase correction [radians]
transformation(2)=mean(OffsetReal); %Mean offset correction.

%==========================================================================
%END OF findparameters SUBROUTINE
%==========================================================================
end

function [KID] = QPnumber(n,KID,ChipInfo,designvalues)
%This function determines the volume of the resonator from design values
%and then calculates the number of thermally excited quasi-particles in the
%active volume. It then stores all relevant information in the KID
%struct array. n is the index of the considered resonator in the KID struct
%array.

% Boltzmann constant
k=8.618*10^-5; %in eV/K
%==========================================================================
% Determine the Volume
%==========================================================================
KIDindexindesignmatrix =(designvalues(:,1)==KID(n).KIDnumber); %determine index of KID in designvalues matrix
KID(n).Area = designvalues(KIDindexindesignmatrix,4); %effective area=col 4; in um^2
KID(n).Thickness = ChipInfo.Thickness; %Thickness [um]
KID(n).Volume = ChipInfo.Thickness*KID(n).Area; %Determine volume in mu^3
%Display useful information to the screen as a check for the user.
disp(['KID ' num2str(KID(n).KIDnumber) ' has A= ' num2str(KID(n).Area) '[um^2] ',...
    'and V= ' num2str(KID(n).Volume) '[um^3] and is made of ' ChipInfo.material ]);
%==========================================================================
% Determine the material and determine Delta and N0
%==========================================================================
if strcmp(ChipInfo.material, 'Al')
    N0=1.70*10^10;%for AL, BEN Mazin p 19 in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*1.6;%From Ben Mazin
elseif strcmp(ChipInfo.material, 'Almod')
    N0=1.70*10^10;%for AL, BEN Mazin p 19 in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*1.81;%By Reinier
elseif strcmp(ChipInfo.material, 'AlBCS')
    N0=1.70*10^10;%for AL, BEN Mazin p 19 in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*1.76;%BCS value
elseif strcmp(ChipInfo.material, 'AlPascale')
    N0=1.70*10^10;%for AL, BEN Mazin p 19 in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*2.01;%Pascale value
elseif strcmp(ChipInfo.material, 'Ta')
    N0=4.08*10^10;% Ta
    delta=ChipInfo.FilmTc*k*1.8;
elseif strcmp(ChipInfo.material, 'Ti')
    N0=1.70*10^10;%Al, so WRONG!!!!
    delta=ChipInfo.FilmTc*k*1.76;%BCS
elseif strcmp(ChipInfo.material, 'TiN')
    N0=8.7*10^9;%for TiN, Leduc APL 2010 in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*1.81;%Escoffier PRL04
elseif strcmp(ChipInfo.material, 'NbTiN')
    N0=3.7*10^10;% Barends in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*1.76; %BCS (can be very wrong)
else
    % Default to BCS Aluminum
    fprintf('Warning: Unknown Material. Assuming BCS Aluminum\n')
    N0=1.70*10^10;%for AL, BEN Mazin p 19 in mu-3 eV-1
    delta=ChipInfo.FilmTc*k*1.76;%BCS value
end
%==========================================================================
% Calculate the number of quasi-particles and write all to data matrix
%==========================================================================
KID(n).Nqp = KID(n).Volume*...      %Quasiparticle Number
    2*N0*sqrt(2*pi*k*delta.*KID(n).Temperature).*exp(-(delta)./(k.*KID(n).Temperature));
KID(n).Delta = delta*1.602*10^(-19); %Gap, eV -> J.
KID(n).Tc = ChipInfo.FilmTc; %Measured Tc
KID(n).Fdesign = designvalues(KIDindexindesignmatrix,2);%design F0
KID(n).Qdesign = designvalues(KIDindexindesignmatrix,3);%design Qc
%==========================================================================
%END OF QPnumber SUBROUTINE
%==========================================================================
end

function  [KID] = findresponse(n,KID,Tcfilm)
%This function takes the KID struct and using the information inside
%determines:
%   Responsivity Method 1:
%   Using a fit of dx=(f-f0)/f0 or B=1/Qi vs Nqp Phase and
%   Amplitude Responsivity, respectively, are determined. The appropriate
%   conversion factors dR/dB and dTheta/ddx, which involve Q,Qi,Qc are used to
%   convert to the proper values. This induces the temperature dependence
%
%   Response:
%   The (change in) phase and radius at the base temperature resonance
%   frequency is determined.
%
%   Responsivity Method 2:
%   Using change in the phase and radius at the base temperature resonance
%   frequency a numerical difference method is used to determine the
%   responsivity for thermal quasiparticles.
%
%NOTE: Makes use of the fact the lowest (base) temperature is always at index 1.
%And that from there the temperature is monotonically increasing.
%==========================================================================
global TCfrac FITF0
Ntemperatures=size(KID(n).Temperature,1); %Number of temperatures / Length of most vectors

%Define Some Variables for easy use.
fres=KID(n).Fres; %resonant frequency
T=KID(n).Temperature; %Temperature
Nqp=KID(n).Nqp; %Quasiparticle Number
dx=KID(n).Fres/KID(n).Fres(1)-1; %dx = (F0(T)-F0(Tbase))/F0(Tbase)
Q=KID(n).Ql; %Q
Qc=KID(n).Qc; %Qc
Qi=KID(n).Qi; %Qi
S21data = KID(n).S21data; %S21data

%==========================================================================
% Method 1: Use a fit to obtain responsivities.
% Phase Responsivity: Fit the slope of dx vs Nqp (Most reliable method found)
% Amplitude Responsivity: Fit the slope of 1/Qi vs Nqp (Nice estimate, not extremely reliable)
%==========================================================================
%To ignore TLS effects only temperatures above Tc*TCfrac are selected for
%the determination of the phase response. Default TCfrac = 1/8;
HighT=find(T>(Tcfilm*TCfrac)); %Use find instead of logicals due to convenience later.
if length(HighT) <= 3 
    %If there are to little temperatures above the specified fraction of Tc.
    %Give a warning.
    fprintf(['WARNING [responsivity determination]: There are ',num2str(length(HighT)),...
        ' Temperatures above ',num2str(TCfrac),'Tc.\n'])
    fprintf('A fit over all temperatures is done to get the responsivity.\n')
    %Fit over all temperatures instead.
    HighT = find(T);
end
fprintf(['Responsivity fitted over ',num2str(T(HighT(1))),' till ',num2str(T(HighT(end))),' [K].\n'])
fprintf(['Requested lower temperature limit: ',num2str(Tcfilm*TCfrac),' [K].\n'])

%Linear fit over the dx(Nqp) data.
%Scale everything to be less subject to internal fitting limits (all values closer to 1).
phasefitparam = polyfit(Nqp(HighT)/1e6,dx(HighT)*1e6,1);
%Convert to Phase Responsivity: dtheta/dnqp=dtheta/ddx*ddx/dnqp=4Q*ddx/dnqp
PhaseResponsivity1 = -4*Q*phasefitparam(1)/1e12;

%Linear fit over the 1/Qi(Nqp) data.
ampfitparam=polyfit(Nqp(HighT),1./Qi(HighT),1);
%Calculate the dR/d(1/Qi)
%This is the numerical variation of radius with 1/Qi (normalized to base T)
drdb = 2*Q;
%Convert to Radius Responsivity: d(1/Qi)dNqp * dRd(1/Qi)
RadiusResponsivity1 = drdb.*ampfitparam(1);

%Save all interesting parameters to the KID struct.
KID(n).ddxdNqp = phasefitparam.*[1e-12 1e-6]; %FIT parameters of the freq. responsivity.
KID(n).dBdNqp = ampfitparam; %FIT parameters for the 1/Qi responsivity
KID(n).ResponsivityM1(:,1) = PhaseResponsivity1; %Phase responsivity from FIT
KID(n).ResponsivityM1(:,2) = RadiusResponsivity1; %Radius responsivity from FIT

%==========================================================================
% Phase Response && Radius Response
% As a function of temperature determine the value of R and Theta at
% frequency = f0(Tbase)
%==========================================================================
Response = zeros(Ntemperatures,8);
%dim 2:
%   1 == Re{S21(T)} @f0(T0)
%   2 == Im{S21(T)} @f0(T0)
%   3 == R{S21(T)} @f0(T0)
%   4 == Theta{S21(T)} @f0(T0)
%   5 == R{S21(T)} - R{S21(T0)} @f0(T0)
%   6 == Theta{S21(T)} - Theta{S21(T)} @f0(T0)
%   7 == Smoothed version of (:,3)
%   8 == Smoothed version of (:,4)

for p=1:Ntemperatures
    %Use linear interpolation between the data points to obtain the value at f0(T0)
    Response(p,1) = interp1(S21data{p}(:,1),S21data{p}(:,4),fres(1),'linear'); %Real in calibrated circle
    Response(p,2) = interp1(S21data{p}(:,1),S21data{p}(:,5),fres(1),'linear'); %Imag in calibrated circle
end
Response(:,3) = abs(Response(:,1)+1i*Response(:,2)); %Radius
Response(:,4) = atan2(Response(:,2),Response(:,1)); %Theta [rad]
%Correction to remove 2pi phase jumps
Response(:,4) = unwrap(Response(:,4));
%OLD METHOD of removing jumps
%PhaseJumps = Response(:,4) <-0.75*pi;
%Response(PhaseJumps,4) = Response(PhaseJumps,4)+2*pi;

%Calculate the response to temperature (change in R and Theta)
Response(:,5) = Response(:,3) - Response(1,3); %R{S21(T)} - R{S21(T0)} @f0(T0)
Response(:,6) = Response(:,4) - Response(1,4); %Theta{S21(T)} - Theta{S21(T)} @f0(T0)

%Save the interesting values to the KID struct
KID(n).ReImF0 = Response(:,1:2);
KID(n).Response = Response(:,5:6);

%==========================================================================
% Method 2: Direct numerical differentiation (after smoothing) to obtain responsivity.
% Phase Responsivity: [Theta(Ti+1) - Theta(Ti)]/[Nqp(Ti+1) - Nqp(Ti)]
%   (Suffers from backbending due to TLS and the fact that at high T you
%   reach the other side of the circle)
% Amplitude Responsivity: [R(Ti+1) - R(Ti)]/[Nqp(Ti+1) - Nqp(Ti)]
%==========================================================================

%Smooth the responses for a nicer responsivity
Response(:,7) = smooth(Response(:,3),3); %Smooth R
Response(:,8) = smooth(Response(:,4),3); %Smooth Theta

%Determine the responsivity by numerical differentiation.
ResponsivityM2 = zeros(Ntemperatures,2);
ResponsivityM2(2:end,1) = (Response(2:end,8)-Response(1:end-1,8))./(Nqp(2:end)-Nqp(1:end-1));%Phase Responsivity.
ResponsivityM2(2:end,2) = (Response(2:end,7)-Response(1:end-1,7))./(Nqp(2:end)-Nqp(1:end-1));%Radius Responsivity.

%Save the method 2 responsivity to KID struct
KID(n).ResponsivityM2 = ResponsivityM2;

%==========================================================================
%END OF findresponse SUBROUTINE
%==========================================================================
end

function KIDParameters = savethedata(KID,MeasurementPath)
%This function takes all the analysis data of KID "n" stored in the
%structarray "KID" and writes a number of csv files to the path containing
%the original S21data. Here it writes two file:
%   S21lowT: S21 curve at the base temperature
%   Tdep: The temperature dependence of all parameters.
%It also writes and overview file with the parameters of Tdep at base
%temperature to compare all KIDs.
%TODO: In addition it writes an file with a summary of the most important
%parameters to the root "MeasurementPath".
%TODO: Declutter the Tdep file (probably make it like Summary), but due to
%backwards compatibility and its use in KIDNEPv1 this is kept for now.
global S21path

%Predefine some headers
LowTheader={'Freq','S21dB','theta','Re Norm', 'Im Norm'}';
DepTheader={'KIDID','T','Q','Qc','Qi','F0','S21@res','power','Pinternal',...
    'response','dtheta/dN FIT','dtheta/dN numerical','qpnumber','dx/dN','V[um3]','Delta',...
    'Re response','Im response','dR/dN','dRdB','grad 1/Qi vs nqp',...
    'Tc','designF0','designQ','KID A','film thickness'}';
SummaryHeader = {'KIDID','T0','F0','Fdesign','Q','Qi','Qc','Qdesign',...
    'Pread [dBm]','Pint [dBm]','ddx/dNqp','dR/dNqp',...
    'Tc [K]','Delta [J]','Film Thickness [um]','Active Area [um^2]'}';

%Initialize the data array of the old overview file.
AllKIDlowToverview = zeros(length(KID),length(DepTheader));
for n=1:length(KID)
    %Get the S21path location and full KID substring from original file name.
    Filenames = KID(n).filename;
    VariableInfo = whos('Filenames');
    if strcmp(VariableInfo.class,'cell') %If cell vector (a file per temperature)
        %KID(n).filename
    else %A string, only a single file per T sweep
        %Determine the root name of this KID.
        IndexNameEnd = strfind(KID(n).filename,'.dat');
        RootName = KID(n).filename(1:IndexNameEnd-1);
    end
    
    %Lowest Temperature data is always located at index 1 due to sort in the
    %import_data. So the full S21 transmission is given by
    KIDS21lowT = KID(n).S21data{1,1}(:,:);
    KIDS21lowT(:,2) = 20*log10(KIDS21lowT(:,2)); %convert |S21| to dB
    
    %Write the base temperature S21 to csv file
    resultsfile = [RootName,'S21lowT.csv'];
    
    %Construct a T dependent parameter array
    TdepParameters = zeros(length(KID(n).Temperature),length(DepTheader));
    TdepParameters(:,1)=KID(n).KIDnumber;
    TdepParameters(:,2)=KID(n).Temperature;
    TdepParameters(:,3)=KID(n).Ql;
    TdepParameters(:,4)=KID(n).Qc;
    TdepParameters(:,5)=KID(n).Qi;
    TdepParameters(:,6)=KID(n).Fres*1e9; %[GHz -> Hz for backwards compatibility]
    TdepParameters(:,7)=KID(n).S21min; %[dB]
    TdepParameters(:,8)=KID(n).ReadPower; %[dBm]
    TdepParameters(:,9)=KID(n).InternalPower; %[dBm]
    TdepParameters(:,10)=KID(n).Response(:,2); %Phase Response
    TdepParameters(:,11)=KID(n).ResponsivityM1(:,1); %Phase Responsivity [FIT]
    TdepParameters(:,12)=KID(n).ResponsivityM2(:,1); %Phase Responsivity [Numerical]
    TdepParameters(:,13)=KID(n).Nqp;
    TdepParameters(:,14)=KID(n).ddxdNqp(1); %Frequency Responsivity [FIT slope]
    TdepParameters(:,15)=KID(n).Volume;
    TdepParameters(:,16)=KID(n).Delta; %[J]
    TdepParameters(:,17)=KID(n).ReImF0(:,1); %Re{S21(T)} @F0(T0)
    TdepParameters(:,18)=KID(n).ReImF0(:,2); %Im{S21(T)} @F0(T0)
    TdepParameters(:,19)=KID(n).ResponsivityM1(:,2); %Amplitude Responsivity [FIT]
    TdepParameters(:,20)=... %dRd(1/Qi) [Modulation of the fit parameter]
        2*KID(n).Ql;
    TdepParameters(:,21)=KID(n).dBdNqp(1); %1/Qi responsivity [FIT slope]
    TdepParameters(:,22)=KID(n).Tc;
    TdepParameters(:,23)=KID(n).Fdesign;
    TdepParameters(:,24)=KID(n).Qdesign;
    TdepParameters(:,25)=KID(n).Area;
    TdepParameters(:,26)=KID(n).Thickness;
    
    %Write the Temperature dependent parameters of the sweep to file
    resultsfile = [RootName,'Tdep.csv'];
    WriteSRONcsv(resultsfile,TdepParameters,DepTheader,'%.12e');
    
    %Take the parameters at the lowest temperature for an overview file.
    AllKIDlowToverview(n,:)=TdepParameters(1,:);
end

% %Identify the indexes of interesting parameters in TdepParameters.
ParamIndexes = [1,2,6,23,... %KIDID, T0, Fres, Fdesign
    3,5,4,24,... %Q's
    8,9,14,19,... %Powers and responsivities
    22,16,26,25]; %Review utilities
KIDParameters = AllKIDlowToverview(:,ParamIndexes);

%==========================================================================
%END OF savethedata SUBROUTINE
%==========================================================================
end
        
function plotTdependence(n,KID,fitTbase,transformation,figurehandle)
% This function takes KID struct array and plots all the available
% information inside. This should give insight in the analysis done by
% other subroutines in the HybridS21 analysis program.
%
% INPUT
% n == number of the KID in the KID struct
% KID == KID struct containing all the data
% fitTbase == Nx2 array [F,|S21|(dB)] that is the evaluated function of the
%             fit result given by the FitS21_X routine.
% transformation == mean (in T) transformation applied to calibrate the
%                   resonance circle
% figurehandle == absolute value of this number will become the number of
%                   the figure. If negative, the figure will be closed
%                   after saving. 0 is not allowed.
%
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%modified slightly to show three methods of responsivity 20-2-2007
%Modified aug 28 2006 JochemB; minor cosmetic changes, 
%addition of saturation point in 3.3.6; resp. points in 3.3.3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set format for printing to screen
format('long','e');
%Generate a color rainbow for temperature sweeps
Ntemperatures = length(KID(n).Temperature);
[Tcolors,~] = GenerateColorsFromMap(Ntemperatures,'RainbowReinier');

%==========================================================================
% Figure == KIDnumber: Overview of all temperature dependencies of the KID
%==========================================================================
figure('Name',num2str(abs(figurehandle)),'visible','off');
clf %clear current figure

%SUBPLOT 3.3.1 KID RESONANCE at base temperature
subplot(3,3,1)
hold on
%Plot the data, fit and resonance frequency
plot(KID(n).S21data{1}(:,1),KID(n).S21data{1}(:,2),'bx','MarkerSize',5) %The measured resonance
plot(fitTbase(:,1),10.^(fitTbase(:,2)/20),'k-','LineWidth',1) %Plot the base Temperature fit
plot(KID(n).Fres(1),10^(KID(n).S21min(1)/20),'kd',...
    'MarkerFaceColor','r','MarkerSize',6)%The determined Fres(0) and |S21(Fres(0))|
%Make a nice layout.
xlim([min(KID(n).S21data{1}(:,1)),max(KID(n).S21data{1}(:,1))])
ylim([0,1.1]) %Due to being in magnitude 0<|S21|<= 1
legend('hide') %Otherwise an automatic legend pops up from the plot fitresult combination.
xlabel('F [GHz]')
ylabel('|S21| [mag]')
title(['KID ',num2str(KID(n).KIDnumber),' @ ',num2str(KID(n).ReadPower(1)),' [dBm]'])
hold off

%SUBPLOT 3.3.2 KID RESONANCES as a function of temperature
subplot(3,3,2)
hold on
for p=1:Ntemperatures
    plot(KID(n).S21data{p}(:,1),20*log10(KID(n).S21data{p}(:,2)),'.-',...
        'color',Tcolors(p,:),'MarkerSize',2) %Plot the resonances in dB space
end
%Two stage loop so that resonance freq dots are on top.
for p=1:Ntemperatures
    plot(KID(n).Fres(p),KID(n).S21min(p),'ko','MarkerSize',4,'MarkerFaceColor','k') %Plot Fres
end
%Make a nice layout
xlabel('F [GHz]')
ylabel('S21 [dB]')
title('KID resonances')
axis tight;
hold off

%SUBPLOT 3.3.3 Q FACTORS
subplot(3,3,3)
semilogy(KID(n).Temperature,KID(n).Ql,'ko-',... %Loaded Q
    'MarkerSize',4,'LineWidth',1)
hold on %Hold on after first semilogy to keep y-axis in log scale
semilogy(KID(n).Temperature,KID(n).Qc,'ro-',... %Coupling Q
    'MarkerSize',4,'LineWidth',1)
semilogy(KID(n).Temperature,KID(n).Qi,'bo-',... %Internal Q
    'MarkerSize',4,'LineWidth',1)
%Make a nice layout
legend('Q_l','Q_c','Q_i')
xlabel('T [K]')
ylabel('Q')
title('Q factors')
hold off

%SUBPLOT 3.3.4 KID RESONANCES in complex plane (As Measured)
subplot(3,3,4)
hold on
for p=1:Ntemperatures
    %Plot the resonance as measured
    MeasuredRe = KID(n).S21data{p}(:,2) .* cos(KID(n).S21data{p}(:,3));
    MeasuredIm = KID(n).S21data{p}(:,2) .* sin(KID(n).S21data{p}(:,3));
    plot(MeasuredRe,MeasuredIm,'o-',...
        'color',Tcolors(p,:),'MarkerSize',2);
end
%Two stage loop to have points on top.
for p=1:Ntemperatures
    %Recalculate the resonance as measured
    MeasuredRe = KID(n).S21data{p}(:,2) .* cos(KID(n).S21data{p}(:,3));
    MeasuredIm = KID(n).S21data{p}(:,2) .* sin(KID(n).S21data{p}(:,3));
    %Determine the measurement point closest to resonance
    [~,F0index] = min(abs(KID(n).Fres(p) - KID(n).S21data{p}(:,1))); 
    plot(MeasuredRe(F0index),MeasuredIm(F0index),'ko',...
        'MarkerSize',5,'MarkerFaceColor','k')
    %Apply the inverse transformation to the (Re,Im)@F0(T0) value.
    RealF0 = KID(n).ReImF0(p,1)+transformation(2); %Reverse transformation on Real axies
    Phase = atan2(KID(n).ReImF0(p,2),RealF0)+transformation(1); %Calculate Phase & Reverse transformation theta
    Radius = abs(RealF0+i*KID(n).ReImF0(p,2)); %Calculate new radius.
    %Plot the trace of F0(T0) in the complex plane
    plot(Radius*cos(Phase),Radius*sin(Phase),'kd',...
        'MarkerSize',5,'MarkerFaceColor','r')
end
%Make a nice layout
grid on
xlabel('Re')
ylabel('Im')
title('Measured KID resonances in complex plane')
axis([-1 1 -1 1])
hold off

%SUBPLOT 3.3.5 KID RESONANCES in complex plane (Calibrated)
subplot(3,3,5)
hold on
for p=1:Ntemperatures
    %Plot the calibrated resonance circles
    plot(KID(n).S21data{p}(:,4),KID(n).S21data{p}(:,5),'.-',...
        'color',Tcolors(p,:),'MarkerSize',1); %Calibrated Measurent in Circle
end
for p=1:Ntemperatures
    %Determine the measurement point closest to resonance
    [~,F0index] = min(abs(KID(n).Fres(p) - KID(n).S21data{p}(:,1))); 
    plot(KID(n).S21data{p}(F0index,4),KID(n).S21data{p}(F0index,5),'ko','MarkerSize',5,'MarkerFaceColor','k'); %Location Fres(T)
    %Plot the F0(T0) trace in the complex plane
    plot(KID(n).ReImF0(p,1),KID(n).ReImF0(p,2),'kd','MarkerSize',5,'MarkerFaceColor','r'); %Location Fres(Tbase)
end
%Make a nice layout
grid on
xlabel('Re')
ylabel('Im')
title('Calibrated KID resonances in complex plane')
axis([-1 1 -1 1]*0.5)
hold off

%SUBPLOT 3.3.6 KID RESPONSIVITY

%Predetermine some axis limits (incl errorcatching for max(Responsivity) <= 0)
Ylims = zeros(2,2); %dim 1=[min,max], dim 2=[Phase, Amp]
%PHASE
Ylims(2,1) = ceil(max(KID(n).ResponsivityM1(:,1)*1e6));
if Ylims(2,1) < 0
    Ylims(1,1) = Ylims(2,1);
    Ylims(2,1) = -0.2*Ylims(1,1);
elseif Ylims(2,1) == 0
    Ylims(2,1) = 1;
    Ylims(1,1) = -0.2;
else
    Ylims(1,1) = -0.2*Ylims(2,1);
end
%AMPLITUDE
Ylims(2,2) = ceil(max(KID(n).ResponsivityM1(:,2)*1e6));
if Ylims(2,2) < 0
    Ylims(1,2) = Ylims(2,2);
    Ylims(2,2) = -0.2*Ylims(1,2);
elseif Ylims(2,2) == 0
    Ylims(2,2) = 1;
    Ylims(1,2) = -0.2;
else
    Ylims(1,2) = -0.2*Ylims(2,2);
end

%Make the actual plot
AX9_1 = subplot(3,3,6);
hold on
%Plot the first axis containing the phase response
set(AX9_1,'Box','off','YColor','b','Ylim',Ylims(:,1))
ylabel(AX9_1,'\delta \theta / \delta N_{qp} [rad/10^6 qp]')
%Plot in these axes the Phase Responsivity
PhaseM1 = line(KID(n).Temperature,KID(n).ResponsivityM1(:,1)*1e6,'Parent',AX9_1); %Phase responsivity by FIT
PhaseM2 = line(KID(n).Temperature,-1*KID(n).ResponsivityM2(:,1)*1e6,'Parent',AX9_1); %Phase responsivity by numerical differentiation
set(PhaseM1,'LineStyle','-','Marker','o','Color','b')
set(PhaseM2,'LineStyle','--','Marker','.','Color','b')
%Create the second axis and put in the amplitude response
AX9_2 = axes('Position',get(AX9_1,'Position'),...
    'color','none','YColor','r',...
    'YAxisLocation','right',...
    'XTick',[],...
    'Ylim',Ylims(:,2),...
    'Box','off');
ylabel(AX9_2,'\delta R / \delta N_{qp} [1/10^6 qp]')
%Plot in these axes the Amplitude Responsivity
AmpM1= line(KID(n).Temperature,KID(n).ResponsivityM1(:,2)*1e6,'Parent',AX9_2); %Amplitude Responsivity by FIT
AmpM2 = line(KID(n).Temperature,-1*KID(n).ResponsivityM2(:,2)*1e6,'Parent',AX9_2); %Amplitude Responsivity by numerical differentiation
set(AmpM1,'LineStyle','-','Marker','o','Color','r')
set(AmpM2,'LineStyle','--','Marker','.','Color','r')
%General Stuff
xlabel(AX9_1,'T [K]')
title('KID responsivity')
legend([PhaseM1,PhaseM2,AmpM1,AmpM2],'Phase (FIT)','Phase (NUM)','Amp (FIT)','Amp (NUM)','location','SouthWest')
hold off

%SUBPLOT 3.3.7 KID 1/Qi fit vs nqp
subplot(3,3,7)
hold on
plot(KID(n).Nqp,1./KID(n).Qi,'k.-',... %1/Qi as function of Nqp
    'MarkerSize',4,'LineWidth',1)
plot(KID(n).Nqp,polyval(KID(n).dBdNqp,KID(n).Nqp),'r-',... %Fit to above data
    'MarkerSize',4,'LineWidth',1)
%Make a nice layout
xlabel('N_{qp}')
ylabel('1/Q_i')
title('Q_i factor fit');

%SUBPLOT 3.3.8 Frequency versus Nqp
subplot(3,3,8)
hold on
DX = (KID(n).Fres - KID(n).Fres(1))/KID(n).Fres(1);
plot(KID(n).Nqp,DX,'k.-',... %(F-F0)/F0 as function of Nqp
    'MarkerSize',4,'LineWidth',1)
plot(KID(n).Nqp,polyval(KID(n).ddxdNqp,KID(n).Nqp),'r-',... %Fit to above data
    'MarkerSize',4,'LineWidth',1)
%Make a nice layout
xlabel('N_{qp}')
ylabel('dx=(F-F0)/F0')
title('Frequency vs Nqp');

%SUBPLOT 3.3.9 KID RESPONSE
AX9_1 = subplot(3,3,9);
hold on
%Plot the first axis containing the phase response
set(AX9_1,'Box','off','YColor','b')
ylabel(AX9_1,'\delta \theta [degree]')
%Plot in these axes the Phase Response
PhaseLine = line(KID(n).Temperature,KID(n).Response(:,2)*180/pi,'Parent',AX9_1); %Phase response
set(PhaseLine,'LineStyle','-','Marker','.','Color','b')
%Create the second axis and put in the amplitude response
AX9_2 = axes('Position',get(AX9_1,'Position'),...
    'color','none',...
    'YColor','r',...
    'YAxisLocation','right',...
    'XTick',[],...
    'Box','off');
ylabel(AX9_2,'\delta R')
%Plot in these axes the Amplitude Response
PhaseLine = line(KID(n).Temperature,KID(n).Response(:,1),'Parent',AX9_2); %Amplitude response
set(PhaseLine,'LineStyle','-','Marker','.','Color','r')
%General Stuff
xlabel(AX9_1,'T [K]')
title('KID response')
hold off


%SAVE the figure
Figfile=[KID(n).filename(1:end-3)];
MakeGoodFigure(15,17,13,Figfile,1);

%Clean up if desired.
if figurehandle > 0
    %keep plot on the screen
else
    %close the figure
    close(num2str(abs(figurehandle)));
end

%==========================================================================
%END OF plotTdependence SUBROUTINE
%==========================================================================
end