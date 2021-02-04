function [TotalPbb,Filtertransmission,NEP,method]=blackbody_14(Tbb,method,plotdata,pad,pOs)
% Update 15-12: Error in Bunching term corrected (was ok up to V8))
% updated 6-12-2016 for Deshima filters and added optional pad variable. It
% is then 1:1 the same as the function used in the MUX.
% NEW VERSION 12-12-2015 Jochem. REMOVE STUPID REINIER COLORS
% gets the filter transmisison (from datafiles) and uses this + the throughput, to
% calculate the power as fy of BB Temperature. The sctipt uses
% numerical integartion from the gap frequency to the
% max freq set in the function. Filter transmission is set to 0 outside
% of this range. At the min freq a realistic datapoint is added
% artificially, to get a reasonable estimate *using interplation) of the
% lekage at low F.
%
%INPUT:
% Tbb = A 1xN vector of temperatures of the blackbody.
% pad = path where the filterfiles are, can be left out, pad = pwd in that case. it will look in [pad,filesep,'filterfiles',filesep];%
% pOs : set to 0 to diable any comments on the screen 
% method. A struct containing information about the used setup including:
%         filters, lens size, and aperture determination. The following
%         struct elements are required:
%
%   method.pol = Number of polarizations used (equal to either 1 or 2)
%
%   method.filter='1_6THz'              % 1.6 THz SPICA SAFARI stack
%   method.filter='350 GHz 4Filters'    %350 Cardiff bandpass with additional LPF
%   method.filter='350 GHz Deshima'     %350 GHz LPF (2x) on cold box, and K1817BPF and K1785 LPG on chip and 1 Thz LPF on BB%
%   method.filter = '650 GHz'           %650 GHz BPF SH (B768), LTbox BPF (B768) and 1 Thz LPF on BB%%
%   method.filter='850 GHz'             %850 GHz filterstack used first with LT010
%
%   method.tp='lambda^2'                %using lambda^2 throughput over the entire filter band%
%   method.tp='Geometrical'             %Using geometrical calculation for throughput determination.
%   
%   ONLY FOR method.tp='Geometrical'%
%   method.lensdiameter=1 %lens diameter in millimeters. 
%   method.opening_angle=10 %angle wrt optical propagation (i.e. total angle is 2x larger) in degrees%
%
%   ONLY FOR method.tp='lambda^2'%
%   method.eta_c = 1.0              %[Optional, Default = 1.00] total SO couplig (CST), only for lambda^2%
%       
%   OPTIONAL    
%   method.GR = 4;                  % prefactor in the GR noise term of the NEP. default = 4
%                                   % following Flanigan et al (APL 108, 083502 (2016). %
%                                   % We used GR=2 in the past (in Pieter/Reineir/SY work%
%   method.freqresolution = 2       %[Optional, Default = 2] frequency resolution in GHz used
%                                       to integrate the blackbody spectrum with all filters.
%                                       NOTE: for the 350 GHz and 325 GHz filters 2 GHz is
%                                       adviced as minimum.
%   method.MaxFreq = 5e12           %[Optional, Default= 5e12] max Freq. of integration%
%   method.Delta = 45.6             %[Optional, Default = 45.6] value of Delta (half
%                                       the superconducting gap) in GHz
%   method.eta_pb = 0.57            %[Optional, Default = 0.57] value of the photon
%                                       pair breaking efficiency. Default from Kozorezov
%                                       et.al. PRB 2000
%
%  plotdata. (either 1 [true] or 0 [false]) If true, the routine plots the
%            photon noise NEP's of Al resonators. 
% 
%OUTPUT: original nput Tbb = 1xN vector
% TotalPbb.             % 1xN vector of blackbody power received by the lens/antenna ateach blackbody temperature%
% method.Etendue.       % Optical throughput of the system at the center frequency of the band %
% method.centrefreq     % central frequency
% method.filterBW       % effective BW
% Filtertransmission.   %2xM vector of total filter transmission [frequency (Hz),transmission].
%                       frequency is given with the resolution specified in method.freqresolution  
% NEP.                  %struct containing four 1xN vectors that for each blackbody temperature give the expected NEP due to:%
%   NEP.g_r         %1xN vector Generation Recombination Noise (using method.Delta) 
%   NEP.poisson     %1xN vector Poissonian photon noise
%   NEP.wave        %1xN vector Wave Bunching Photon Noise
%   NEP.totphoton   %1xN vector Total NEP due to Photon, GR and wavebunching noise.
%
%SUBROUTINES:
% getfilterform (included below)
% plotresults (included below)
%
%REQUIRED FILES (filters, located in /filterfiles subdir):
% W969 37cmLPESCUBAII.txt
% W1275_350GHz.dat
% totalLbandfilters.txt
% B386 18cm LPE.dat
% W1052 14cm LPE SCUBAII.dat
% H20coupling efficiency upto1_3THz.txt
% H10coupling efficiency.txt
%
%
%VERSION: 12.0
%   V1.0 original blackbody script by Stephen Yates
%   V2.0 overall clarity update by Jochem Baselmans
%   V2.1 update of the program documentation / comments
%   V3.0 further comments update. Change to frequency space as primary
%   variable instead of 1/lambda space for integration. Note: most was
%   almost always the case, but specified awkwardly.
%   V3.1 (2012-09-21,JB/PdV): Added the SPICA Lband filter at 1.6 THz and
%                             1.6 THz full ADR aperture
%   V4.0 (2012-11-09,RJ): Made references to filterfiles in seperate
%   subdir.
%   V5.0 (2013-03-11,RJ): If a beam pattern calculated entendue is used for
%   for the radiation efficiency, then the coupling efficiency is used as
%   an additional filter. This means the total antenna efficiency is
%   included for a co-aligned lens-antenna and pinhole.
%   V6.0 (2013-05-23,RJ): Included the new filterstack for the 350 GHz
%   bandpass experiment.
%   V6.1 (2013-05-23,RJ): Included some additional figures in the plotting
%   subroutine.
%   V6.2 (2013-06-05,RJ): Updated the beam pattern etendues based on
%   refined mesh integrations and up-to-date pinhole heights.
%   V8.0 (2013-12-16,RJ): Included the 850 GHz filterstack
%   V8.1 (2014-02-20,RJ): Included lambda^2 TP for 850 GHz
% V12 simplified to only lambda^2 and geometrical, all filter/lens issues
% are now included in the \eta_c parameter
%DATE: Dec. 12 2015
%AUTHOR: Reinier Janssen/Jochem Baslemans

%==========================================================================
% Set overal values and correct for missing optional inputs.
%==========================================================================
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
if nargin == 3
    pad = pwd;
    pOs = 1;
end

if nargin == 4
    pOs = 1;
end

% Define some global parameters and constants of nature
global c;
h = 6.6262e-34;		% J.s
c = 2.9998e8;		% m/s
k = 1.3806e-23;		% J/K
%==========================================================================
% Set the parameters for the numerical integration. 
%==========================================================================


if isfield(method,'pol') == 0 || isfield(method,'filter') == 0 || isfield(method,'tp') == 0 
    error('method struct not complete');
end

% Check if the polarization was specified correctly
if method.pol~=1 && method.pol~=2
    error('ERROR: Polarization should be either 1 or 2')
end

%check lambda^2 case
if strcmp(method.tp,'lambda^2')
    method.lensdiameter=[];method.opening_angle=[];
    if isfield(method,'eta_c') == 0
        error('No coupling efficiency given.')
    end
end
%check Geometrical case
if strcmp(method.tp,'Geometrical')
    if isfield(method,'lensdiameter') == 0 || isfield(method,'opening_angle') == 0
        error('Geometrical throughput chosen, but no .lensdiameter or .opening_angle defined');
    end
    method.eta_c=1;
    method.solidangle=1*pi*(1-cosd(method.opening_angle)^2); %including  effective area reduction at large angles
end
%==========================================================================
% Check if the optional parameters are given as input.
%==========================================================================
if isfield(method,'MaxFreq') == 0
    method.MaxFreq = 5e12;
end
if isfield(method,'GR') == 0
    method.GR = 4;
end
if isfield(method,'Delta') == 0
    method.Delta = 45.6;     % gap estimate Aluminium in GHz, from H10
end
delta = method.Delta*1e9*h; % gap value in J

% Pair breaking efficiency
if isfield(method,'eta_pb') == 0
   method.eta_pb = 0.57;   % 0.57 pair breaking efficiency by Kozorezov et.al. RPB 2000
   %went to the value of 0.4 using the results from Guruswamy et al. (2015)
end

% Frequency Resolution used for integration
if isfield(method,'freqresolution') == 0
    method.freqresolution = 2;     % Frequency resolution default of 2 GHz
end

%correct method for throughput
if isfield(method,'tp') == 0 || ~(strcmp(method.tp,'lambda^2') || strcmp(method.tp,'Geometrical'))
    error('No thropughput defined')     % Perfect alignment
end

% Show the specified method as a double check for the user.
if pOs ~=0
    fprintf('Specified Method for Blackbody Power calculation:\n')
    disp(method);
end
%==========================================================================
% Start script
%==========================================================================

% Set the parameters for the numerical integration. 
minfrequency = 2*method.Delta*1e9; %minimum frequency [Hz]: gap frequency in Hz.
freq = minfrequency:(method.freqresolution*1e9):method.MaxFreq;
nowave = length(freq);

% Get the filters and interpolate of obtain the values at the desired frequencies.%
[filterform, legendstrr]=getfilterform(method,minfrequency,pad,pOs);% filterform is a cell array with the raw filter specs

%Interpolate the data from file to obtain transmission at the desired integration frequencies.%
filter = ones(1,nowave);        %start the overall filter with full transmission
for I=1:nowave                  %loop over each frequency
    for n=1:length(filterform)  %loop over all consecutive filters
        if freq(I)/c > max(filterform{n}(:,1)) || freq(I)/c < min(filterform{n}(:,1))
            filter(I) = 0;% Outside the filter range. No transmission is assumed.
        else
            % Inside the filter range. Use interpolation to obtain value at
            % specified frequency. Multiply the current transmission with
            % the interpolation value.
            filter(I)=filter(I)*interp1(filterform{n}(:,1),filterform{n}(:,2),freq(I)/c,'pchip','extrap');
        end
    end
end
%Prepare the filter transmission for output.
filter((filter<0))=0;
Filtertransmission = [freq',filter']';
method.centrefreq=sum(freq.*filter)/sum(filter);                    % central frequency
method.filterBW=sum(filter*(method.freqresolution*1e9))/max(filter); % effective BW
%==========================================================================
% Numerical Integration of the Planck function over all frequencies for all
% temperatures requested.
%==========================================================================
% initialize arrays
N_Tbb=length(Tbb); %number of blackbody temperatures
irradiation=cell(1,length(Tbb));powerbla2=cell(1,length(Tbb));
TotalPbb=zeros(1,N_Tbb); %Total power arriving at lens/antenna
NEP.g_r=zeros(1,N_Tbb); %GR noise
NEP.poisson=zeros(1,N_Tbb); %Poissonian photon noise
NEP.wave=zeros(1,N_Tbb); %Wavebunching photon noise
NEP.totphoton=zeros(1,N_Tbb); %Total noise due to photon fluctuations

%Numerical integration of Planck spectral radiance for all blackbody temperatures
%the hroughput is included here for each individual frequency
for p=1:length(Tbb)
    %Photon Occupation Number (Bose-Einstein distribution):
    occupation=1./(exp(h*freq/(k*Tbb(p)))-1);
    %Planck brillance W/(m^2 str Hz) [B_{\nu}(T_{bb})]:
    brilliance=(method.pol*h*freq).* (freq.^2/(c^2)) .*occupation;%also called brightness
    % total irradiation with fixed throughput but without filters in W/Hz
    % [B_{\nu}(T_{bb})*A\Omega]:
    if strcmp(method.tp,'lambda^2')
            Etendue=(c./freq).^2; %array
            method.Etendue=(c/method.centrefreq)^2;
        elseif strcmp(method.tp,'Geometrical')
            Etendue=method.solidangle*pi*(method.lensdiameter/2)^2+zeros(1,length(freq));
            method.Etendue=method.solidangle*pi*(method.lensdiameter/2)^2;
        else
            error('no correct string for thoughput option');
    end
    irradiation{p}=brilliance.*Etendue; 
    % total irradiation @ lens front in W/Hz
    % [B_{\nu}(T_{bb})*A\Omega*F_{\nu}]:
    filterirad2=irradiation{p}.*filter.*method.eta_c;
    % Power per frequency bin in W.
    % [B_{\nu}(T_{bb})*A\Omega*F_{\nu} d\nu]:
    powerbla2{p}=filterirad2*method.freqresolution*1e9;
    % Total received power [W] Summation over power in all F bins
    TotalPbb(p)=sum(powerbla2{p});
    
    % Calculation of the various NEP's
    NEP.g_r(p)=sqrt(sum(method.GR*powerbla2{p}*delta/method.eta_pb)); %g-r noise, only recombination
    NEP.poisson(p)=sqrt(sum(2*powerbla2{p}.*h.*freq)); %Poisson noise term
    
    lambda2=(c./freq).^2; %lamdba^2
    % wave bunching: 2Phf        * opt coupling source - detector * occupation%
    wave=(2*powerbla2{p}.*h.*freq).*(filter.*method.eta_c).*occupation; 
    NEP.wave(p)=sqrt(sum(wave)); %wave contribution
end
NEP.totphoton=sqrt(NEP.g_r.^2+NEP.poisson.^2+NEP.wave.^2); %total KID NEP

%==========================================================================
% Plotting and warp up.
%==========================================================================
if plotdata==1
    %Plot the total power, filter transmission and NEPs if desired.
    plotresults(Tbb,NEP,filter,TotalPbb,filterform,freq,powerbla2,legendstrr);
end

%END OF MAIN PROGRAM
end

function plotresults(Tbb,NEP,filter,powersum,filterform,freq,powerbla2,legendstrr)
%This plotting function shows the results of the blackbody main
%program. It requires 6 variables from this function as input.
%creates one figure
global c;

%FIGURE 
figure(123451);
%clf %clears figure window

%Total received power [W] as a function of Tbb
subplot(2,2,1)
plot(Tbb,powersum)
grid on
hold on
legend('Integrated power')
ylabel('Recieved Power [W]')
xlabel('T_{BB} [K]')

%Overview of the calculated NEP's use to photon induced pair breaking
subplot(2,2,2)
loglog(Tbb,NEP.poisson)
hold on
loglog(Tbb,NEP.wave,'g')
loglog(Tbb,NEP.g_r,'r')
loglog(Tbb,NEP.totphoton,'k')
legend('Poisson NEP','wave NEP','g-r NEP','NEPtot')
ylabel('Photon NEP [W/Hz^1/2]')
xlabel('T_{BB} [K]')
title(['NEP calc. overview ']);
axis tight
grid on

%filter T
subplot(2,2,3)
kolors = colormap(jet(length(filterform)));
semilogy(freq(:)/1e12,filter(:),'color','k','LineWidth',3);hold on
for Nfilters=1:length(filterform)
    semilogy(c.*filterform{Nfilters}(:,1)/1e12,filterform{Nfilters}(:,2),'color',kolors(Nfilters,:),'LineWidth',2)
end
ylabel('Transmission')
xlabel('Frequency (THz)')
axis([0 2 1e-5 1])
legend([{'all'},legendstrr]);
hold off
grid on

%Overview of the power received at each wavelength
subplot(2,2,4)
Tsweep = colormap( jet( ceil(length(Tbb)/10) ) );


nmn=1;

for Nbb=1:10:length(Tbb)
    semilogy(freq(:)/1e12,powerbla2{1,Nbb}(:),'color',Tsweep(nmn,:),'LineWidth',2);hold on
    plottedTBB(nmn)=Tbb(Nbb);
    nmn=nmn+1;
end
ylim([1e-40 1]);grid on
ylabel('B_{\nu}(T_{bb})*A\Omega*F_{\nu} d\nu')
xlabel('Frequency (THz)')
axis tight
hold off
colormap(Tsweep);

h=colorbar('Ticks',[0, 1],...
         'TickLabels',{plottedTBB(1) , plottedTBB(end)});
h.Label.String='T_{BB}[K]';


%END OF plotresults
end

function [filterform, legendstrr]=getfilterform(method,minfrequency,pad,pOs)
% This function uses the "method.filter" component of the "method" struct to
% determine the used filters. It then loads their filter characteristics
% from file. And returns the cell-array "filterform" of length N. Here N is
% the number of filters in the setup. Each cell contains an 2xM_i array with
% [1/lambda filtertransmission]. NB: at the lowest F point a single
% datapoint is patched to all filters: 0 for BPF/HPF and 1 for LPF. In the
% integration interplation will be used between the last real datapoint and
% this point.

% The filter combinations of the following setups are available.
% method.filter,'1_6THz') % 1.6 THz SPICA SAFARI stack) 
% method.filter,'350 GHz 4Filters') %350 GHz BP from cardiff with additional LPF
% method.filter,'850 GHz') %850 GHz BP from cardiff with additional LPF from spacekids%
% method.filter,'350 GHz Deshima' Deshima setup Utrecht ADR, see help above
% in BB11

%Determine path to the filterfiles.
FilterPath = [pad,filesep,'filterfiles',filesep];
if pOs ~= 0
    disp(['BlackBody: Looking for filter files in: '])
    disp([FilterPath]);
end

global c;
% Minimum frequency as 1/lambda, required to patch some filters.
mininvwave = c/minfrequency;


if strcmp(method.filter,'1_6THz') % 1.6 THz SPICA SAFARI stack) 
    % Added JB/PdV 21_9_2012
    bba=flipdim(dlmread([FilterPath,'totalLbandfilters.txt'],'\t'),1);
    filterform{1}(:,1)=bba(:,1)*100;    % from cm-1 to m-1
    filterform{1}(:,2)=bba(:,2);        % include optical transmission estimate EP
    legendstrr='All filters combined';
    
elseif strcmp(method.filter,'350 GHz 4Filters') %350 GHz BP from cardiff with additional LPF
    %BPF holder
    bba=flipdim(dlmread([FilterPath,'W1275_350GHz.dat'],'\t'),1);%BPF Sample Holder
    filterform{1}(:,1)=[bba(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{1}(:,2)=[bba(:,2)' 0];%patcghing 0 transmission at lowest F as this is a full BPF stack
    %LPFs on cold box
    bba=flipdim(dlmread([FilterPath,'W1052 14cm LPE SCUBAII.dat'],'\t'),1);
    filterform{2}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{2}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    bba=flipdim(dlmread([FilterPath,'B386 18cm LPE.dat'],'\t'),1);
    filterform{3}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{3}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    %4K LPF
    bba=flipdim(dlmread([FilterPath,'W969 37cmLPESCUBAII.txt'],'\t'),1);%KLPF BB
    filterform{4}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{4}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    legendstrr={'W1275 BPF holder','W1052 LPF box','B386 LPF box','W696 LPF BB'};

elseif strcmp(method.filter,'350 GHz Deshima') %350 GHz BP from cardiff with additional LPF
    %BPF+LPF holder
    bba=flipdim(dlmread([FilterPath,'K1817_BPFDeshima.txt'],'\t'),1);%BPF Sample Holder
    filterform{1}(:,1)=[bba(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{1}(:,2)=[bba(:,2)' 0];%patcghing 0 transmission at lowest F as this is a full BPF stack
    bba=flipdim(dlmread([FilterPath,'K1785_LPFDeshima.txt'],'\t'),1);%BPF Sample Holder
    filterform{2}(:,1)=[bba(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{2}(:,2)=[bba(:,2)' 0];%patcghing 0 transmission at lowest F as this is a full BPF stack
    %LPFs on cold box
    bba=flipdim(dlmread([FilterPath,'W1052 14cm LPE SCUBAII.dat'],'\t'),1);
    filterform{3}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{3}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    bba=flipdim(dlmread([FilterPath,'B386 18cm LPE.dat'],'\t'),1);
    filterform{4}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{4}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    %4K LPF
    bba=flipdim(dlmread([FilterPath,'W969 37cmLPESCUBAII.txt'],'\t'),1);%KLPF BB
    filterform{5}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{5}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    legendstrr={'K1817 BPF holder','K1785 LPF holder','W1052 LPF box','B386 LPF box','W696 LPF BB',};
    
elseif strcmp(method.filter,'850 GHz') %850 GHz BP from cardiff with additional LPF
    %Added RJ 2013/12/16
    K1979=flipdim(dlmread([FilterPath,'K1979 855GHz BPF.dat'],'\t'),1);
    K1981=flipdim(dlmread([FilterPath,'K1981 1140GHz LPF.dat'],'\t'),1);
    K1980=flipdim(dlmread([FilterPath,'K1980 990GHz LPF.dat'],'\t'),1);
    %B624=flipdim(dlmread([FilterPath,'B624 660GHz HPF.dat'],'\t'),1); We
    %do not have rthis filter, was here for hystorical reasons and is very
    %similar to B588
    B588=flipdim(dlmread([FilterPath,'B588 HPF.dat'],'\t'),1);
    
    %Filters on the 4K blackbody
    filterform{1}(:,1)=[K1980(:,1)'*100 mininvwave];  % Add 1 datapoint at lowest F
    filterform{1}(:,2)=[K1980(:,2)' 1];%LPF, add 1 out of range at low F
    filterform{2}(:,1)=[B588(:,1)'*100 mininvwave];
    filterform{2}(:,2)=[B588(:,2)' 0];%HPF, add 0 out of range at low F
    filterform{3}(:,1)=[K1979(:,1)'*100 mininvwave];
    filterform{3}(:,2)=[K1979(:,2)' 0];%BPF, add 0 out of range at low F
    
    %Filters on the 100mK outer box
    filterform{4}(:,1)=[K1981(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{4}(:,2)=[K1981(:,2)' 1];%LPF, add 1 out of range at low F
    filterform{5}(:,1)=[K1979(:,1)'*100 mininvwave]; 
    filterform{5}(:,2)=[K1979(:,2)' 0];%BPF, add 0 out of range at low F
    
    %Filters on the 100mK sample
    filterform{6}(:,1)=[K1980(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{6}(:,2)=[K1980(:,2)' 1];%LPF, add 1 out of range at low F
    filterform{7}(:,1)=[B588(:,1)'*100 mininvwave];
    filterform{7}(:,2)=[B588(:,2)' 0];%HPF, add 0 out of range at low F
    filterform{8}(:,1)=[K1979(:,1)'*100 mininvwave];
    filterform{8}(:,2)=[K1979(:,2)' 0];%BPF, add 0 out of range at low F
    
    legendstrr={'K1980 LPF BB','B588 HPF BB','K1979 BPF BB','K1981 LPF Box','K1979 BPF Box','K1980 LPF holder','B588 HPF holder','K1979 BPF holder'};

elseif strcmp(method.filter,'650 GHz') %650 GHz BP from cardiff 
    %Added JB 2017/04/05
    B768=flipdim(dlmread([FilterPath,'B768 650GHz BPF.dat'],'\t'),1); 
    K2136 = flipdim(dlmread([FilterPath,'K2136_450um_bp.txt'],'\t'),1); 
    
    %4K LPF (from 350 ghz setup)
    bba=flipdim(dlmread([FilterPath,'W969 37cmLPESCUBAII.txt'],'\t'),1);%KLPF BB
    filterform{1}(:,1)=[bba(:,1)'*100 mininvwave];
    filterform{1}(:,2)=[bba(:,2)' 1];%LPF, add 1 out of range at low F
    
    %Filter on the LT box
    filterform{2}(:,1)=[B768(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{2}(:,2)=[B768(:,2)' 1];%BPF, add 1 out of range at low F
    
    %Filter on the SH
    filterform{3}(:,1)=[B768(:,1)'*100 mininvwave];  % from cm-1 to m-1
    filterform{3}(:,2)=[B768(:,2)' 1];%BPF, add 1 out of range at low F
    
    
    
    legendstrr={'W969 LPF BB','K2136 650GHz BPF','B768 650GHz BPF'};


else
    error('No good filter selection')
end

end
