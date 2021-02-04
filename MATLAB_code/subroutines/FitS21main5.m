function [fres,Q,S21min,fitevaluation,errors] = FitS21main5(S21data,FITF0,UserGuesses)
% This function takes S21data and fits the desired functions to it to
% obtain the resonance frequency, Q factors, maximum dip depth. The fitting
% is done to the S21data using the fitting function specified by FITF0.
%
% The fits currently implemented / reserved are:
% 0 == use minimum value of measured |S21| to determine fres and S21min.
%       Fit is only used to determine Q.
% 1 == use a simple pure Lorentzian Fit to obtain Q Fres and S21min. (FitS21_3 by Pieter de Visser)
%       It uses a square space fit to get the ini parameters nicely.
% 2 == use the function by Khalil et.al. to fit to |S21|(dB) (FitS21_2)
% 3 == [FUTURE] use the function by Khalil et.al. to fit to S21 (complex number, FitS21_3)
% 4 == [FUTURE] make a point-by-point fit in the complex plane. AND/OR make
%               use of theta(F) to determine bias point.
%
% NOTE: To add your own fitting function, just write the desired fitting
% subroutine at the end of the document. And add and FITF0 value in the
% ifelse chain in the main routine.
%
%INPUT:
%   S21data == A Mx3 double array containing the transmission data of the
%   resonator. [F,|S21|,phase(S21)] |S21| should be given in magnitude space.
%   phase(S21) should be given in radians.
%   FITF0 = selects the fitting function. See above.
%   UserGuesses = [Optional] A 5 element vector with [fres,S21min,Ql,Qi,Qc].
%                  Note that the units must be the same as S21data.
%
%OUTPUT:
%   Fres == resonance frequency from the fit
%   Qvector == A 3 element vector containing: [Ql,Qi,Qc]
%   S21min == |S21| @ Fres. As per the evaluation of the fit (dB)
%   fitevaluation == an Nx2 array containing the evaluation of the fit made
%                    in this function. [Freq,|S21|] of the fit made.
%   errors == [u(fres),u(S21min),u(Ql),u(Qi),u(Qc)]
%
%SUBROUTINES:
%   EstimateInitials
%   FitS21_SRON
%   FitS21_Khalil
%
%REQUIRED FILES:
%
%VERSION: 1.0
%   V1.0 (2012-10-24,RJ): Made this into a combination code from which
%                         various fitting modes are called.
%   V2.0 (2012-11-15,RJ): Modified routine to get the common estimation of
%                         initial values for all routines equal. Renamed
%                         some subroutines for clarity
%   V3.0 (2013-02-14,RJ): Changed the Yates/Baselmans fit to a simple
%                         lorenztian with no skewing etc.
%   V4.0 (2013-11-20,RJ): Included an estimate of the error in all derived parameters. 
%   V5.0 (2014-02-18):      Corrected YdB error (FITFO=1 case) and added
%                           (with comments) a simple p[lot for debugging. Uncommenting switches it
%                           on. (line 250
%DATE: October 24, 2012
%AUTHOR: Reinier Janssen
%==========================================================================
%Set output format to screen
format('long','e')

%Catch an a FITF0 value defined outside useful range. Set to default
%(FITF0 == 1).
if (FITF0 < 0)||(FITF0 > 2)
    FITF0 = 1;
end

%Evaluate the base Temperature fit at 100 times measurement resolution
dF = (S21data(end,1)-S21data(1,1))/100/length(S21data(:,1));
F = S21data(1,1):dF:S21data(end,1);

if nargin == 2
    %No initial guesses specified by the user. Internally estimate them.
    [UserGuesses] = EstimateInitials(S21data);
end

if FITF0 == 0
    %A FIT IS ONLY APPLIED TO DETERMINE Q. FIT RESULTS FOR Fres AND
    %S21min ARE NOT USED. Especially useful for shallow dips.
    [~,Q,~,fitresult,FitErrors] = FitS21_PdVLorentz(S21data(:,1:3),UserGuesses);
    YdB = feval(fitresult,F);
    
    %write the resonance frequency and S21min from the raw data values
    [S21min,F0index] = min(S21data(:,2)); %get location and value of S21min (resonance)
    S21min = 20*log10(S21min); %Convert to dB.
    fres = S21data(F0index,1); %Resonance frequency
    
    errors = FitErrors;
    errors(1,1) = abs(S21data(F0index-1,1)-S21data(F0index+1,1))/2;
    uS21min = abs(0.5*(S21data(F0index-1,2)+S21data(F0index+1,2))-S21data(F0index,2));
    errors(2,1) = 20*uS21min/S21data(F0index,2)/log(10);
    errors(3:5,1) = FitErrors(3:5,1);
    %==================================================================

elseif FITF0 == 1
    %SIMPLE LORENTZIAN MULTISTAGE FIT BY DE VISSER. Only a narrow band fit
    %around the dip generates the output, broad band fitis used only to get
    %the guesses right. Fit is done in dB's!
    [fres,Q,S21min,fitresult,FitErrors] = FitS21_PdVLorentz(S21data(:,1:3),UserGuesses);
    YdB = feval(fitresult,F);
    
    errors = FitErrors;
    %==================================================================

elseif FITF0 == 2
    %FIT TO |S21| BY KHALIL ET.AL.
    [fres,Q,S21min,fitresult,FitErrors] = FitS21_Khalil(S21data(:,1:3),UserGuesses);
    YdB = feval(fitresult,F);
    
    errors = FitErrors;
    %==================================================================
elseif FITF0 == 3
    %FIT TO S21 BY KHALIL ET.AL.
    %==================================================================
end

%Prepare the evaluated fit function for export
fitevaluation = zeros(length(F),2);
fitevaluation(:,1) = F;
fitevaluation(:,2) = YdB;

%======================================================================
% END OF FitS21main MAIN ROUTINE
%======================================================================
end

function  [Fres,Qvector,S21min,resultFres,FitErrors] = FitS21_PdVLorentz(S21data,UserGuesses)

% Similar as S21_1 (FitS21_SRON) but with a pure Lorentzian fit, without the skewing
%
% This function takes S21data and fits a pure Lorentzian resonance dip to 
% this data to obtain Q. A parabolic fit is then performed to obtain a more
% accurate vaule of Fres and S21min.
%
%
%INPUT:
%   S21data == A Mx2 double array containing the transmission data of the
%   resonator. [F,|S21|] |S21| should be given in magnitude space.
%   UserGuesses == [Optional] A 3 element vector containing guesses for [Fres,Q,S21min]
%
%OUTPUT:
%   Fres == resonance frequency from the fit
%   Qvector == A 3 element vector containing: [Ql,Qi,Qc]
%   S21min == |S21| @ Fres. As per the evaluation of the fit (dB)
%   resultFres == Direct output of the "fit" function. The fit for the Q is
%   given as an output here. It also contains the slope.
%   FitErrors = [u(Fres),u(S21min),u(Ql),u(Qi),u(Qc)] errors of output
%   parameters
%
%SUBROUTINES:
%
%REQUIRED FILES:
%
%VERSION: 2.2
%   Original version 1 by Stephen Yates / Jochem Baselmans. Fitted
%   completely different function
%   V2.0 (2012-08-24/28,RJ): addition of many comments, general code cleaning 
%                            without functionality change, added section
%                            that it can guess its own initial values.
%   V2.1 (2012-12-13, PdV):  - Removed the stretch parameter. 
%                            - The slope is now divided out before fitting, which saves an additional fit parameter.
%                            - All three outputs: Q,fres and S21min are now based on the fit over the small range around fres, to have
%                              a set of consistent fit parameters that correspond to one pure Lorentzian.
%                            - The slope is now divided instead of subtracted like in FitS21_1 (which is wrong!!!!)
%   V2.2 (2013-11-20, RJ): Modified the code to output errors on fit
%                          parameters. Note to self: Need to finish for
%                          Khalil fit.
%DATE: December 12, 2012
%AUTHOR: Pieter de Visser
%==========================================================================
% Set output format to screen
format('long','e');
% If bandwidth = Fres/Ql
% Limit the data to which a big fit is made (for Q) to Fres +/- 0.5*bandwidth*Nbw1
Nbw1 = 2;
% Limit the data to which a small fit is made (for Fres) to Fres +/- 0.5*bandwidth*Nbw2
Nbw2 = 0.4;

%==========================================================================
% Check if initial guesses are present and if not make them internally.
%==========================================================================
if nargin == 1 %Only S21data
    UserGuesses = EstimateInitials(S21data);
end

fres_g = UserGuesses(1);
S21min_g = UserGuesses(2);
Q_g = UserGuesses(3);

%==========================================================================
% Select the bandwidth limited data for the fits
%==========================================================================
%Determine the approximate bandwidth of the resonator.
BW_g = fres_g/Q_g;

%Select data for the Lorentzian fit.
BigFitData = S21data(:,1) >= fres_g - 0.5*BW_g*Nbw1 & S21data(:,1) <= fres_g + 0.5*BW_g*Nbw1;
%For the fit of Q a minimum of 5 points is set as requirement
if sum(BigFitData) < 5
    %Not enough points. Revert to bigger bandwidth fit.
    fprintf('Warning FitS21_PdVLorentz: bandwidth guess does not leave enough points for fitting\n');
    fprintf('Using entire S21 curve for fitting Q instead.\n');
    BigFitData(:,:) = 1;
end

%Select data for the parabolic fit
SmallFitData = S21data(:,1) >= fres_g - 0.5*BW_g*Nbw2 & S21data(:,1) <= fres_g + 0.5*BW_g*Nbw2;
%For the fit of Fres a minimum of 10 points is set as requirement
if sum(SmallFitData) < 10
    %Not enough points. Revert to bigger bandwidth fit.
    fprintf('Warning FitS21_PdVLorentz: bandwidth guess does not leave enough points for fitting\n');
    fprintf('Using entire S21 curve for fitting Fres and |S21|(Fres) instead.\n');
    SmallFitData(:,:) = 1;
end

%==========================================================================
% (PdV 13-12-12: this fit is now not really necessary anymore, I have only 
% used it now to guess the Q for the next fit)
%==========================================================================

%Guess a linear background slope (drift)
slope_g = (S21data(1,2)^2 - S21data(end,2)^2)/(S21data(1,1) - S21data(end,1));%the ^2 here is again because of power space

%Set fit options
s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[fres_g Q_g S21min_g],...
    'Lower',[fres_g*0.9 Q_g/2 S21min_g/4],...
    'Upper',[fres_g*1.1 Q_g*2 S21min_g*2],...
    'MaxFunEvals',1000);

%Fit Lorentzian in power space with linear background already divided out
%I checked that this is the right equation for |S21|^2 (PdV)
ftype = fittype(['(1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2)).*(' num2str(S21data(1,2).^2) '+' num2str(slope_g) '*(x-' num2str(S21data(1,1)) '))'],'options', s); 
%Perform the actual fit
FitResult=fit(S21data(BigFitData,1),S21data(BigFitData,2).^2,ftype); %the ^2 is because of the 'power space'

%=============================================================================
% Fit a simple Lorentzian in dB space to find F0, S21min and Q
% Also now get Q from this method (you want to know the Q close to resonance!)
%=============================================================================
%Specify fit options
s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[fres_g FitResult.Q S21min_g],...
    'MaxFunEvals',1000);
%Fit Lorentzian in logspace close to resonance (SmallFitData), with linear
%background divided out (slope is here a predetermined number, not a fit parameter
%JB:mach 2014: Checked against the Mazin Eqn and gives the same result
ftype = fittype(['10*log10((1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2)).*(' num2str(S21data(1,2).^2) '+' num2str(slope_g) '*(x-' num2str(S21data(1,1)) ')))'],'options', s);
%Perform the actual fit
resultFres=fit(S21data(SmallFitData,1),20*log10(S21data(SmallFitData,2)),ftype);

%Help plotting
% figure(10000)
% plot(S21data(BigFitData,1),S21data(BigFitData,2).^2,'ro');hold on;
% plot(FitResult,'r');
% plot(S21data(BigFitData,1),20*log10(S21data(BigFitData,2)),'ro');hold on;
% plot(S21data(SmallFitData,1),20*log10(S21data(SmallFitData,2)),'ko','MarkerFaceColor','k');hold on;
% plot(resultFres,'k');
% legend('data', 'broad fit (used for input finding)', 'Realfit, done with narrow range data');
% hold off;

%==========================================================================
% Determine interesting parameters form fits and prepare for output
%==========================================================================
%Now only use the values from the 'small' fit range to have a consistent
%set of values from one Lorentzian fit:
Fres = resultFres.Fr;
Q = resultFres.Q;
S21min = resultFres.Smin;

Qi = Q/S21min;
Qc = Qi*Q/(Qi-Q);
Qvector = [Q,Qi,Qc];

%==========================================================================
% Determine the errors and prepare for output
%==========================================================================
ConfIntFit = confint(resultFres);
ufres = (Fres-ConfIntFit(1,1))/1.96;
uQ = (Q-ConfIntFit(1,2))/1.96;
uS21min = (S21min-ConfIntFit(1,3))/1.96;

uQi = Qi*sqrt((uQ/Q)^2+(uS21min/S21min)^2);
uQc = Qc^2*sqrt((uQ/Q^2)^2+(uQi/Qi^2)^2);

FitErrors = zeros(5,1);
FitErrors(1,1) = ufres;
FitErrors(2,1) = 20*uS21min/S21min/log(10);
FitErrors(3:5,1) = [uQ,uQi,uQc];

%==========================================================================
%Prepare S21min for output in dB
%==========================================================================

S21min = 20*log10(S21min); %output in dB

end%function

function  [Fres,Qvector,S21min,FitResult,FitErrors] = FitS21_SRON(S21data,UserGuesses)

% This function takes S21data and fits a skwed Lorentzian resonance dip to 
% this data to obtain Q. A parabolic fit is then performed to obtain a more
% accurate vaule of Fres and S21min.
%
%
%INPUT:
%   S21data == A Mx2 double array containing the transmission data of the
%   resonator. [F,|S21|] |S21| should be given in magnitude space.
%   UserGuesses == A 3 element vector containing guesses for [Fres,Q,S21min]
%
%OUTPUT:
%   Fres == resonance frequency from the fit
%   Qvector == A 3 element vector containing: [Ql,Qi,Qc]
%   S21min == |S21| @ Fres. As per the evaluation of the fit (dB)
%   FitResult == Direct output of the "fit" function. The fit for the Q is
%   given as an output here.
%
%SUBROUTINES:
%
%REQUIRED FILES:
%
%VERSION: 2.0
%   Original version 1 by Stephen Yates / Jochem Baselmans. Fitted
%   completely different function
%   V2.0 (2012-08-24/28,RJ): addition of many comments, general code cleaning 
%                            without functionality change, added section
%                            that it can guess its own initial values.
%   V2.1 (2012-11-15,RJ): Make use of EstimateInitials Routine
%
%DATE: August 28, 2012
%AUTHOR: Reinier Janssen
%==========================================================================
% Set output format to screen
format('long','e');
% If bandwidth = Fres/Ql
% Limit the data to which a big fit is made (for Q) to Fres +/- 0.5*bandwidth*Nbw1
Nbw1 = 2;
% Limit the data to which a small fit is made (for Fres) to Fres +/- 0.5*bandwidth*Nbw2
Nbw2 = 0.3;

%==========================================================================
% Check if initial guesses are present and if not make them internally.
%==========================================================================
if nargin == 1 %Only S21data
    UserGuesses = EstimateInitials(S21data);
end

fres_g = UserGuesses(1);
S21min_g = UserGuesses(2);
Q_g = UserGuesses(3);


%==========================================================================
% Select the bandwidth limited data for the fits
%==========================================================================
%Determine the approximate bandwidth of the resonator.
BW_g = fres_g/Q_g;

%Select data for the skew Lorentzian fit.
BigFitData = S21data(:,1) >= fres_g - 0.5*BW_g*Nbw1 & S21data(:,1) <= fres_g + 0.5*BW_g*Nbw1;
%For the fit of Q a minimum of 5 points is set as requirement
if sum(BigFitData) < 5
    %Not enough points. Revert to bigger bandwidth fit.
    fprintf('Warning FitS21_1: bandwidth guess does not leave enough points for fitting\n');
    fprintf('Using entire S21 curve for fitting Q instead.\n');
    BigFitData(:,:) = 1;
end

%Select data for the parabolic fit
SmallFitData = S21data(:,1) >= fres_g - 0.5*BW_g*Nbw2 & S21data(:,1) <= fres_g + 0.5*BW_g*Nbw2;
%For the fit of Fres a minimum of 10 points is set as requirement
if sum(SmallFitData) < 10
    %Not enough points. Revert to bigger bandwidth fit.
    fprintf('Warning FitS21_1: bandwidth guess does not leave enough points for fitting\n');
    fprintf('Using entire S21 curve for fitting Fres and |S21|(Fres) instead.\n');
    SmallFitData(:,:) = 1;
end

%==========================================================================
% Perform a skew Lorentzian fit to obtain Q
%==========================================================================

%Guess a linear background slope (drift)
slope_g = (S21data(1,2)^2 - S21data(end,2)^2)/(S21data(1,1) - S21data(end,1));

%Estimate the skewterm (incl error catching)
subF3dBpoint = find(S21data(:,1)>(fres_g-0.5*BW_g),1);
if isempty(subF3dBpoint)
    subF3dBpoint = 1;
end
superF3dBpoint = find(S21data(:,1)>(fres_g+0.5*BW_g),1);
if isempty(superF3dBpoint)
    superF3dBpoint = length(S21data(:,1));
end
%actual estimate skew parameter
stretch=(S21data(subF3dBpoint,2)^2-S21data(superF3dBpoint,2)^2)/BW_g; 

%Set fit options
s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[fres_g Q_g S21min_g slope_g stretch],...
    'Lower',[fres_g*0.9 Q_g/2 S21min_g/4 -25*abs(slope_g) -25*abs(stretch)],...
    'Upper',[fres_g*1.1 Q_g*2 S21min_g*2 abs(slope_g)*25 +25*abs(stretch)],...
    'MaxFunEvals',1000);

%Fit Lorentzian like in power space with linear background and stretch
ftype = fittype('(1-(1-Smin^2)*(1+stretch*(x-Fr)/Fr)/(1+(2*Q*(x-Fr)/Fr)^2)+l*(x-Fr)/Fr)','options', s);  
%Perform the actual fit
FitResult=fit(S21data(BigFitData,1),S21data(BigFitData,2).^2,ftype);

%==========================================================================
% Fit a simple Lorentzian in dB space to find F0 and S21min
%==========================================================================
%Specify fit options
s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[fres_g Q_g S21min_g],...
    'Lower',[fres_g*0.9 Q_g/2 S21min_g/4],...
    'Upper',[fres_g*1.1 Q_g*2 S21min_g*2],...
    'MaxFunEvals',100);
%Specify function to fit
ftype = fittype('10*log10(1-(1-Smin^2)/(1+(2*Q*(x-Fr)/Fr)^2))','options', s);
%Perform the actual fit
[resultFres]=fit(S21data(SmallFitData,1),20*log10(S21data(SmallFitData,2)),ftype);
%==========================================================================
% Determine interesting parameters form fits and prepare for output
%==========================================================================
Fres = resultFres.Fr;
Q=FitResult.Q;
S21min = resultFres.Smin;

Qi = Q/S21min;
Qc = Qi*Q/(Qi-Q);
Qvector = [Q,Qi,Qc];

S21min = 20*log10(S21min); %output in dB
end

function [Fres,Qvector,S21min,FitResult,FitErrors] = FitS21_Khalil(S21data,UserGuesses)

% This function takes S21data and fits a Lorentzian resonance dip to this
% data. The fit allows for assymmetry due to mismatches in the transmission
% line. As defined by Khalil et.al., J APPL PHYS 111, 054510 (2012)
% This routine fits equation 12 of this paper to the full S21data set.
% From the fit it returns key resonantor parameters as well as the fit itself.
%
% In order to perform the fit a number of guesses for resonance frequency
% and quality factors are required. These can be either user supplied or
% are calculated internally.
%
% NOTE: The assymmetric is not appropriate for the effect induced by
% overdriving the resonator.
% NOTE: In addition to correcting for assymmetry a correction for an
% overall linear slope is made. This can be detrimental if the S21data does
% not contain a lot of transmission besides the resonance dip.
%
%INPUT:
%   S21data == A Mx2 double array containing the transmission data of the
%   resonator. [F,|S21|] |S21| should be given in magnitude space.
%   UserGuesses == [Optional] A 3 element vector containing guesses for [Fres,Q,Qc]
%
%OUTPUT:
%   Fres == resonance frequency from the fit
%   Qvector == A 3 element vector containing: [Ql,Qi,Qc]
%   S21min == |S21| @ Fres. As per the evaluation of the fit (dB)
%   FitResult == Direct output of the "fit" function.
%
%SUBROUTINES:
%
%REQUIRED FILES:
%
%VERSION: 2.0
%   Original version 1 by Stephen Yates / Jochem Baselmans. Fitted
%   completely different function
%   V2.0 (2012-08-20,RJ): Large rewrite, incl. addition of many comments,
%                         changed function that is fitted. Added section
%                         that it can guess its own initial values.
%   V2.1 (2012-11-15,RJ): Make use of EstimateInitials Routine
%
%DATE: August 20, 2012
%AUTHOR: Reinier Janssen
%==========================================================================
format('long','e'); %Set the way things are written to screen

%==========================================================================
% Check if initial guesses are present and if not make them internally.
%==========================================================================
if nargin == 1 %Only S21data
    UserGuesses = EstimateInitials(S21data);
end

fres_g = UserGuesses(1);
Q_g = UserGuesses(3);
Qc_g = UserGuesses(5);

%Determine the approximate bandwidth of the resonator.
BW_g = fres_g/Q_g;

%Guess the linear slope using the first and last point.
slope_g = (S21data(end,2)-S21data(1,2))/(S21data(end,1)-S21data(1,1));

%==========================================================================
% Performing the actual NLLSq fit using standard matlab routines
%==========================================================================
% Set fit options: tollerances and max,min,guess value for all fit
% parameters.
fitopt = fitoptions('method','NonlinearLeastSquares',...
    'TolX',1e-16,'TolFun',1e-16,...
    'Lower',[0.1*Q_g 0.1*Qc_g -0.5*pi fres_g-3*BW_g -25*abs(slope_g)],...
    'Upper',[10*Q_g 10*Qc_g 0.5*pi fres_g+3*BW_g 25*abs(slope_g)],...
    'Startpoint',[Q_g Qc_g 0 fres_g slope_g]);
% %There are abs(slope) to avoid that LowerBound > UpperBound

% Define the function that needs to be fit
fittp = fittype('20*log10(abs(  1- ((Q/Qe*exp(1i*phi))/(1+2i*Q*(freq-fres)/fres)) +slope*(freq-fres) ))',...
    'dependent',{'Magnitude'},'independent',{'freq'},...
    'coefficients',{'Q','Qe','phi','fres','slope'});

%Perform the actual fit
FitResult = fit(S21data(:,1),20*log10(abs(S21data(:,2))),fittp,fitopt);

%==========================================================================
% Use the FitResult to determine key resonator parameters.
%==========================================================================
%Extract the values of the fit parameters
ParValues = coeffvalues(FitResult); %[Q,|Qe|,phi,Fres,slope] (CORRECT ????????)

%Resonance Frequency
Fres = ParValues(4);

%Quality Factors
Ql = ParValues(1);
Qc = ParValues(2);
Qi = 1/((1/Ql)-(1/Qc)*cos(ParValues(3)));
Qvector = [Ql Qi Qc];

%S21min
S21min = abs(1-Ql/Qc*exp(1i*ParValues(3))); %mag

%==========================================================================
% Determine the errors and prepare for output
%==========================================================================
ConfIntFit = confint(resultFres);
ufres = (fres-ConfIntFit(1,1))/1.96;
uQ = (Q-ConfIntFit(1,2))/1.96;
uS21minMag = (S21min-ConfIntFit(1,3))/1.96;

uQi = Qi*sqrt((uQ/Q)^2+(uS21min/S21min)^2);
uQc = Qc^2*sqrt((uQ/Q^2)^2+(uQi/Qi^2)^2);

FitErrors = zeros(5,1);
FitErrors(1,1) = ufres;
FitErrors(2,1) = 20*uS21min/S21minMag/log(10);
FitErrors(3:5,1) = [uQ,uQi,uQc];

%==========================================================================
%Prepare S21min for output in dB
%==========================================================================

S21min = 20*log10(S21min); %output in dB
end

function [Guesses] = EstimateInitials(S21data)
%This function estimates resonance frequency, resonance dip depth and
%loaded quality factor based on the provided resonance data.
%Guesses will be a vector containing: [fres,S21min,Ql,Qi,Qc]
%Requires S21data in magnitude space. Output (fres,S21min) will be in same
%units as input.

%2013-02-14, RJ: Removed minor filtering from guestimation.

%==========================================================================
% Estimate the resonance frequency using max dip depth
%==========================================================================
[S21min_g,F0index] = min(S21data(:,2)); %get S21 at resonance in magnitude space
fres_g=S21data(F0index,1); %get resonance frequency guess in GHz!!!
%==========================================================================
% Estimate the initial value for loaded Q.
%==========================================================================
% Determine the index at which the resonance dip is half its depth.
% See also Ben Mazin, Thesis page 29, eq 2.44
[~,bwindex]=min(abs(S21data(:,2).^2-((S21min_g.^2+1)/2)));
if bwindex-F0index==0 %half dip depth is same as resonance.
    %Assume the half dip depth is located 1 index away from resonance.
    Q_g=fres_g/(abs(2*(S21data(F0index+1,1)-S21data(F0index,1)))); %Q=f/df, df whole dip bandwitdh
else
    %Loaded Q related to "FWHM" of resonce dip.
    Q_g=fres_g/(abs(2*(S21data(bwindex,1)-S21data(F0index,1)))); %Q=f/df, df whole dip bandwitdh
end
%==========================================================================
% Estimate the initial value of the coupling Q
%==========================================================================
%First guess Qi using Ben Mazin, Thesis page 29, eq 2.42 and below 2.43
Qi_g = Q_g/S21min_g;
%Then guess Qc using Qi_g and Q_g
Qc_g = (Q_g*Qi_g)/(Qi_g - Q_g);

%==========================================================================
%Prepare output
%==========================================================================
Guesses = zeros(5,1);
Guesses(1,:) = fres_g;
Guesses(2,:) = S21min_g;
Guesses(3,:) = Q_g;
Guesses(4,:) = Qi_g;
Guesses(5,:) = Qc_g;

end
