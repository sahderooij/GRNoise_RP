function [tau,level,setupnoise,taumin,taumax] = crossfit(fdata,noisedata,GRlevel_start,taumax)
%fit spectrum with single lorentzian roll-off
% fdata = frequency data in Hz
% noisedata is crossPSD LINEAR not log
% GRnoiselevek is the mean level of the GR noise (used as fit start)
% tau start in sec is the lifetime estimate used as fit start
if nargin == 3
    taumax = 1e-3;
end
%define bounds for the fit
snf = 1.5;% factor to define the fit range setuop noise: Lower = setupnoise_start/snf, Upper = setupnoise_start*snf%
GRnf = 1.5;% factor to define the fit range of GR noise: Lower = GRlevel_start/snf, Upper = GRlevel_start*snf%
taunf = 300;% factor to define the fit range of lifetime: Lower = taumax/taunf, Upper = taumax%
%data prepare: normalize noisedata to GRnoise level
noisedata = noisedata/GRlevel_start; 
setupnoise_start = mean(noisedata(end-10:end)); %setupnoise estimate
s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',   [1, setupnoise_start, taumax],...
    'Lower',        [1/GRnf, setupnoise_start-snf*abs(setupnoise_start), taumax/taunf],...
    'Upper',        [1*GRnf, setupnoise_start+snf*abs(setupnoise_start), taumax],...
    'MaxFunEvals',200);
ftype = fittype('level./(1 + (2*pi*tau*x).^2)+ setupnoise ','options',s);%coeffnames(ftype) = level, setupnoise, tau
% Fit
result = fit(fdata,noisedata,ftype);
% mappoing result, and renormalizing the levels to the input data
tau = result.tau;
if tau == taumax || tau == taumax/taunf %bounded
    result.tau = NaN;
end
cfi = confint(result);
taumin = cfi(1,3);taumax = cfi(2,3);
level = result.level*GRlevel_start;
setupnoise = result.setupnoise*GRlevel_start;


% % %test plot
% % semilogx(fdata,noisedata,'r');hold on;
% % semilogx(fdata,  1./(1+(2*pi*taumax*fdata).^2) + setupnoise_start,'-k'  );
% % semilogx(fdata,result(fdata),'b')
end