function [tau, noiselevel, setupnoise] = gettau_noise(fdata,noisedata,indfref,blawindow,maxtau,plotfigure,taustart)
% tries to get the lifetime, if not ok, retursm tau=taumax
% fdata is post detection frequency in Hz, 1 col
% noisedata is NoisePSD in dBc/Hz,1 col
% blawindow is amount of ponts (+- blawindiow_ around indfref to use for
% Ref noise
% indfref is the index of ref value to get noise level
% tau in msec.
%close all
%plotfigure = 1;
if nargin == 6
    taustart = 0.5e-3;
end
%get phase noise PSD photons
if blawindow <=1
    blawindow=5;
    disp('NEP ref range set to default (5 points)')
end
if ~isempty(indfref) && indfref-blawindow > 0 % fit can proceed
    tnd = 10.^(noisedata((indfref-blawindow):(indfref+blawindow))/10);
    noiselevel=mean(tnd); % ref Hz value
else
    tau = maxtau;
    return
end
%get setup level
rtuF=find(fdata < 3e5 & fdata > 1000);
[~,thf2]=min(noisedata(rtuF));%we take the miniumu value up to 3e5Hz
highfind2=rtuF(thf2);
setupnoise = 10^(noisedata(highfind2)/10);
if noiselevel / setupnoise < 1.1
    tau = maxtau;
    return
else
    noiselevel = noiselevel - setupnoise;
end

%fit
s = fitoptions('Method','NonlinearLeastSquares',...
    'Startpoint',[taustart],...
    'MaxFunEvals',200,...
    'Lower',0);
ftype = fittype(['10*log10(' num2str(noiselevel) './(1 + (2*pi*tau*x).^2)+ ' num2str(setupnoise) ')'],'options',s);
result = fit(fdata,noisedata,ftype);
if result.tau <= maxtau
    tau = result.tau;
else
    tau = maxtau;
end

if plotfigure == 1
    figure(2001)
    semilogx(fdata,noisedata,'b');hold on;
    %semilogx(fdata,10*log10(noiselevel./(1 + (2*pi*taustart*fdata).^2)+setupnoise),'--k');
    semilogx(fdata,10*log10(noiselevel./(1 + (2*pi*tau*fdata).^2)+setupnoise),'-k');
    disp('bla');

end