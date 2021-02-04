
function[powersum,NEP,BBcal,method]=blackbody_int(TTT,BBcal,method)
% calls external function blackbody_2.m if BBcal.T=0 (first KID in main
% program). backbody.m Calculates the powers vs black body T at many
% temperatures and puts the results in BBcal struct.
% afterwards, and also if BBcal.T~=0 the script inperpoaltes the BBcal
% struct ti get parameters.

% TTT col of temperatures
% BBcal struct with tempertaure, power and photon noise NEP
% pol polarisation 1 or 2
% method: string describing the optics of setup (see blackbody.m)
if isequal(BBcal.T,0) % if BB.cal is empty call blackbody.m
    BBcal.T=[2:0.05:30 31:1:60]';
    [BBcal.power,~,BBcal.NEP,method]=blackbody_14(BBcal.T,method,1); % [TotalPbb,Filtertransmission,NEP,method]
end

% interpolate BBcal struct, log space
powersum=10.^(interp1(log10(BBcal.T),log10(BBcal.power),log10(TTT)));

NEP.poisson=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.poisson),log10(TTT)));%2phF
NEP.g_r=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.g_r),log10(TTT)));%2PDelta/eta
NEP.wave=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.wave),log10(TTT)));%wave term
NEP.totphoton=10.^(interp1(log10(BBcal.T),log10(BBcal.NEP.totphoton),log10(TTT)));
end
