function [fDSB,pDSB]=SSB2DSB(fSSB,pSSB)
%returns DSB value, expects f to be 0 at first index and near sampling
%frequency at max index
%WORKS FOR CROSS-COVARIANCE.

%example: N=3 @ 1 Hz sampling freq
%fSSB=[0 1/3 2/3] Hz, pSSB=[p0 p1 p2]
%fDSB=[0 1/3] Hz, pDSB=[p0 p1+p2]

%example: N=4 @ 1 Hz sampling freq
%fSSB=[0 1/4 2/4 3/4] Hz, pSSB=[p0 p1 p2 p3]
%fDSB=[0 1/4 2/4] Hz, pDSB=[p0 p1+p3 p2]

%example: N=5 @ 1 Hz sampling freq
%fSSB=[0 1/5 2/5 3/5 4/5] Hz, pSSB=[p0 p1 p2 p3 p4]
%fDSB=[0 1/5 2/5] Hz, pDSB=[p0 p1+p4 p2+p3]

N=length(fSSB);


m=fix(N/2)+1; %number of elements to return in DSB

fDSB=fSSB(1:m);

t=flipdim(pSSB,2);
t=[0 t];
u=[pSSB 0];

pDSB=u+conj(t); %for cross covariance: psd(w)=psd(-w)*
pDSB=pDSB(1:m);

%we did do one thing to much for the even case:
if mod(m,2)==0
    pDSB(m)=pSSB(m);
end