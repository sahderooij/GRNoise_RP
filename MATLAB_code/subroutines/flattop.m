function [flattop]=flattop(n)
%generates flattop window
x=(0:n-1) /(n-1);
t=2*pi*x;

a0=0.215578948;
a1=0.41663158;
a2=0.277263158;
a3=0.083578947;
a4=0.006947368;

flattop=(a0-a1*cos(t)+a2*cos(2*t)-a3*cos(3*t)+a4*cos(4*t))';