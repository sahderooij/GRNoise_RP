close all
clear all
clc

method.filter='350 GHz WT1275'; %Cardiff filter
method.lensdiameter=2; %1 or 2 %lens size
method.aperture='ADR small'; %small aperture in the ADR, 0.071rad is the angle to 1 mm radius hole, theta=atan(1mm/14mm), angle to normal, TP=2pi(1-costheta)
method.eta_pb = 0.57; %pairbreaking efficiency.

method6 = method;
method6.filter='350 GHz 4Filters';

c=3.0e8;
%==========================================================================
% Check the change in filter transmission due to incorporating the coupling
% efficiency.
%==========================================================================
Tbb = [3:0.05:50]';
[TotalPbb5b,~,FilterTransmission5b,~]=blackbody_5b(Tbb,1,method,1);
[TotalPbb5,~,FilterTransmission5,~]=blackbody_5(Tbb,1,method,1);
[TotalPbb4,~,FilterTransmission4,~]=blackbody_4(Tbb,1,method,1);
[TotalPbb6,~,FilterTransmission6,~]=blackbody_6(Tbb,1,method6,1);

figure
axes('LineWidth',2,'FontSize',16)
hold on
plot(FilterTransmission4(1,:)/1e9,FilterTransmission4(2,:),'k-')
plot(FilterTransmission5(1,:)/1e9,FilterTransmission5(2,:),'r-')
plot(FilterTransmission5b(1,:)/1e9,FilterTransmission5b(2,:),'b-')
plot(FilterTransmission6(1,:)/1e9,FilterTransmission6(2,:),'g-')
xlabel('Frequency [GHz]')
ylabel('Filter Transmission')
hold off

figure
axes('LineWidth',2,'FontSize',16,'Box','on')
hold on
plot(Tbb,TotalPbb5./TotalPbb4,'k-','LineWidth',2)
plot(Tbb,TotalPbb5b./TotalPbb4,'b-','LineWidth',2)
plot(Tbb,TotalPbb6./TotalPbb4,'g-','LineWidth',2)
xlabel('T [K]')
ylabel('P5/P4 [W/W]')
legend('Original Coupling','Extended Coupling','New Filters')
hold off

Nonzero = find(FilterTransmission4(2,:));
MeanCoupling = sum(FilterTransmission5(2,Nonzero))/sum(FilterTransmission4(2,Nonzero));

%Check the new function of BB7
[TotalPbb7p100,~,Filtertransmission,~]=blackbody_7(Tbb,1,method,1);
method.EpsTrans = 0.70;
[TotalPbb7p70,~,Filtertransmission,~]=blackbody_7(Tbb,1,method,1);

TotalPbb7p70./TotalPbb7p100

%==========================================================================
% Calculate the power for an absorber
%==========================================================================
[TotalPbbEST,~,~,~]=blackbody_6(Tbb,1,method,1);

methodabs = method;
methodabs.aperture='ADR small Geometrical';
[TotalPbbABS,~,~,~]=blackbody_6(Tbb,2,methodabs,1);
%TotalPbbABS = TotalPbbABS /(pi*1e-6)*(4e-6); %Circular to Square Absorber

[Tbb,TotalPbbEST',TotalPbbABS']

