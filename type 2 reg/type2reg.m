clear all; close all;

%% Regulator values
R1 = 200e3; 
R2 = 8.06e6;
C1 = 100e-9;
C3 = 1e-9;

%% Boost Values

Vbatt = 16.4;
Vin = 9; % 9v is about mppt
I = 0.5; %charge current (Available current at MPP)
L = 4.7e-6; %% 4.7uH
C = 49.4e-6; %% ~50 uF output capacitance
R = Vbatt/I; %% Effective resistance

fs = 55e3;
D = (1 - (Vin/Vbatt))
K = (2*L)/((Vbatt/I)*(1/fs))
M =  Vbatt/Vin %(1+sqrt(1+(4*(D^2/K))) ) / 4
D_prime = 1 - D; % 16% duty, boosting from 9 to 10.5 (Vbat low) 82% is Vbat high (16.4)

%%PWM values:
PWM = 4.5

%% bode options
options = bodeoptions;
options.FreqUnits = 'Hz'; % or 'rad/second', 'rpm', etc

figure(1)

%% CCM Boost model
subplot(2,3,2)

Gd0_ccm = Vin/D_prime; %% DC gain
wz_ccm = (D_prime^2*R)/L;
w0_ccm = D_prime/(sqrt(L*C));
Q_ccm = D_prime*R*sqrt(C/L);

Gvd_ccm = tf([-1/wz_ccm 1], [(1/w0_ccm)^2, 1/(Q_ccm*w0_ccm), 1]);
Gpwm = 1/PWM;


%%open loop
Tu_ccm = Gd0_ccm*Gvd_ccm*Gpwm;
handle = bodeplot(Tu_ccm,options);
p=getoptions(handle);    
p.Ylim{1}= [-40 50]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Ylim{2}= [-225 10]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Title.String = 'CCM Uncompensated Loop Response';
p.Title.FontWeight = 'Bold';
%p.Title.FontSize = 20;
p.PhaseMatching = 'on';
setoptions(handle,p); %update your plot

%% Battery model
Batt = tf([0.01*0.093 1],[0.093 0]); %% model of a battery. 1F / 150Ah capacity. 3500mAH * 4 = 14Ah. Hence, 1/150 * 14 = 0.093F. 10mohms in parallel


%% DCM Boost model

P = I*(Vbatt - Vin);
Re = (Vin)^2/P;
D = sqrt((2*L*P*fs)/(Vin^2))
M = Vbatt/Vin;

subplot(2,3,5)

Gd0_dcm = ((2*Vbatt)/D) * ((M-1)/((2*M)-1));

fp = (2*M-1)/(2*pi()*(M-1)*R*C);
Gp = tf([1], [1/(2*pi()*fp) 1]);

Gvd_dcm = Gd0_dcm*Gp;

Tu_dcm = Gvd_dcm * Gpwm;

handle = bodeplot(Tu_dcm,options);
p=getoptions(handle);    
p.Ylim{1}= [-40 50]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Ylim{2}= [-225 10]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Title.String = 'DCM Uncompensated Loop Response';
p.Title.FontWeight = 'Bold';
%p.Title.FontSize = 20;
p.PhaseMatching = 'on';
setoptions(handle,p); %update your plot


%% Regulator
subplot(2,3,[1,4])

%%Continous TF
Gc_s = tf([C1*R2 , 1], [R1*C1*C3*R2, (C3+C1) * R1, 0]); %Derived from first principles

%%Discrete TF
Gc_z = c2d(Gc_s, 0.000001, 'tustin'); %100Khz sample rate

%%plot CS and DS
handle = bodeplot(Gc_s, Gc_z,options);
title('LS2 Regulator');
p=getoptions(handle);
p.Ylim{1}= [-90 45]; %Setting the y-axis limits
p.Xlim{1}= [1 10^3]; %Setting the x-axis limits
p.Title.FontWeight = 'Bold';
setoptions(handle,p); %update your plot
legend('Continous Time', 'Discrete Time');

%% CCM compensated

subplot(2,3,3)

Tc_ccm_s = Gd0_ccm*Gvd_ccm*Gpwm*Gc_s; %% compensated Loop gain in s
%Tc_ccm_z = Gd0*Gvd_ccm*Gpwm*Gc_z; %% compensated Loop gain in z

handle = bodeplot(Tc_ccm_s, options);
p=getoptions(handle);    
p.Ylim{1}= [-40 50]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Ylim{2}= [-225 10]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Title.String = 'CCM Compensated Loop Response';
p.Title.FontWeight = 'Bold';
%p.Title.FontSize = 20;
p.PhaseMatching = 'on';
setoptions(handle,p); %update your plot
%legend('Continous Compensated Response', 'Discrete Compensated Response');

%figure(3)
%bode(Batt);
% 
% figure(4)
% subplot(1,2,1)
% handle = bodeplot(Gvd*Gdo, options);%% Erickson
% p=getoptions(handle);    
% p.Ylim{2}= [-180 180]; %Setting the y-axis limits
% p.Xlim{1}= [1 10^4]; %Setting the x-axis limits
% p.Title.String = 'Erickson';
% setoptions(handle,p); %update your plot
% 
% Gdo = Vin/D_prime^2;
% 
% wz = (D_prime^2*R)/L;
% w0 = (1/sqrt(L*C))* D_prime
% Q = w0/(1/(C*R))
% 
% GVd = tf([-1/wz 1], [(1/w0)^2 1/(w0*Q) 1])
% 
% subplot(1,2,2)
% handle = bodeplot(Gvd*Gdo, options);%% Erickson
% p=getoptions(handle);    
% p.Ylim{2}= [-180 180]; %Setting the y-axis limits
% p.Xlim{1}= [1 10^4]; %Setting the x-axis limits
% p.Title.String = 'TI';
% setoptions(handle,p); %update your plot


%% Compensated DCM
subplot(2,3,6)

Tc_dcm_s = Gvd_dcm*Gpwm*Gc_s; %% compensated Loop gain in s
%Tc_dcm_z = Gd0*Gvd_dcm*Gpwm*Gc_z; %% compensated Loop gain in z

handle = bodeplot(Tc_dcm_s,options);
p=getoptions(handle);    
p.Ylim{1}= [-40 50]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Ylim{2}= [-225 10]; %Setting the y-axis limits
p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
p.Title.String = 'DCM Compensated Loop Response';
p.Title.FontWeight = 'Bold';
%p.Title.FontSize = 20;
p.PhaseMatching = 'on';
setoptions(handle,p);
%legend('Continous Compensated Response', 'Discrete Compensated Response');

% subplot(3,2,1)
% handle = bodeplot(Batt,options);
% p=getoptions(handle);    
% 
% p.Ylim{1}= [-40 120]; %Setting the y-axis limits
% p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
% p.Ylim{2}= [-225 10]; %Setting the y-axis limits
% p.Xlim{1}= [1 10^5]; %Setting the x-axis limits
% p.Title.String = 'Battery';
% p.Title.FontWeight = 'Bold';
% p.Title.FontSize = 20;
% p.PhaseMatching = 'on';
% setoptions(handle,p); %update your plot


%% New Regulator Design
figure();

%%Desired Values
fc = 4e3 %% set the cut off frequency to x10 below the switching frequency
pm = 90;  %% this is the phase margin at the desired cutoff (5Khz) this differs from Erickson because we need to reduce pm not increase it
pmd = 60; %% we want 60

%PD reg

f0 = 1/(2*pi*sqrt(L*C)); %% this is used to calc the DC gain, but I did this manually instead, so it's not used.

fz = fc * sqrt( (1-sin(pmd - pm)) / (1+sin(pmd - pm)) ) %% from erickson
Gz = tf([1 0], [2*pi*fz]) + tf([1], [1]);
fp = fc * sqrt( (1+sin(pmd - pm)) / (1-sin(pmd - pm)) )%% from erickson
Gp = tf([1 0], [2*pi*fp]) + tf([1], [1]);

Gdc_pd = 4.47;%%(fc/f0)^2 * (1/dcgain(Tu_dcm)) * sqrt(fp/fz) This was the difference in gain from the 0dB -> gain at 5Khz   =----- I could flip this upside down too

%I = tf([1], [1 0]); %% add a pole at the origin to remove ss error

Gpd = (Gp/Gz)*Gdc_pd*50*20;  %% flipped this upside down, so that we're subtracking phase instead of adding

Tc_dcm_pd = Gpd*Tu_dcm; %% convolve in laplac

subplot(1,3,1)

%Gpd = tf([2724 2.512e7], [0.003272 350 0])

%Tc_dcm_pd = Gpd*Tu_dcm;

bode(Tu_dcm, options)
title('Uncompensated DCM Loop Response', 'fontweight', 'bold');
subplot(1,3,2)
bode(Gpd,options)
title('Proposed Compensator', 'fontweight', 'bold');
subplot(1,3,3)
bode(Tc_dcm_pd,options)
title('Compensated DCM Loop Response', 'fontweight', 'bold');


%% Comparison of old and new design + TH

figure();

%%TH PI controlle
subplot(2,3,3);
P = tf([1], [1]);
I = tf([0.5], [1 0]);
Gpi = P + I;
Tc_dcm_pi = Tu_dcm*Gpi; 
step(Tc_dcm_pi/(Tc_dcm_pi+1), 1e-2);
title('PI Controller Step Response', 'fontweight', 'bold')
subplot(2,3,6);
step(Tc_dcm_pi*Batt/((Tc_dcm_pi*Batt)+1), 1e-2);
title('PI Controller Step Response with Battery', 'fontweight', 'bold');
%%MS PD controlle
subplot(2,3,1);
step(Tc_dcm_pd/(Tc_dcm_pd+1), 1);
title('New Controller Step Response', 'fontweight', 'bold');
subplot(2,3,4); 
step((Tc_dcm_pd*Batt)/((Tc_dcm_pd*Batt)+1), 1);
title('New Controller Step Response with Battery', 'fontweight', 'bold');
%%LS controller
subplot(2,3,2); 
step(Tu_dcm*Gc_s/(Tu_dcm*Gc_s+1), 1e-2);
title('Old (LS2) Controller Step Response', 'fontweight', 'bold');
subplot(2,3,5); 
step(Tu_dcm*Gc_s*Batt/((Tu_dcm*Gc_s*Batt)+1), 1e-1)
title('Old (LS2) Controller Step Response with Battery', 'fontweight', 'bold');

%% PZ map
figure();
pzmap(Tc_dcm_pd)

%% Discretize the controller
figure();
Gpd_z = c2d(Gpd, 0.000001, 'tustin') %100Khz sample rate
bode(Gpd, Gpd_z, options);
title('Digital Vs Continous', 'fontweight', 'bold');
[n,d,t] = tfdata(Gpd_z)
Gpd_z_m = tf(n,d,t,'variable','z^-1')


%% PD 
p_want = 55
p_now = 1.5

phi = p_want - p_now;

a2 = (-sin(phi)-1 )/ ( sin(phi) -1 )

tau = 1/(1e3*2*pi()*sqrt(a2));

LEAD = tf([a2*tau 1], [tau 1])