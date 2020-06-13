clear all;
    
Vin = 12;       %% Ratio
Vref = 1.8;       %% Volts
L = 150e-6;       %% 1 mH
C = 107e-6;     %% 250 uF
R = 20;          %% 30 ohms
Vout = 5.0;       %% Volts
H = Vref/Vout;  %% Ratio
D = Vout/Vin;   %% Ratio
Vm = 3.3;         %% Volts
fs = 100e3;     %% herts, switching frequency

%% construct the open loop transfer function

num = 1;
den =[L*C, L/R, 1];

f0 = 1/(2*pi*sqrt(L*C))        %% Corner/natural frequency (this must be greater than the cut off frequency)

Tua = tf(num, den) %% AC model
Tu0 = (H*Vout)/(D*Vm);  %% Open loop DC gains

Tu = Tu0*Tua;          %% Uncompensated Loop gain TF

%% Plot the uncompensated response
options = bodeoptions;
options.FreqUnits = 'Hz';
figure(1);
subplot(4,3,[1,4]);
bode(Tu, options);
title('Uncompensated Loop Gain','fontweight','bold');

%% Calculate the PD compensator
desired_phase = 85   %% degrees
fc = 5000            %% Hz

fz = fc * sqrt( (1-sind(desired_phase)) / (1+sind(desired_phase)) )
fp = fc * sqrt( (1+sind(desired_phase)) / (1-sind(desired_phase)) )

%% Gain needed to keep the current cut off frequency
Gc0 = sqrt(fz/fp)*(1/Tu0)*(fc/f0)^2;
%% Generate the PD compensator TF
Gpd_num = [1/(2*pi*fz), 1];
Gpd_den = [1/(2*pi*fp), 1];
Gpd = tf(Gpd_num, Gpd_den);

%% plot the PD bode plot
subplot(4,3,10);
bode(Gpd,options)
title('PD compensator','fontweight','bold');

%% Gain needed to keep the current cut off frequency
Gc_inf = fc/(Tu0*f0);
%% Generate the PI compensator TF
fl = 0.1*fc;            %% Set the inverted 0 to 1 tenth the cross over frequency 

Gpi_num = [1/(2*pi*fl), 1];       %% the inverted zero will be placed at one tenth the ...
Gpi_den = [1/(2*pi*fl), 0];                        %% ...switching frequency
Gpi = tf(Gpi_num, Gpi_den);

%% plot the PI bode plot
subplot(4,3,11);
bode(Gpi,options)
title('PI compensator','fontweight','bold');

%% plot the PID bode plot
subplot(4,3,12);
Gpid = Gpi*Gpd*Gc0;
bode(Gpid,options)
title('PID compensator','fontweight','bold');

%% Add the HF pole
f_gbw = 4e6;        %% gain bandwidth product of op amp
A_cl = 20*log10(Gc0*(fp/fz))  %% Closed loop gain
f_cl = f_gbw/A_cl           %% maximum frequency at which we can keep this gain

Gh_num = [1];
Gh_den = [1/(2*pi*f_cl), 1];
Gh = tf(Gh_num,Gh_den);
Gc_pid_h = Gc0*Gpd*Gpi*Gh;

subplot(4,3,8);
bode(Gc_pid_h, options);
title('PID Compensator with HF pole','fontweight','bold');

%%plot the PID compensated loop gain
subplot(4,3,[2,5]);
Gc_pid = Gc0*Gpd*Tu*Gpi;
bode(Gc_pid, options);
title('PID Compensated Loop Gain','fontweight','bold');

%% Plot the PID Compensated Loop Gain WITH high frequency pole
subplot(4,3,[3,6]);
bode(Gc_pid_h*Tu, options);
title('PID Compensated Loop Gain with HF pole','fontweight','bold');
hold off;
%% Calculate Component Values:
R2 = 100e3;         %% This is set

R1 = R2/Gc0
R2
C2 = (2*pi*R2*fl)^-1;
C1 = (2*pi*R1*fz)^-1;
C1
C2
R3 = (2*pi*C1*fp)^-1
Rs1 = R3/H
Rs2 = (H*Rs1)/(1-H)

%% Digital Regulator Stuff

figure(2)

t_cntl = 3e-6; %%input('Control time: ');
t_delay = D*(1/fs);

t_total = t_cntl + t_delay;

Tud = Tu * tf(1,1,'InputDelay',t_total) %% delayed input, can now be discretized

subplot(1,2,1)

Tud = c2d(Tud,fs^-1, 'method', 'Tustin');
Gpid_d = c2d(Gpid, fs^-1, 'method', 'Tustin');

bode(Tud, Tu, options)
legend('Discrete', 'Continous')

subplot(1,2,2)
bode(Tud*Gpid_d, Tu*Gpid, options)
legend('Discrete', 'Continous')

%[n d t] = tfdata(Gpid_d);

%Gpid_d_m = tf(n,d,t, 'variable', 'z^-1')
Gpid_d.variable='z^-1'
