function [r] = pupil_model(Is, Ip, dt)
%
% Fan, Yao (2011). Modeling Transient Pupillary Light Reflex Induced by a
% Short Light Flash. IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING 58(1),  36-42
%
% Is: sympathetic input
% Ip: para-sympathetic input
% dt: time step

%% parameters
% human
Kc = 0.06; % elastic constant of constrictor
Kd = 1.25; % elastic constant of dilator
L0c = 0.9; % radius of constrictor at restdt
L0d = 3.6; % radius of dilator at rest
P0 = -0.76; % static force at rest
D = 3.78; % viscosity
taup = 0.26; % parasympathetic system delay (s)
taus = 0.64; % sympathetic system delay (s)
% taup = 0.2;

% time steps
% dt = 0.01; % time step (s)
% t = 0:dt:10; % time line
t = (0:(length(Is)-1))*dt;

% pupil diameter range and limits
pmin = 1; % pupil min diameter
pmax = 7; % pupil max diameter
k = 5; % softclip hardness

%% retinal input (light stimulation)
% stimL = 5;
% fs = zeros(size(t)); % sympathetic input
% fs(round((1/dt:(1+stimL)/dt)+taus/dt)) = 2;
% fp = zeros(size(t)); % parasympathetic input
% % fp(round((1/dt:(1+stimL)/dt)+taup/dt)) = 3;
fp = Ip; fs = Is;
fp1 = fp;
fs1 = fs;

%% initial conditions
r(1) = 2.724;
% r(1) = 2.674;
r0(1) = r(1);
v(1) = 0;

%% oscillator analysis (strongly damped)
% solution (homogeneous): r(t)=r1*exp(lam1*t)+r2*exp(lam2*t)
% bet=D/2;
% om2=2*Kc*(L0c-r(1))+2*Kd*(L0d-r(1));
% lam1=-bet+sqrt(bet^2-om2);
% lam2=-bet-sqrt(bet^2-om2);

%% simulation
for i = 1:length(t)-1,
    fp1(min(i+round(taup/dt),length(t))) = fp(min(i+round(taup/dt),length(t)))*log(r(i)^2)/log(r(1)^2);
    fs1(min(i+round(taus/dt),length(t))) = fs(min(i+round(taus/dt),length(t)))*log(r(1)^2)/log(r(i)^2);
    a(i) = -Kc*(-L0c+r(i))^2 + Kd*(L0d-r(i))^2 - D*v(i) - fp1(i) + fs1(i) + P0;
    v(i+1) = v(i) + dt*a(i);
    r(i+1) = softclip(r(i) + dt*v(i), k, [pmin pmax]);
%     r0(i+1) = r0(i) + dt*v(i);
end

%% plot results
% figure; hold on;
% plot(t,r,'r-')
% % plot(t,r0,'g:')
% plot([1 1],[r(1)-.1 r(1)+.1],'k:')
% plot([1+stimL 1+stimL],[r(1)-.1 r(1)+.1],'k:')
% % set(gca, 'yscale', 'log')
% axis([1 6 2 3.5])
% xlabel('Time (s)')
% ylabel('Pupil radius (mm)')
% figure
% plot(r,v)
% figure;
% plot(r(1:end-1),a)
