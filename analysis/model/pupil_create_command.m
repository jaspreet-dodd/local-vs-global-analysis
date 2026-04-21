function [fs, fp] = pupil_create_command(luminance, arrousal, SCinput, taup, taus, dt)
%
% translates luminance, arrounsal and SC activiy into neural commands for
% the sympathetic (fs) and the parasympathetic (fp) pathways

a = .04;
b = 3.5;
dl = [diff(luminance) 0];

t = -1:0.01:1;
% y = exp(-t.^2/(0.2/(dt*1000))^2); y = y./sum(y);
y = exp(-t.^2/(0.1)^2); y = y./sum(y);
dls = conv(dl, y, 'same');

fp = a*luminance + b*dls;
fp = [zeros(1,round(taup/dt)) fp(1:end-round(taup/dt))];
fp = max(fp,0);

fs = zeros(size(fp));