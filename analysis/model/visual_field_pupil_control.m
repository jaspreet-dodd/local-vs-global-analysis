% visual_field_pupil_control.m

%% let's start
warning off
clear
addpath '/Users/jaspreetdodd/Desktop/local analysis/code';
addpath '/Users/jaspreetdodd/Desktop/local analysis/data';

%% Load data
% load('subj01_retinalcoor2.mat')
% load('subj02_retinalcoor.mat')
% load('subj03_retinalcoor.mat')
% load('subj04_retinalcoor.mat')
% load('subj05_retinalcoor.mat')
load('subj06_retinalcoor.mat')
% options
op = 2; % 1: pupil-lumincance correlation
        % 2: pupil-simulated pupil (based on luminance) correlation

%% Set up data
if op == 1,
    d = 20; % frame delay (multiple of 16.7ms)
    ig = 120; % number of frames at beginning of each clip to ignore
elseif op == 2,
    d = 0; % frame delay (multiple of 16.7ms)
    ig = 120; % number of frames at beginning of each clip to ignore
end

X = []; Y = [];
for i = 1:length(matrix),
    X = [X; [ones(size(global_screen{i}(ig+1:end-d))) global_screen{i}(ig+1:end-d) matrix{i}(ig+1:end-d,:)]];
    % without global luminance
    % X = [X; [ones(size(global_screen{i}(ig+1:end-d))) matrix{i}(ig+1:end-d,:)]];
    % without offset or global luminance
    % X = [X; [matrix{i}(ig+1:end-d,:)]];

    % interpolate NaNs
    y = pupil_mm{i};
    nanx = isnan(y);
    t = 1:numel(y);
    y(nanx) = interp1(t(~nanx), y(~nanx), t(nanx), 'linear', 'extrap');

    % detrend
    Yd = detrend(y); % detrend pupil

    % filter
    ff = normpdf(-10:10,0,.1); ff = ff./sum(ff);
    Yd = conv(Yd, ff,'same');
    Y = [Y; Yd((ig+1+d):end)];
end

%% Convert luminance to pupil size
if op == 2,
    dt = 0.0167; % time step (s) - Gunnar put this down
    % taup = 0.26; % parasympathetic system delay (s)
    taus = 0.64; % sympathetic system delay (s)
    taup = 0.15;

    P = ones(size(X));
    arrousal = zeros(size(X,1),1);
    SCinput = zeros(size(X,1),1);
    wb = waitbar(0,'Getting started...');
    t = 1:numel(X(:,1));
    for i = 2:size(X,2),
        nanx = isnan(X(:,i));
        X(nanx,i) = interp1(t(~nanx), X(~nanx,i), t(nanx), 'linear', 'extrap');
        [fs, fp] = pupil_create_command(X(:,i)', arrousal', SCinput', taup, taus, dt);
        P(:,i) = pupil_model(fs, fp, dt);
        waitbar(i/size(X,2),wb,[num2str(100*i/size(X,2),'%1.1f')])
    end
    close(wb)
    % save('luminance2pupil.mat','P','Y')
end

%% Cross-validate regularlized linear regression model
Lambda = logspace(-5,5,30);
if op == 1,
    CVMdl = fitrlinear(X,Y,'ObservationsIn','rows','KFold',5,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
    numCLModels = numel(CVMdl.Trained)
elseif op == 2,
    CVMdl = fitrlinear(P,Y,'ObservationsIn','rows','KFold',5,'Lambda',Lambda,...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
    numCLModels = numel(CVMdl.Trained)
end
mse = kfoldLoss(CVMdl);

% ind = 7;
% find best model
ind = find(mse == min(mse));
if op == 1,
    Mdl = fitrlinear(X,Y,'ObservationsIn','rows','Lambda',Lambda(ind),...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
elseif op == 2,
    Mdl = fitrlinear(P,Y,'ObservationsIn','rows','Lambda',Lambda(ind),...
        'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
end

% compute R^2
if op == 1,
    Yp = predict(Mdl,X);
elseif op == 2,
    Yp = predict(Mdl,P);
end
SSR = nansum((Yp - Y).^2);
TSS = nansum(((Y - nanmean(Y)).^2));
Rsquared = 1 - SSR/TSS;

%% Plot results
% % optimal regularization hyperparameter plot
% figure
% plot(log10(Lambda),log10(mse)); 
% ylabel('log_{10} MSE')
% xlabel('log_{10} Lambda')

% plot best model weights
figure
imagesc(reshape(Mdl.Beta(3:end),49,62))
c=colorbar;
caxis(max(abs(c.Limits))*[-1 1])
title(['R2 = ' num2str(Rsquared)])

%% debug pupil model parameters...
figure;
subplot(2,2,1)
plotyy(1:size(P,1),Y,1:size(P,1),P(:,2))
subplot(2,2,3)
plotyy(1:size(P,1),X(:,2),1:size(P,1),P(:,2))
subplot(2,2,2)
plot(Y,P(:,2),'.')
[b,bint,r,rint,stats] = regress(Y,[ones(size(P(:,2))) P(:,2)]);
title(['R^2 = ' num2str(stats(1))])
% cross-correlation between pupil and pupil prediction based on global
% luminance
[c,lags] = xcorr(P(:,2),Y);
subplot(2,2,4); plot(lags,c)
mi = find(c == max(c)); 
title(['optimal lag = ' num2str(lags(mi)) ' frames'])