%% ================================================================
%  retinal_extract_all_movies.m
%  Build gaze-centered 2x retinal windows and 1°x1° patch luminance
%  for ALL movies for one subject.
%
%  OUTPUT: subj##_retinalcoor.mat  (matrix + global luminance + metadata)
%
%  Author: Jaspreet Dodd (24cb7 @ queensu.ca)
%  Last updated: 12/01/25
%  ================================================================

%% USER INPUTS
subjRange   = 7:51;
lumBasePath = '/Users/jaspreetdodd/Desktop/luminance_frames'; % change this to where luminance frames live

movies = 1:10;  % movie01 ... movie10

% Visual params
degPerPatch = 1;         % patch = 1° × 1°
borderScale = 2.0;       % retinal window = 2× screen size
bgLum       = 0.12;      % This is the photometer measured average luminance around the screen

for s = subjRange
    subjID   = sprintf('subj%02d', s);
    baseDir  = '/Users/jaspreetdodd/Desktop/local analysis';   % <- NEW BASE
    subjDir  = fullfile(baseDir, subjID);                     % <- subjXX folder

    if ~isfolder(subjDir)
        error('Subject folder not found: %s', subjDir);
    end

    PPD         = FV{1,s}.res.PPD;

%% Determine grid size
% We load the first available luminance frame from movie01
testLum = load(fullfile(lumBasePath,'movie01','lum_frame0001.mat'));
[LH, LW] = size(testLum.lum_map);      % LH = 1024, LW = 1280

pxPerPatch = round(PPD * degPerPatch); % ~41 px per 1° patch

% Target retinal window in pixels (2× size)
Wbg_target = borderScale * LW;
Hbg_target = borderScale * LH;

% Number of whole patches that fit in each dimension
nX = floor(Wbg_target / pxPerPatch);   % horizontal patches
nY = floor(Hbg_target / pxPerPatch);   % vertical patches

% Final enforced retinal window dimensions
Wbg = nX * pxPerPatch;
Hbg = nY * pxPerPatch;

halfW = floor(Wbg/2);
halfH = floor(Hbg/2);

fprintf('Retinal grid: %d × %d patches (nX × nY)\n', nX, nY);

%% Preallocate structure
% We will store each movie as a separate cell
matrix = cell(1, numel(movies));           % each cell: [frames × nPatches]
global_screen = cell(1, numel(movies));    % each cell: [frames × 1]
frame_numbers = cell(1, numel(movies));    % each cell: [frames × 1]
gaze_px = cell(1, numel(movies));          % each cell: [frames × 2]

nPatches = nX * nY;

%% Main loop over movies
for m = movies
    movieID = sprintf('movie%02d', m);

    fprintf('\n=== Processing %s ===\n', movieID);

    % Load aggregated data
    aggFile = sprintf('%s_%s_aggData.mat', subjID, movieID);
    aggPath = fullfile(subjDir, aggFile);

    if ~isfile(aggPath)
        error('AggData file not found: %s', aggPath);
    end

    load(aggPath, 'clipAgg');

    % Combine all clips into one table for this movie
    T = [];
    for c = 1:numel(clipAgg)
        T = [T ; clipAgg(c).aggData]; 
    end

    N = height(T);
    Xmov   = nan(N, nPatches);
    Gmov   = nan(N,1);
    Fmov   = nan(N,1);
    GMmov  = nan(N,2);
    Pup    = nan(N,1);

    % Pad size
    padX = halfW;
    padY = halfH;

    % Loop frames
    for f = 1:N
        gx = T.gaze_x(f);
        gy = T.gaze_y(f);
        pup = T.pupil_mm(f);

        Fmov(f) = T.frame(f);
        GMmov(f,:) = [gx gy];
        Pup(f,1) = pup;

        % If gaze is NaN then store NaNs for the whole row (including pupil)
        if isnan(gx) || isnan(gy)
            continue;
        end

        % Load luminance frame
        lumFile = fullfile(lumBasePath, movieID, ...
            sprintf('lum_frame%04d.mat', T.frame(f)));
        S = load(lumFile);
        L = S.lum_map;

        % Global screen luminance (original 1280×1024)
        Gmov(f) = mean(L(:),'omitnan');

        % Pad the frame with background luminance
        Lpad = padarray(L, [padY padX], bgLum, 'both');

        % Gaze coords in padded space
        x0 = round(gx) + padX;
        y0 = round(gy) + padY;

        % Crop retinal window of size Hbg × Wbg
        x1 = x0 - halfW + 1;
        x2 = x0 + halfW;
        y1 = y0 - halfH + 1;
        y2 = y0 + halfH;

        % Clamp to padded boundaries
        x1 = max(1, x1);  y1 = max(1, y1);
        x2 = min(size(Lpad,2), x2);
        y2 = min(size(Lpad,1), y2);

        R = Lpad(y1:y2, x1:x2);

        % Enforce exact dims (may crop/pad by ±1)
        R = R(1:min(end,Hbg), 1:min(end,Wbg));
        if size(R,1) < Hbg
            R(end+1:Hbg, :) = bgLum;
        end
        if size(R,2) < Wbg
            R(:, end+1:Wbg) = bgLum;
        end

        % Reshape into patches: [nY x nX]
        R4 = reshape(R, pxPerPatch, nY, pxPerPatch, nX);
        R4 = permute(R4, [2 4 1 3]);   % [nY x nX x px x px]
        patchMeans = mean(mean(R4,4,'omitnan'),3,'omitnan');

        % Flatten into vector with Gunnar's indexing
        vec = nan(1, nPatches);
        for i = 1:nX
            for j = 1:nY
                k = (i-1)*nY + j;
                vec(k) = patchMeans(j, i);
            end
        end

        Xmov(f,:) = vec;
    end

    % Store into output
    matrix{m} = Xmov;
    global_screen{m} = Gmov;
    frame_numbers{m} = Fmov;
    gaze_px{m} = GMmov;
    pupil_mm{m} = Pup;

end

%% Metadata
metadata = struct();
metadata.nX = nX;
metadata.nY = nY;
metadata.nPatches = nPatches;
metadata.pxPerPatch = pxPerPatch;
metadata.degPerPatch = degPerPatch;
metadata.PPD = PPD;
metadata.Wbg = Wbg;
metadata.Hbg = Hbg;
metadata.borderScale = borderScale;
metadata.indexFcn = 'k = (i-1)*nY + j;  i=1..nX, j=1..nY';
metadata.fields = { ...
    'matrix: cell{movie}{frame,k} → patch means', ...
    'global_screen: cell{movie}{frame} → avg luminance', ...
    'frame_numbers: cell{movie}{frame} → movie frame index', ...
    'gaze_px: cell{movie}{frame,2} → gaze [x,y] in px', ...
    'pupil_mm: cell{movie}{frame} → pupil diameter mm' ...
};
%% Save
saveFile = sprintf('%s_retinalcoor.mat', subjID);
save(saveFile, 'matrix', 'global_screen', 'frame_numbers', ...
               'gaze_px', 'pupil_mm', 'metadata', '-v7.3');

fprintf('\nSaved: %s\n', saveFile);

end
