%% This script loads videos frames and computes luminance per pixel 
% 
clc;
clear all;

% 1. Load video reader to extract frames 
videoPath = '/Users/jaspreetdodd/Desktop/analysis/FV_m4v/mtv3clip10_fixed.mp4';  % change this to movie path
vid = VideoReader(videoPath);      % will be used to load frames by index

% Screen geometry 
screen_width_px = 1280;
screen_height_px = 1024;

nFrames = floor(vid.Duration * vid.FrameRate);  % Total frames

% luminance = zeros(screen_height_px, screen_width_px, nFrames, 'single');

outputDir = '/Users/jaspreetdodd/Desktop/analysis/luminance_frames/movie10'; % change this per movie - should match the clip you are grabbing
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for f = 1:nFrames
    frame_rgb = read(vid, f);
    rgb_vals = double(reshape(frame_rgb, [], 3)); % Flatten to Nx3
    rgb_vals = round(rgb_vals); % Clamp RGB for LUT
    rgb_vals = min(max(rgb_vals, 1), 255); % Clamp RGB for LUT

    lum_rgb = rgb2lum(rgb_vals, 1); % Apply gamma LUT
    pixel_lum = sum(lum_rgb, 2);  % Combine R+G+B channels
    lum_map = reshape(pixel_lum, screen_height_px, screen_width_px);

    %Save as individual .mat
    save(fullfile(outputDir, sprintf('lum_frame%04d.mat', f)), 'lum_map', '-v7.3');

    if mod(f, 50) == 0
        fprintf('Saved frame %d / %d\n', f, nFrames);
    end
end
