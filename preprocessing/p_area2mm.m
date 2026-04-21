function [pupil_in_mm] = p_area2mm(pupil_area, viewing_distance)
%[pupil_in_mm] = p_area2mm(pupil_area)
%   Converts Eye link pupil area (pixels) to pupil diameter in mm
%   from artificual eye measurements from Don:
%  mm [2 4 6 8 10], area [402 1425 3049 5130 8240]
% note the smallest measure seemed a bit innaccurate so I left it out of
% the calc
%
% B. White, Dec 2024

% if we develop the function to compute based on different view distances
if nargin == 1
    viewing_distance = 600; % mm  
end
% Don's calcs
p_mm = [2 4 6 8 10];
p_area = [402 1425 3049 5130 8240];

% p = polyfit(p_mm,sqrt(p_area),1); % linear seems right here when taking sqrt(area)
%
% x = 0.5 : 0.1 : 10;
% y = p(1) .* x + p(2)

PP = spline(p_area, p_mm);% fit spline
pupil_in_mm = ppval(PP,[pupil_area]); % convert pixels area to mm

end