function results = segmentParticipantClips(FV, vidInfo, participantIndex)

screen_width  = 1280;
screen_height = 1024;

participant = FV{1, participantIndex};
nMovies     = length(participant.ftimes);
viewOrder   = participant.mclipnum;
blinkData   = participant.blinks;

for m = 1:nMovies
    movieID       = viewOrder(m);
    frameTimes    = participant.ftimes{m,1};
    raw           = participant.rawtrials{m,1};
    clipStarts    = vidInfo.clipTimes{movieID}(:);
    nFrames       = vidInfo.nFrames(movieID);
    clipEnds      = [clipStarts(2:end) - 1; nFrames];
    nClips        = length(clipStarts);
    allFrameIdx   = frameTimes(:,1);
    allTimestamps = frameTimes(:,2);

    % CORRECT BLINK MAPPING: use viewing order index "m"
    movieBlinks = blinkData(blinkData.trial == m, :);

    clipData = struct();
    clipAgg  = struct();

    for c = 1:nClips
        sIdx = clipStarts(c);
        eIdx = clipEnds(c);

        sMatch = find(allFrameIdx == sIdx, 1);
        eMatch = find(allFrameIdx == eIdx, 1);
        if isempty(sMatch) || isempty(eMatch), continue; end

        t_start = allTimestamps(sMatch);
        t_end   = allTimestamps(eMatch);
        timeCol = raw{:,1};

        mask      = timeCol >= t_start & timeCol <= t_end;
        clipRaw   = raw(mask, :);
        clipAbsIdx = find(mask);

        % Blink mask
        clipBlinkMask = false(height(clipRaw),1);
        for b = 1:height(movieBlinks)
            blinkStart = movieBlinks.FullWidthIND(b,1);
            blinkEnd   = movieBlinks.FullWidthIND(b,2);
            if blinkEnd < clipAbsIdx(1) || blinkStart > clipAbsIdx(end), continue; end
            blinkStartClip = max(blinkStart, clipAbsIdx(1));
            blinkEndClip   = min(blinkEnd, clipAbsIdx(end));
            localIdx = (blinkStartClip:blinkEndClip) - clipAbsIdx(1) + 1;
            clipBlinkMask(localIdx) = true;
        end

        time     = clipRaw{:,1};
        x        = clipRaw{:,2};
        y        = clipRaw{:,3};
        pupil_mm = p_area2mm(clipRaw{:,4});

        off_screen_mask = (x < 0) | (y < 0) | ...
                          (x > screen_width) | (y > screen_height);

        % Frame mapping
        frameCol = zeros(size(time));
        for i = 1:length(time)
            [~, idx] = min(abs(double(allTimestamps) - double(time(i))));
            frameCol(i) = allFrameIdx(idx);
        end

        % RAW output
        clipTable = table(time, frameCol, x, y, pupil_mm, ...
            clipBlinkMask, off_screen_mask, ...
            'VariableNames', {'time','frame','x','y','pupil_mm','isBlink','isOffScreen'});
        clipData(c).data = clipTable;

        % AGG MASKING
        xAgg     = x;
        yAgg     = y;
        pupilAgg = pupil_mm;

        badSamples = clipBlinkMask | off_screen_mask | isnan(x) | isnan(y);

        xAgg(badSamples)     = NaN;
        yAgg(badSamples)     = NaN;
        pupilAgg(badSamples) = NaN;

        % AGGREGATE
        [frame_groups, frame_vals] = findgroups(frameCol);

        tstamp = splitapply(@median, time, frame_groups);
        gaze_x = splitapply(@(v) median(v,'omitnan'), xAgg, frame_groups);
        gaze_y = splitapply(@(v) median(v,'omitnan'), yAgg, frame_groups);
        pupil  = splitapply(@(v) median(v,'omitnan'), pupilAgg, frame_groups);

        aggTable = table(frame_vals, tstamp, gaze_x, gaze_y, pupil, ...
            'VariableNames', {'frame','time','gaze_x','gaze_y','pupil_mm'});
        clipAgg(c).aggData = aggTable;
    end

    baseDir = '/Users/jaspreetdodd/Desktop/data request';
    saveDir = fullfile(baseDir, sprintf('subj%02d', participantIndex));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end

    save(fullfile(saveDir, sprintf('subj%02d_movie%02d_raw.mat', participantIndex, movieID)), 'clipData');
    save(fullfile(saveDir, sprintf('subj%02d_movie%02d_aggData.mat', participantIndex, movieID)), 'clipAgg');

    results(participantIndex).movie(movieID).clip    = clipData;
    results(participantIndex).movie(movieID).aggData = clipAgg;
    results(participantIndex).viewingOrder = viewOrder;
end
end
