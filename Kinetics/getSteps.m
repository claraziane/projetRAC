%% This function extracts participants' step onsets from left and right kinetic data
%
% Input variables:
% -kineticR: z force time-series recorded for the right foot (vector of length = total number of frames)
% -kineticL: z force time-series recorded for the left foot (vector of length = total number of frames)
% -Freq:     acquisition frequency (single value variable)
%
% Output variables:
% -stepOnsetR:   right step onset time values expressed in frames (vector of length = number of right steps)
% -stepOnsetL:   left step onset time values expressed in frames (vector of length = number of left steps)
% -stepOnsetAll: left and right step onset time values expressed in frames (vector of length = total number of steps)
%
% C. Ziane

function[stepOnsetR, stepOnsetL, stepOnsetAll] = getSteps(kineticR, kineticL,Freq)

%% Right foot
figure; plot(kineticR); title('Right Foot'); hold on;

% Low-pass filter audio signal at 5 Hz to get signal envelop
[f,e] = butter(2,2*5/Freq);
kineticFiltR = filtfilt(f,e,abs(kineticR));
% plot(kineticFiltR)

% Find envelop peaks
peakThreshold = 3;
[pksFilt, locsFilt] = findpeaks(kineticFiltR);

% Find first stepOnset and remove peaks before first stepOnset
[minPksFilt, minIndexPksFilt] = min(pksFilt(1:3));
if minIndexPksFilt > 1
    locsFilt(1:minIndexPksFilt-1) = [];
    pksFilt(1:minIndexPksFilt-1) = [];
end

% Only keep one peak per step
pksFiltTemp = pksFilt;
pksFiltTemp(pksFilt < peakThreshold) = 0;
pksFiltTemp(pksFilt > peakThreshold) = 1;
pks2Keep = [];

% if one zero value is missing
pksSingle = [];
for iPksFilt = 1:length(pksFilt)-3
    if mean(pksFiltTemp(iPksFilt:iPksFilt+3)) == 1
        pksFiltTemp(iPksFilt+3:end+1) = pksFiltTemp(iPksFilt+2:end);
        pksFiltTemp(iPksFilt+2) = 0;

        pksFilt(iPksFilt+3:end+1) = pksFilt(iPksFilt+2:end);
        pksFilt(iPksFilt+2) = 0;

        locsFilt(iPksFilt+3:end+1) = locsFilt(iPksFilt+2:end);
        locsFilt(iPksFilt+2) = locsFilt(iPksFilt+2)-1;
    elseif pksFiltTemp(iPksFilt) == 1 && pksFiltTemp(iPksFilt-1) == 0 && pksFiltTemp(iPksFilt+1) == 0
        pksSingle = [pksSingle; locsFilt(iPksFilt)];
    end
end

for iPksFilt = 1:length(pksFilt)
    if pksFiltTemp(iPksFilt) == 0 && iPksFilt ~= length(pksFilt)
        pks2Keep = [pks2Keep; iPksFilt+1];
    end
end
locsFilt = locsFilt(pks2Keep);
pksFilt = pksFilt(pks2Keep);

% Remove peaks below peakThreshold
locsFilt(pksFilt < peakThreshold) = [];
pksFilt(pksFilt < peakThreshold) = [];
% plot(locsFilt, kineticFiltR(locsFilt), 'b*')

% Find peaks corresponding to beat onsets
minPeak = -0.2;
[pks,locs] = findpeaks(kineticR, 'MinPeakHeight', minPeak);

stepOnsetR = []; stepValueR = [];
singleIndex = 1;
for iPksFilt = 1:length(pksFilt)
    nPeaks = 200; %Number of peaks to include before trigger
    [M, I] = min(abs(locs-locsFilt(iPksFilt)));
    if iPksFilt == 1 && I < nPeaks
        nPeaks = I-1;
    elseif ~isempty(pksSingle) && singleIndex <= length(pksSingle) && locsFilt(iPksFilt) == pksSingle(singleIndex)
        nPeaks = nPeaks*2;
        singleIndex = singleIndex+1;
    end
    tempKineticR = kineticR(locs(I-nPeaks:I));
    tempFrames = locs(I-nPeaks:I);
    kineticRound = round(tempKineticR,1);
    if kineticRound(1) > 0
        for iKineticRound = 1:length(kineticRound)
            if kineticRound(iKineticRound) <= 0 && kineticRound(iKineticRound+1) <= 0
                kineticRoundIndex = iKineticRound-1;
                break;
            end
        end
        kineticRound(1:kineticRoundIndex) = [];
        tempFrames(1:kineticRoundIndex)   = [];
        tempKineticR(1:kineticRoundIndex) = [];       
    end
    kineticRound(kineticRound<=0) = 0;
    kineticRound(kineticRound>0) = 1;
    %     plot(tempFrames, tempKineticR, 'k*')

    for i = 1:length(kineticRound)-1
        if kineticRound(i) == 0 && mean(kineticRound(i+1:end)) == 1
            tempFrames(1:i-1) = [];
            break;
        end
    end
    if mean(kineticRound) == 1
        [maxTempKineticR, indexTempKineticR] = max(abs(diff(tempKineticR)));
        tempFrames(2) = tempFrames(indexTempKineticR+1);
    end

    kineticIndex = tempFrames(2);
    while kineticR(kineticIndex) >= kineticR(kineticIndex-1)
        kineticIndex = kineticIndex-1;
    end
    stepOnsetR = [stepOnsetR; kineticIndex];
    stepValueR = [stepValueR; kineticR(kineticIndex)];

end
plot(stepOnsetR, stepValueR, 'r*')

clear pks locs pksFilt locsFilt pksSingle

%% Left foot
figure; plot(kineticL); title('Left Foot'); hold on;

% Low-pass filter audio signal at 5 Hz to get signal envelop
[h,g] = butter(2,2*5/Freq);
kineticFiltL = filtfilt(h,g,abs(kineticL));
% plot(kineticFiltL)

% Find envelop peaks
peakThreshold = 3;
[pksFilt, locsFilt] = findpeaks(kineticFiltL);

% Find first stepOnset and remove peaks before first stepOnset
[minPksFilt, minIndexPksFilt] = min(pksFilt(1:3));
if minIndexPksFilt > 1
    locsFilt(1:minIndexPksFilt-1) = [];
    pksFilt(1:minIndexPksFilt-1) = [];
end

% Only keep one peak per step
pksFiltTemp = pksFilt;
pksFiltTemp(pksFilt < peakThreshold) = 0;
pksFiltTemp(pksFilt > peakThreshold) = 1;
pks2Keep = [];

% if one zero value is missing
pksSingle = [];
for iPksFilt = 2:length(pksFilt)-3
    if mean(pksFiltTemp(iPksFilt:iPksFilt+3)) == 1
        pksFiltTemp(iPksFilt+3:end+1) = pksFiltTemp(iPksFilt+2:end);
        pksFiltTemp(iPksFilt+2) = 0;

        pksFilt(iPksFilt+3:end+1) = pksFilt(iPksFilt+2:end);
        pksFilt(iPksFilt+2) = 0;

        locsFilt(iPksFilt+3:end+1) = locsFilt(iPksFilt+2:end);
        locsFilt(iPksFilt+2) = locsFilt(iPksFilt+2)-1;
    elseif pksFiltTemp(iPksFilt) == 1 && pksFiltTemp(iPksFilt-1) == 0 && pksFiltTemp(iPksFilt+1) == 0
        pksSingle = [pksSingle; locsFilt(iPksFilt)];
    end
end

for iPksFilt = 1:length(pksFilt)
    if pksFiltTemp(iPksFilt) == 0 && iPksFilt ~= length(pksFilt)
        pks2Keep = [pks2Keep; iPksFilt+1];
    end
end
locsFilt = locsFilt(pks2Keep);
pksFilt = pksFilt(pks2Keep);

% Remove peaks below peakThreshol
locsFilt(pksFilt < peakThreshold) = [];
pksFilt(pksFilt < peakThreshold) = [];
% plot(locsFilt, kineticFiltL(locsFilt), 'b*')

% Find peaks corresponding to beat onsets
minPeak = -0.2;
[pks,locs] = findpeaks(kineticL, 'MinPeakHeight', minPeak);

stepOnsetL = []; stepValueL = [];
singleIndex = 1;
for iPksFilt = 1:length(pksFilt)
    nPeaks = 200; %Number of peaks to include before trigger
    [M, I] = min(abs(locs-locsFilt(iPksFilt)));
    if iPksFilt == 1 && I < nPeaks
        nPeaks = I-1;
    elseif ~isempty(pksSingle) && singleIndex <= length(pksSingle) && locsFilt(iPksFilt) == pksSingle(singleIndex)
        nPeaks = nPeaks*2;
        singleIndex = singleIndex+1;
    end
    tempKineticL = kineticL(locs(I-nPeaks:I));
    tempFrames = locs(I-nPeaks:I);
    kineticRound = round(tempKineticL,1);
    if kineticRound(1) > 0
        for iKineticRound = 1:length(kineticRound)
            if kineticRound(iKineticRound) <= 0 && kineticRound(iKineticRound+1) <= 0
                kineticRoundIndex = iKineticRound-1;
                break;
            end
        end
        kineticRound(1:kineticRoundIndex) = [];
        tempFrames(1:kineticRoundIndex)   = [];
        tempKineticL(1:kineticRoundIndex)   = [];
    end
    kineticRound(kineticRound<=0) = 0;
    kineticRound(kineticRound>0) = 1;
%     plot(tempFrames, tempKineticL, 'k*')

    for i = 1:length(kineticRound)-1
        if kineticRound(i) == 0 && mean(kineticRound(i+1:end)) == 1
            tempFrames(1:i-1) = [];
            break;
        end
    end
    if mean(kineticRound) == 1
        [maxTempKineticL, indexTempKineticL] = max(abs(diff(tempKineticL)));
        tempFrames(2) = tempFrames(indexTempKineticL+1);
    end

    kineticIndex = tempFrames(2);
    while kineticL(kineticIndex) >= kineticL(kineticIndex-1)
        kineticIndex = kineticIndex-1;
    end
    stepOnsetL = [stepOnsetL; kineticIndex];
    stepValueL = [stepValueL; kineticL(kineticIndex)];

end
plot(stepOnsetL, stepValueL, 'r*')

% Make sure you only have one value per step and that no steps are missing
stepOnsetDiffR = diff(stepOnsetR);
for iStepR = 1:length(stepOnsetDiffR)
    if stepOnsetDiffR(iStepR) > mean(stepOnsetDiffR) + 250
        warning(' !!! Seems like at least one right step is missing !!!' );
    elseif stepOnsetDiffR(iStepR) < mean(stepOnsetDiffR) - 250
        warning(' !!! Seems like there are too many right steps !!!' );
    end
end

stepOnsetDiffL = diff(stepOnsetL);
for iStepL = 1:length(stepOnsetDiffL)
    if stepOnsetDiffL(iStepL) > mean(stepOnsetDiffL) + 250
        warning(' !!! Seems like at least one left step is missing !!!' );
    elseif stepOnsetDiffL(iStepL) < mean(stepOnsetDiffL) - 250
        warning(' !!! Seems like there are too many left steps !!!' );
    end
end

%% Both feet
stepOnsetAll = vertcat(stepOnsetR, stepOnsetL);
stepOnsetAll = sort(stepOnsetAll);

end