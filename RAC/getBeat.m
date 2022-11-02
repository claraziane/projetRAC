function[beatFreq, BPM, IOI, beatOnset] = getBeat(Audio, Freq, bpmInitial)

% Centering the signal around 0
%Audio = Audio-mean(Audio);
audioTemp = Audio;
windowLength = Freq/4;
for iWindow = 1:length(Audio)-windowLength
    audioTemp(iWindow:iWindow+windowLength) = audioTemp(iWindow:iWindow+windowLength)-mean(audioTemp(iWindow:iWindow+windowLength));
end
Audio = mean(horzcat(Audio,audioTemp),2);

% Inverting the signal
Audio = -Audio;

% Plot audio signal
figure;
plot(Audio); hold on;

% Low-pass filter audio signal at 40 Hz to get signal envelop
[f,e] = butter(2,2*40/Freq);
audioFilt = filtfilt(f,e,abs(Audio)); plot(audioFilt)

% Find envelop peaks
minDistance = (Freq/(bpmInitial/60))-100;
[pksFilt, locsFilt] = findpeaks(audioFilt, 'MinPeakDistance', minDistance); plot(locsFilt, audioFilt(locsFilt), 'b*')

% Find peaks corresponding to beat onsets
minPeak = 0;
[pks,locs] = findpeaks(Audio, 'MinPeakHeight', minPeak);
%peakDistance = diff(locs);

beatOnset = []; beatValue = [];
nPeaks = 6; %Number of peaks to include before trigger
for iPksFilt = 2:length(pksFilt)

    if iPksFilt < 44
        [M, I] = min(abs(locs-locsFilt(iPksFilt)));
        tempAudio = Audio(locs(I-nPeaks:I));
        tempFrames = locs(I-nPeaks:I); plot(tempFrames, tempAudio, 'k*')
        [minAudio, minIndex] = min(tempAudio);
        audioRound = round(tempAudio);

        for i= 1:length(audioRound)-1
            if tempAudio(i) < 1.2 && tempAudio(i+1) > 1.3 %% audioRound(i) == 0 && audioRound(i+1) ~= 0 % && abs(tempAudio(i+1) - tempAudio(i)) > 0.5
                tempAudio(1:i-1) = [];
                tempFrames(1:i-1) = [];
                break;
            end
        end

    elseif iPksFilt > 43 && iPksFilt < 140
        nPeaks = 10;
        [M, I] = min(abs(locs-locsFilt(iPksFilt)));
        tempAudio = Audio(locs(I-nPeaks:I));
        tempFrames = locs(I-nPeaks:I); plot(tempFrames, tempAudio, 'k*')
        [minAudio, minIndex] = min(tempAudio);
        audioRound = round(tempAudio);

        for i= 1:length(audioRound)-1
            if tempAudio(i) < 1.2 && tempAudio(i+1) > 1.3 %% audioRound(i) == 0 && audioRound(i+1) ~= 0 % && abs(tempAudio(i+1) - tempAudio(i)) > 0.5
                tempAudio(1:i-1) = [];
                tempFrames(1:i-1) = [];
                break;
            end
        end

    elseif iPksFilt > 139 && iPksFilt < 164
        nPeaks = 5;
        [M, I] = min(abs(locs-locsFilt(iPksFilt)));
        tempAudio = Audio(locs(I-nPeaks:I));
        tempFrames = locs(I-nPeaks:I); plot(tempFrames, tempAudio, 'k*')
        [minAudio, minIndex] = min(tempAudio);
        audioRound = round(tempAudio);

        for i= 1:length(audioRound)-1
            if tempAudio(i+1) > 1.1 %% audioRound(i) == 0 && audioRound(i+1) ~= 0 % && abs(tempAudio(i+1) - tempAudio(i)) > 0.5
                tempAudio(1:i-1) = [];
                tempFrames(1:i-1) = [];
                break;
            end
        end

    elseif iPksFilt == 164
        nPeaks = 5;
        [M, I] = min(abs(locs-locsFilt(iPksFilt)));
        tempAudio = Audio(locs(I-nPeaks:I));
        tempFrames = locs(I-nPeaks:I); plot(tempFrames, tempAudio, 'k*')
        tempFrames(2) = 152751;

    elseif iPksFilt > 164 %&& iPksFilt < 166
        nPeaks = 6;
        [M, I] = min(abs(locs-locsFilt(iPksFilt)));
        tempAudio = Audio(locs(I-nPeaks:I));
        tempFrames = locs(I-nPeaks:I); plot(tempFrames, tempAudio, 'k*')
        [minAudio, minIndex] = min(tempAudio);
        audioRound = round(tempAudio);

        for i= 1:length(audioRound)-1
            if tempAudio(i+1) > 0.6 && tempAudio(i+2) > 0.6%% audioRound(i) == 0 && audioRound(i+1) ~= 0 % && abs(tempAudio(i+1) - tempAudio(i)) > 0.5
                tempAudio(1:i-1) = [];
                tempFrames(1:i-1) = [];
                break;
            end
        end
      
%         elseif iPksFilt > 165 %&& iPksFilt < 166
%         nPeaks = 6;
%         [M, I] = min(abs(locs-locsFilt(iPksFilt)));
%         tempAudio = Audio(locs(I-nPeaks:I));
%         tempFrames = locs(I-nPeaks:I); plot(tempFrames, tempAudio, 'k*')
%         [minAudio, minIndex] = min(tempAudio);
%         audioRound = round(tempAudio);
% 
%         for i= 1:length(audioRound)-1
%             if tempAudio(i+1) > 1 %% audioRound(i) == 0 && audioRound(i+1) ~= 0 % && abs(tempAudio(i+1) - tempAudio(i)) > 0.5
%                 tempAudio(1:i-1) = [];
%                 tempFrames(1:i-1) = [];
%                 break;
%             end
%         end
    end

    beatIndex = tempFrames(2);
    while Audio(beatIndex) > Audio(beatIndex-1)
        beatIndex = beatIndex-1;
    end
    beatOnset = [beatOnset; beatIndex];
    beatValue = [beatValue; Audio(beatIndex)];

end

% Plot beat onsets as red stars on open figure
plot(beatOnset, beatValue, 'r*')

% Extract metronome IOI, frequency & BPM
nBeats   = length(beatOnset);
for iBeat = 1:length(beatOnset)-1
    IOI(iBeat) = beatOnset(iBeat+1)-beatOnset(iBeat);
end

% Make sure no beat is missing
IOIDistance = diff(IOI);
for iIOI = 1:length(IOIDistance)
    if abs(IOIDistance(iIOI)) > 30
        warning(' !!! Seems like at least one beat is missing !!!' );
    end
end

beatFreq = nBeats / ((beatOnset(end) - beatOnset(1))/Freq);
BPM      = beatFreq * 60;

end
