 clear;
 close all
 clc;

% Input participant's code name, condition, cadence
Participant = input('What is the participant''s code name ?', 's');
if ~exist(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant '/Audio'], 'dir')
    mkdir(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant '/Audio'])
    BPM = str2double(input('What is the participant''s preferred cadence (steps/min) ?', 's'));
    save(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant '/Audio/preferredBPM.mat'], 'BPM');
else
    % Load preferred cadence if already input
    load(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant '/Audio/preferredBPM.mat']);
end
Condition = input('What is the condition ?', 's');

if strcmp(Condition, 'RAC_preferredWalk')
elseif strcmp(Condition, 'RAC_slowWalk')
    BPM = BPM-(BPM*0.1);
elseif strcmp(Condition, 'RAC_fastWalk')
    BPM = BPM+(BPM*0.1);
else
    warning('Condition does not exist. Please run the code again.');
    (pause);
end
bpm2sec = (60/BPM);

% Load drum sound and convert .wav to vector
[audio, Freq] = audioread('C:\Users\lim\Desktop\Scripts_NeuroBiomec\Metronome\Sounds\projetRAC_C3.wav');
[audioDev, FreqDev] = audioread('C:\Users\lim\Desktop\Scripts_NeuroBiomec\Metronome\Sounds\projetRAC_B3.wav');

% Randomize when random tones will be heard
Time = 300*Freq;
onsetDeviant = randperm(BPM,5); %Randomize which beat will be deviant for each minute 
Beats = zeros(Time+1,1);
Beats(1:bpm2sec*Freq:Time) = 1; Beats(1) = [];
% Set one deviant tone per minute
iBeat = 0;
iDeviant = 1;
for iTime = 1:length(Beats)
    if Beats(iTime) == 1
        iBeat = iBeat+1;
        if iBeat == onsetDeviant(iDeviant)
            Beats(iTime) = 2;
            if iDeviant < 5
                iDeviant = iDeviant + 1;
                onsetDeviant(iDeviant) = onsetDeviant(iDeviant) + BPM*(iDeviant-1);
            end
        end
    end
end

% Create an audio band for which drum beats occur at cadence-matching frequency
audioBand = zeros(Time,2) ;
for iTime = 1:Time
    
    if Beats(iTime) == 1
        audioBand(iTime : iTime+length(audio)-1,:) = audio;        
    elseif Beats(iTime) == 2
        audioBand(iTime : iTime+length(audio)-1,:) = audioDev;
    end
    
end

% Save audio file of metronome
audiowrite('C:\Users\lim\Desktop\Scripts_NeuroBiomec\Metronome\titi.wav',audioBand,Freq)

% Save other audio-related variables
save(['/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/Behavioural/' Participant '/Audio/' Condition '.mat'], 'audioBand', 'onsetDeviant', 'BPM');

% Play metronome
% sound(audioBand, Freq);