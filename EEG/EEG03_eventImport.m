%% Importing events
% -Import all right step onsets within EEG structure
% -Import all left step onsets within EEG structure
% -Import all beat onsets within EEG structure

close all;
clear all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
    pathData   = '\\10.89.24.15\q\Projet_RAC\DATA\Preprocessed\EEG\';
    pathEEGLab = 'C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0';
    load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\RAC\All\RAC.mat');
    load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\Kinetics\All\stepData.mat');
else
%     pathData   = '/Volumes/Clara/UdeM/S2M/projetRAC/DATA/Preprocessed/EEG/';
    pathEEGLab = '/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1';
    load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/RAC/All/RAC.mat');
    load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/Kinetics/All/stepData.mat');
end

addpath(pathEEGLab)
extensionRoot  = '.set';
extensionFinal = '_events.set';

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'};
Conditions = {'RAC_slowWalk';   'RAC_fastWalk'; ...
            'noRAC_slowWalk'; 'noRAC_fastWalk'};
% Conditions = {'RAC_preferredWalk';   'RAC_slowWalk';   'RAC_fastWalk';   'RAC_preferredRest'; ...
%             'noRAC_preferredWalk'; 'noRAC_slowWalk'; 'noRAC_fastWalk'; 'noRAC_preferredRest'}; 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 1:length(Participants)
    directory = fullfile(pathData, Participants{iParticipant}, '/');

    for iCondition = 1:length(Conditions)

        fileRead  = [Conditions{iCondition} extensionRoot];
        fileWrite = [Conditions{iCondition} extensionFinal];

        EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','on'); % Edits/saves EEG dataset structure information

        %% Foot strikes

        if strcmpi(Conditions{iCondition}(end-3:end), 'Walk')

            % Extract events' acquisition frequency
            stepRate = stepData.([Participants{iParticipant}]).(Conditions{iCondition}).sRate;

            % Extract events from structure
            stepOnsetR   = stepData.([Participants{iParticipant}]).(Conditions{iCondition}).stepOnsetR;
            stepOnsetL   = stepData.([Participants{iParticipant}]).(Conditions{iCondition}).stepOnsetL;

            % Interpolate values to fit EEG acquisition frequency
            stepOnsetR   = round(stepOnsetR * (EEG.srate/stepRate)); %interp1(1:length(stepOnsetR), stepOnsetR, linspace(1,length(stepOnsetR), EEG.srate)) ;
            stepOnsetL   = round(stepOnsetL * (EEG.srate/stepRate));

            nEvents = length(EEG.event);
            for iEvent=1:length(stepOnsetR)
                EEG.event(nEvents+iEvent).type = 'stepR' ;
                EEG.event(nEvents+iEvent).latency = stepOnsetR(iEvent) ;
                EEG.event(nEvents+iEvent).duration = 1 ;
                EEG.event(nEvents+iEvent).urevent = nEvents+iEvent  ;
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

            nEvents = length(EEG.event);
            for iEvent=1:length(stepOnsetL)
                EEG.event(nEvents+iEvent).type = 'stepL' ;
                EEG.event(nEvents+iEvent).latency = stepOnsetL(iEvent) ;
                EEG.event(nEvents+iEvent).duration = 1 ;
                EEG.event(nEvents+iEvent).urevent = nEvents+iEvent  ;
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

        end

        %% Beat onsets

        if strcmpi(Conditions{iCondition}(1:1+2), 'RAC')

            % Extract events' acquisition frequency
            beatRate = RAC.([Participants{iParticipant}]).(Conditions{iCondition}).acquisitionFrequency;

            % Extract events from structure
            beatOnset = RAC.([Participants{iParticipant}]).(Conditions{iCondition}).beatOnset;

            % Interpolate values to fit EEG acquisition frequency
            beatOnset = round(beatOnset * (EEG.srate/beatRate));

            nEvents = length(EEG.event);
            for iEvent=1:length(beatOnset)
                EEG.event(nEvents+iEvent).type = 'RAC' ;
                EEG.event(nEvents+iEvent).latency = beatOnset(iEvent) ;
                EEG.event(nEvents+iEvent).duration = 1 ;
                EEG.event(nEvents+iEvent).urevent = nEvents+iEvent  ;
            end
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

        end

        % Save new _event.set file in preprocessed folder
        EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', directory); %pop_saveset saves one or more EEG dataset structures
        ALLEEG = [];
        clear stepOnsetR stepOnsetL beatOnset;

    end
end