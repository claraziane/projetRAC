clear all;
% close all;
clc;

Computer = 'Windows'; %'Windows' or 'Macintosh'

if strcmp(Computer,'Macintosh')
    pathData = '/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/RAC/All/RAC';
    pathRaw = '/Volumes/10.89.24.15/Projet_RAC/DATA/RAW/';
    pathResults = '/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/DATA/RAC/All/RAC';
    pathFunctions = '/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Scripts/RAC/Functions';
    load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/RAC/All/RAC.mat');
elseif strcmp(Computer,'Windows')
    pathData = '\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\RAC\All\RAC';
    pathRaw = '\\10.89.24.15\j\Projet_RAC\DATA\RAW\';
    pathResults = 'C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\DATA\RAC\All\RAC';
    pathFunctions = 'C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\dataAnalysis\Scripts\RAC\Functions';
    load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\RAC\All\RAC.mat');
    addpath '\\10.89.24.15\j\Projet_RAC\dataAnalysis\Scripts\Clara\EEG\Functions'
end

addpath(pathFunctions)
extensionRoot  = '.c3d';

Participants = {'P08'};
Conditions   = {'RAC_preferredWalk'; 'RAC_slowWalk'; 'RAC_fastWalk'; 'RAC_preferredRest'};

for iParticipant = 1:length(Participants)
    load([pathRaw Participants{iParticipant} '/Audio/preferredBPM.mat'])

    for iCondition = 1:length(Conditions)

        % Define BPM
        if strcmp(Conditions{iCondition}, 'RAC_slowWalk')
            bpmInitial = round(preferredBPM * 0.9); %10% slower
        elseif strcmp(Conditions{iCondition}, 'RAC_fastWalk')
            bpmInitial = round(preferredBPM * 1.1); %10% faster
        elseif strcmp(Conditions{iCondition}(1:13), 'RAC_preferred')
            bpmInitial = preferredBPM;
        end

        % Extact audio data from structure
        Audio = RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).stimData(:,1);
        Freq  = RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).acquisitionFrequency(1,1);

        % Extract beat frequency, BPM, and IOI
        [beatFreq, BPM, IOI, beatOnset] = getBeat_FFT(Audio, Freq, bpmInitial);

        % Store data in structure
        RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).beatOnset(:,1)     = beatOnset; % Store beat onsets in structure      
        RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).beatFrequency(1,1) = beatFreq; % Store frequency in structure (other method)
        RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).BPM(1,1)           = BPM; % Store BPM in structure

        % Save structure
        save(pathData, 'RAC');
        save(pathResults, 'RAC');
        close all;

    end

end

