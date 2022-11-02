%% Preprocessing: phase 2
% - Run ICA (once for all conditions)
% - Save .set file with ICAs in preprocessed folder

close all;
clear all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
    pathData = '\\10.89.24.15\q\Projet_RAC\DATA\Preprocessed\EEG\';
    pathEEGLab = 'C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0';
else
%     pathData = '/Volumes/Clara/UdeM/S2M/projetRAC/DATA/Preprocessed/EEG/';
    pathEEGLab = '/Users/claraziane/Documents/AcadeÌ?mique/Informatique/MATLAB/eeglab2021.1';
end

addpath(pathEEGLab)
extensionRoot  = '_chanClean.set';
extensionFinal = '_ica.set';

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P10'; 'P11'};
Conditions = {'RAC_fastWalk';  ...
            'noRAC_fastWalk'};
% Conditions = {'RAC_preferredWalk';   'RAC_slowWalk';   'RAC_fastWalk';   'RAC_preferredRest'; ...
%             'noRAC_preferredWalk'; 'noRAC_slowWalk'; 'noRAC_fastWalk'; 'noRAC_preferredRest'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 1:length(Participants)
    directory   = fullfile(pathData, Participants{iParticipant}, '/');
    directoryResults = fullfile('Q:\Projet_RAC\DATA\Preprocessed\EEG\', Participants{iParticipant}, '/');

    for iCondition = 1:length(Conditions)
        fileRead = [Conditions{iCondition} extensionRoot];
        fileWrite = [Conditions{iCondition} extensionFinal];

        EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); % Edits/saves EEG dataset structure information

        % Run ICAs
        EEG = pop_runica(EEG, 'extended',1,'pca', size(EEG.data,1)-1,'interupt','on'); %EEG = pop_runica(EEG, 'icatype','runica','concatcond','on','options',{'extended',1});
        EEG = eeg_checkset(EEG);

        % Save
        EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', directoryResults); %pop_saveset saves one or more EEG dataset structures

        ALLEEG = [];

    end

end