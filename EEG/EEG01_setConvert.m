%% This script converts .eeg files to .set files
% -Load .eeg file
% -Save as .set file
% Note: if file name was changed after recording, must change EEG structure accordingly and resave file
close all;
clear all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
    directoryRoot = '\\10.89.24.15\j\Projet_RAC\DATA\';
    addpath('C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0');
    addpath('C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0\plugins\Fieldtrip-lite20220917')
    chanStr = 'C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp';
else
    directoryRoot = '/Volumes/10.89.24.15/Projet_RAC/DATA/';
    addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/fieldtrip-20211102')
    addpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1')
    chanStr = '/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp';
end

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'};
Conditions = {'RAC_preferredWalk';   'RAC_slowWalk';   'RAC_fastWalk';   'RAC_preferredRest'; ...
            'noRAC_preferredWalk'; 'noRAC_slowWalk'; 'noRAC_fastWalk'; 'noRAC_preferredRest'};

extensionRaw  = sprintf('.eeg');
extensionSet  = sprintf('.set');
 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for iParticipant = 1:length(Participants)

    for iCondition = 1:length(Conditions)

        directoryRaw = fullfile(directoryRoot,'RAW/', Participants{iParticipant}, '/EEG/');

        fileRead  = [directoryRaw Conditions{iCondition} extensionRaw]; %[ directoryRaw Conditions{iCondition} extensionRaw]
        fileWrite = [directoryRaw Conditions{iCondition} extensionSet];

        EEG = pop_fileio(fileRead); % Import .eeg data into EEGLab
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); % Edit/save EEG dataset structure information
        %EEG = pop_chanedit(EEG, 'lookup','C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
        EEG = pop_chanedit(EEG, 'lookup',chanStr);      
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); % Stores specified EEG dataset(s) in the ALLEEG variable

        EEG = pop_saveset( EEG, fileWrite); % Saves one or more EEG dataset structures
        ALLEEG = pop_delset(ALLEEG, 1); %Deletes data from ALLEEG

    end

end
