%% Preprocessing : phase 2
% -Removes ICAs

close all;
clear all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
    pathData = '\\10.89.24.15\q\Projet_RAC\DATA\Preprocessed\EEG\'; %pathData = '\\10.89.24.15\e\Bureau\Clara\UdeM\S2M\projetRAC\DATA\Preprocessed\EEG\';
    pathEEGLab = 'C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0';
    load([pathData 'All\icaReject.mat'])
else
    pathData = '/Volumes/Projet_RAC/DATA/Preprocessed/EEG/'; %pathData = '/Volumes/10.89.24.15/Projet_RAC/DATA/';
    pathEEGLab = '/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1';
    load([pathData 'All/icaReject.mat'])
end

addpath(pathEEGLab)
extensionRoot  = '_ica.set';
extensionFinal = '_icaClean.set';

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P10'; 'P11'}; %'P06'; 'P08'; 'P09'; 'P12'
Conditions = {'RAC_slowWalk';  ...
            'noRAC_slowWalk'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 5:length(Participants)

    for iCondition = 1:length(Conditions)

        directory   = fullfile(pathData, Participants{iParticipant}, '/');
        if ~exist(directory, 'dir')
            mkdir(directory)
        end
        fileRead  = [Conditions{iCondition} extensionRoot];
        fileWrite = [Conditions{iCondition} extensionFinal];

        EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','on'); % Edits/saves EEG dataset structure information

        iICA = 1;
        while(1)
            icaRemove = input('Should some ICAs be removed from the data ?, Y/N [Y]:', 's');
            if icaRemove == 'Y'
                icaReject.([Participants{iParticipant}]).(Conditions{iCondition})(iICA) = input('Which ICA should be removed from the data ?');
                iICA = iICA + 1;
            elseif icaRemove == 'N'
                break
            end
        end

        % Remove identified bad components
        EEG = pop_subcomp(EEG, icaReject.([Participants{iParticipant}]).(Conditions{iCondition})(:), 0); %removes ICA from EEG and subtracts their activities from the data
        eeg_checkset(EEG);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

        % Rereference
        EEG = pop_reref(EEG, []);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

        % Save new .set file in preprocessed folder
        EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', directory); %pop_saveset saves one or more EEG dataset structures
        save([pathData 'All/icaReject'], 'icaReject');

        ALLEEG = [];

    end

end


