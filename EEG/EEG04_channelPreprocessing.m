%% Preprocessing : phase 2
% -Removes bad channels

close all;
clear all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
    pathData   = '\\10.89.24.15\q\Projet_RAC\DATA\Preprocessed\EEG\';
    pathEEGLab = 'C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0';
else
    pathData   = '/Volumes/Projet_RAC/DATA/Preprocessed/EEG/';
    pathEEGLab = '/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1';
end

addpath(pathEEGLab)
extensionRoot  = '_events.set';
extensionFinal = '_chanClean.set';
load([pathData '/All/chanReject.mat'])

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P10'; 'P11'};
Conditions = {'RAC_fastWalk';  ...
            'noRAC_fastWalk'};
% Conditions = {'RAC_preferredWalk';   'RAC_slowWalk';   'RAC_fastWalk';   'RAC_preferredRest'; ...
%             'noRAC_preferredWalk'; 'noRAC_slowWalk'; 'noRAC_fastWalk'; 'noRAC_preferredRest'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 4:length(Participants)

    for iCondition = 1:length(Conditions)

        directory = fullfile(pathData, Participants{iParticipant}, '/');
        fileRead  = [Conditions{iCondition} extensionRoot];
        fileWrite = [Conditions{iCondition} extensionFinal];

        EEG = pop_loadset('filename', fileRead,'filepath', directory); % Loads an EEG dataset
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); % Edits/saves EEG dataset structure information

        % Identify bad channels
        iChan = 1;
        chanReject.([Participants{iParticipant}]).(Conditions{iCondition})(iChan) = NaN;
        
        while(1)
            chanRemove = input('Should some channels be removed from the data ?, Y/N [Y]:', 's');
            if chanRemove == 'Y'
                chanReject.([Participants{iParticipant}]).(Conditions{iCondition})(iChan) = input('Which channel should be removed from the data ?');
                iChan = iChan + 1;
            elseif chanRemove == 'N'
                break
            end
        end

        % Remove identified bad channels
        if ~isnan(chanReject.([Participants{iParticipant}]).(Conditions{iCondition})(1))
            chanReject.([Participants{iParticipant}]).(Conditions{iCondition}) = sort(chanReject.([Participants{iParticipant}]).(Conditions{iCondition}));
            EEG = pop_select( EEG, 'nochannel', chanReject.([Participants{iParticipant}]).(Conditions{iCondition})(1:iChan-1));
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

            % Rereference
            EEG = pop_reref( EEG, []);
            [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');
        end

        % Save new .set file in preprocessed folder
        EEG = pop_saveset(EEG, 'filename', fileWrite, 'filepath', directory); %pop_saveset saves one or more EEG dataset structures
        ALLEEG = [];

        save([pathData '/All/chanReject'], 'chanReject');
    end

end


