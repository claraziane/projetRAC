%% Preprocessing : phase 1
% - Remove data before and after triggers
% - Remove DC offset
% - High pass at 0.1 Hz
% - Rereference according to the common average
% - Save preprocessed data in preprocessed folder

close all;
clear all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')    
    pathData = '\\10.89.24.15\j\Projet_RAC\DATA\';
    pathFinal = '\\10.89.24.15\q\Projet_RAC\DATA\Preprocessed\EEG';
    pathEEGLab = 'C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0';
else
    pathData = '/Volumes/10.89.24.15/Projet_RAC/DATA/';
%     pathFinal = '/Volumes/Clara/UdeM/S2M/projetRAC/DATA/Preprocessed/EEG/';
    pathEEGLab = '/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1';
end

addpath(pathEEGLab)
directoryRoot  = [pathData 'RAW/'];
extensionSet   = '.set';

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'};
% Conditions = {'RAC_preferredWalk'; 'RAC_preferredRest'; ...
%             'noRAC_preferredWalk'; 'noRAC_preferredRest'};
Conditions = {'RAC_slowWalk';   'RAC_fastWalk'; ...
            'noRAC_slowWalk'; 'noRAC_fastWalk'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
for iParticipant = 2:length(Participants)
    directoryRaw = fullfile(directoryRoot, Participants{iParticipant}, '/EEG/');

    for iCondition = 1%:length(Conditions)

        fileRead = [Conditions{iCondition} extensionSet];

        EEG = pop_loadset('filename', fileRead,'filepath', directoryRaw); % Loads an EEG dataset
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); % Edits/saves EEG dataset structure information
        disp(Conditions{iCondition})

        % Rereferencing according to common average
        EEG = pop_chanedit(EEG, 'append',63,'changefield',{64,'labels','FCz'},'changefield',{64,'theta','0'},'changefield',{64,'radius','0.128'},'changefield',{64,'X','0.391'},'changefield',{64,'Y','0'},'changefield',{64,'Z','0.921'},'changefield',{64,'sph_theta','0'},'changefield',{64,'sph_phi','67'},'changefield',{64,'sph_radius','85'},'setref',{'1:64','FCz'});
        EEG = pop_reref(EEG, [],'refloc',struct('labels',{'FCz'},'type',{''},'theta',{0},'radius',{0.128},'X',{0.391},'Y',{0},'Z',{0.921},'sph_theta',{0},'sph_phi',{67},'sph_radius',{85},'urchan',{[]},'ref',{'FCz'},'datachan',{0}));
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');

        % Remove baseline of the signal (must be before filtering)
        EEG = pop_rmbase(EEG, [],[]);

        % Highpass filter at 0.1 Hz
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.1);

        % Keep only signal inbetween start and end triggers
        triggers =  [EEG.event.latency];
        for iTrigger = 1:length(triggers)
            if triggers(iTrigger) ~= 1
                if triggers(iTrigger) == min(triggers(iTrigger:end)) && ~exist('triggerStart','var')
                    triggerStart =  triggers(iTrigger);
                end
                if triggers(iTrigger) < triggers(end) && triggers(iTrigger) > triggerStart * 2 && ~exist('triggerEnd','var')
                    triggerEnd = triggers(iTrigger);
                end
            end
        end

        % Check triggers (must be 5 minutes elapsed between triggerStart and triggerEnd)
        triggerTime = ((triggerEnd-triggerStart)/ALLEEG.srate)/60;
        while round(triggerTime) ~= 5
            triggerEnd = triggerStart + (ALLEEG.srate*60*5);
            triggerTime = ((triggerEnd-triggerStart)/ALLEEG.srate)/60;
        end
        EEG = pop_select(EEG,'time',[triggerStart/ALLEEG.srate triggerEnd/ALLEEG.srate]);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'overwrite', 'on', 'gui','off'); % Edit/save EEG dataset structure information
        
        % Re-reference according to the common average (again)
        EEG = pop_reref( EEG, []);
        EEG = eeg_checkset(EEG);
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off');        
    
        % Save new .set file in preprocessed folder
        EEG = pop_saveset(EEG, 'filename', fileRead, 'filepath', [pathFinal '/' Participants{iParticipant} '/']); %pop_saveset saves one or more EEG dataset structures
        ALLEEG = [];
        clear triggerStart triggerEnd;

    end

end


