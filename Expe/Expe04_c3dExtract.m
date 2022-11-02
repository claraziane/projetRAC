clear all;
close all;
clc;

% Declare paths
pathData = '\\10.89.24.15\j\Projet_RAC\DATA\';
pathEzc3d = 'C:\Users\p1208638\Documents\MATLAB\ezc3d_matlab';
load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\RAC\All\RAC.mat');
load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\Kinetics\All\stepData.mat')

addpath(pathEzc3d)
extensionRoot  = '.c3d';

Participants = {''}; %'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09';'P10'; 'P11'; 'P12'
Conditions   = {'RAC_preferredWalk'; 'RAC_slowWalk'; 'RAC_fastWalk'; 'RAC_preferredRest'; 'noRAC_preferredWalk'; 'noRAC_slowWalk'; 'noRAC_fastWalk'; 'noRAC_preferredRest'};

for iParticipant = 1:length(Participants)

    for iCondition = 1:length(Conditions)

        c3d = ezc3dRead([pathData 'RAW/' Participants{iParticipant} '/Nexus/' Conditions{iCondition} extensionRoot]);

        % Store analog acquisition frequency in structure
        RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).acquisitionFrequency(1,1) = c3d.parameters.ANALOG.RATE.DATA;

        %% Trigger data

        % Find trigger data in c3d file
        triggerFind = strfind(c3d.parameters.ANALOG.LABELS.DATA, 'Synchronization');
        for iIndex = 1:length(triggerFind)
            if cell2mat(triggerFind(iIndex)) == 1
                triggerIndex = iIndex;
            end
        end
        triggerData = c3d.data.analogs(:,triggerIndex);

        [triggerPks,triggerLocs] = findpeaks(triggerData, 'MinPeakHeight', 0, 'MinPeakDistance', c3d.parameters.ANALOG.RATE.DATA*60);
        figure; plot(triggerData); hold on
        plot(triggerLocs, triggerPks, 'r*')

        % Check triggers (must be 5 minutes elapsed between triggerStart and triggerEnd)
        triggerTime = ((triggerLocs(end)-triggerLocs(1))/c3d.parameters.ANALOG.RATE.DATA)/60;
        while round(triggerTime) ~= 5
            triggerLocs(end) = triggerLocs(1) + (c3d.parameters.ANALOG.RATE.DATA*60*5);
            triggerTime = ((triggerLocs(end)-triggerLocs(1))/c3d.parameters.ANALOG.RATE.DATA)/60;
            if triggerLocs(end) >= length(triggerData)
                triggerLocs(end) = length(triggerData);
                triggerLocs(1) = triggerLocs(end) - (5*60*c3d.parameters.ANALOG.RATE.DATA);
            end
        end

        %% Audio data
        
        % Only include RAC conditions
        if strcmp(Conditions{iCondition}(1:3), 'RAC')

            % Find metronome data in c3d file
            metronomeFind = strfind(c3d.parameters.ANALOG.LABELS.DATA, 'Metronome.1');
            for iIndex = 1:length(metronomeFind)
                if cell2mat(metronomeFind(iIndex)) == 1
                    metronomeIndex = iIndex;
                end
            end
            metronomeData = c3d.data.analogs(:,metronomeIndex);

            % Cut metronome trial with triger data
            metronomeData = metronomeData(triggerLocs(1):triggerLocs(end));

            % Store metronome data in structure
            RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).stimData(:,1) = metronomeData;

            % Save extracted data
            save([pathData 'Preprocessed\RAC\All\RAC'], 'RAC');
            save(['C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\DATA\RAC\All\RAC'], 'RAC');

        end

        %% Kinetic data
        stepData.([Participants{iParticipant}]).(Conditions{iCondition}).sRate = c3d.parameters.ANALOG.RATE.DATA;

        % Right foot
        % Find Fz data in c3d file
        rightFind = strfind(c3d.parameters.ANALOG.LABELS.DATA, 'Right_foot.Fz_R');
        for iIndex = 1:length(rightFind)
            if cell2mat(rightFind(iIndex)) == 1
                rightIndex = iIndex;
            end
        end
        kineticDataR = c3d.data.analogs(:,rightIndex);

        % Remove data outside of triggers
        kineticDataR = kineticDataR(triggerLocs(1):triggerLocs(end));
        stepData.([Participants{iParticipant}]).(Conditions{iCondition}).kineticDataR = kineticDataR;

        % Left foot
        % Find Fz data in c3d file
        leftFind = strfind(c3d.parameters.ANALOG.LABELS.DATA, 'Left_foot.Fz_L');
        for iIndex = 1:length(leftFind)
            if cell2mat(leftFind(iIndex)) == 1
                leftIndex = iIndex;
            end
        end
        kineticDataL = c3d.data.analogs(:,leftIndex);
        
        % Remove data outside of triggers
        kineticDataL = kineticDataL(triggerLocs(1):triggerLocs(end));
        stepData.([Participants{iParticipant}]).(Conditions{iCondition}).kineticDataL = kineticDataL;

        save([pathData 'Preprocessed\Kinetics\All\stepData'], 'stepData');
        save(['C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\DATA\Kinetics\All\stepData'], 'stepData');
    end

end

