clear all;
close all;
clc;

Computer = 'Macintosh'; %'Windows' or 'Macintosh'

if strcmp(Computer,'Macintosh')
    pathData  = '/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/';
    pathResults = '/Volumes/10.89.24.15/Projet_RAC/dataAnalysis/Results/Beahavioural/All/';
    addpath('/Users/claraziane/OneDrive - Universite de Montreal/S2M/projetRAC/dataAnalysis/Toolbox/CircStat2012a/');
    load([pathData 'RAC/All/RAC.mat']);
    load([pathData 'Kinetics/All/stepData.mat']);
    load([pathResults 'resultsBehavioural.mat']);
elseif strcmp(Computer,'Windows')
    pathData = '\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\';
    pathResults = 'W:\Projet_RAC\dataAnalysis\Results\Beahavioural\All\';
    addpath('C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\dataAnalysis\Toolbox\CircStat2012a');
    load([pathData 'RAC\All\RAC.mat']);
    load([pathData 'Kinetics\All\stepData.mat']);
    load([pathResults 'resultsBehavioural.mat']);
end

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'}; %'P05';
Conditions   = {'RAC_preferredWalk'};% 'RAC_slowWalk'; 'RAC_fastWalk'};

for iParticipant = 9:10%1:length(Participants)

    for iCondition = 1:length(Conditions)

        %% Estimating period-matching accuracy (i.e., extent to which step tempo matches stimulus tempo) using IBI deviation

        % Extracting beat onsets
        beatOnset = RAC.([Participants{iParticipant}]).([Conditions{iCondition}]).beatOnset;

        % Extracting step onsets
        stepOnset = stepData.([Participants{iParticipant}]).([Conditions{iCondition}]).stepOnsetAll;

        % Matching step onsets to closest beat
        beatMatched = zeros(length(stepOnset),1);
        for iStep = 1:length(stepOnset)
            [minValue matchIndex] = min(abs(beatOnset-stepOnset(iStep)));
            beatMatched(iStep,1) = beatOnset(matchIndex);
            if iStep >= 2 && beatMatched(iStep,1) == beatMatched(iStep-1,1) 
                beatMatched(iStep) = beatOnset(matchIndex+1);
            end
        end

        % Make sure beatMatched(1) and beatMatched(2) are not the same (ex. if first beat was removed but not first step)
        if beatMatched(1) == beatMatched(2)
            beatMatched(1) = [];
            if length(stepOnset) > length(beatOnset)
                stepOnset(1) = [];
            elseif length(stepOnset) < length(beatOnset)
                beatOnset(1) = [];
            end
        elseif beatMatched(end) == beatMatched(end-1)
            beatMatched(end) = [];
            if length(stepOnset) > length(beatOnset)
                stepOnset(end) = [];
            elseif length(stepOnset) < length(beatOnset)
                beatOnset(end) = [];
            end
        end

        % Calculating interstep interval
        stepInterval = diff(stepOnset);

        % Calculating interbeat interval
        racInterval = diff(beatMatched);
        if ~isempty(racInterval(racInterval == 0))
            for iInt = 1:length(racInterval)
                if racInterval(iInt) == 0
                    if racInterval(iInt-1) > racInterval(iInt+1)
                        [M, minIndex] = min(abs(beatOnset-beatMatched(iInt)));
                        beatMatched(iInt) = beatOnset(minIndex-1);
                    elseif racInterval(iInt-1) < racInterval(iInt+1)
                        [M, minIndex] = min(abs(beatOnset-beatMatched(iInt)));
                        beatMatched(iInt) = beatOnset(minIndex+1);
                    end
                    racInterval = diff(beatMatched);
                end
            end
        end
       

        % Calculating IBI deviation
        IBI = mean(abs(stepInterval - racInterval))/mean(racInterval);

        % Storing IBI in results structure
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).IBIDeviation = IBI;


        %% Estimating phase-matching accuracy (i.e., the difference between step onset times and beat onset times) using circular asynchronies

        asynchrony           = stepOnset - beatMatched; % Used formula diff from Leow et al. 2017
        asynchronyNormalized = asynchrony(1:end-1)./stepInterval;
        asynchronyCircular   = asynchronyNormalized * 360;
        asynchronyRad        = asynchronyCircular * pi/180;
        asynchronyMean       = circ_mean(asynchronyRad, [], 1);
        figure; scatter(1,asynchronyCircular)

        % Running Rao's test (a not-significant test means participant failed to synchronize)
        [p U UC] = circ_raotest(asynchronyCircular);

        % Calculating circular variance
        [varianceCircular varianceAngular] = circ_var(asynchronyCircular/(180/pi)); % Degrees are converted to radians

        % Calculating phase angles (error measure of synchronization based on the phase difference between two oscillators)
        phaseAngle = 360*((stepOnset(1:end-1) - beatMatched(1:end-1))./racInterval);
%         if ~isempty(phaseAngle(phaseAngle == -inf))
%             phaseAngle(phaseAngle == -inf) = [];
%         end
        phaseRad = deg2rad(phaseAngle);
        phaseAngleMean = circ_mean(phaseRad, [], 1);

        % Calculating resultant vector length (expresses the stability of the relative phase angles over time)
        resultantLength = circ_r(phaseRad, [], [], 1);

        % Storing results in structure
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).Asynchrony = asynchrony;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).circularAsynchrony = asynchronyCircular;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).asynchronyMean = asynchronyMean;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).circularVariance = varianceCircular;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).pRao = p;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).phaseAngle = phaseAngle;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).phaseAngleMean = phaseAngleMean;
        resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).resultantLength = resultantLength;
        clear asynchrony asynchronyNormalized asynchronyCircular asynchronyMean varianceCircular p phaseAngle phaseAngleMean resultantLength beatOnset

    end

end

save([pathResults 'resultsBehavioural.mat'], 'resultsBehavioural');