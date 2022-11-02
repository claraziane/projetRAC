clear all;
close all;
clc;

Computer = 'Macintosh'; %'Windows' or 'Macintosh'

if strcmp(Computer,'Macintosh')
    load('/Volumes/10.89.24.15/Projet_RAC/dataAnalysis/Results/Beahavioural/All/resultsBehavioural.mat');
    load('/Volumes/Seagate/Clara/UdeM/projetRAC/dataAnalysis/Results/EEG/All/flanking_noEpoch_icaClean/resultsEEG.mat');
    load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/Kinetics/All/stepData.mat');
elseif strcmp(Computer,'Windows')
%     load('C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\dataAnalysis\Results\Behavioural\All\resultsBehavioural.mat');
%     load('C:\Users\p1208638\OneDrive - Universite de Montreal\S2M\projetRAC\dataAnalysis\Results\EEG\All\resultsEEG.mat');
end

Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'};
Conditions   = {'RAC_preferredWalk'};

for iCondition = 1

    for iParticipant = 1:length(Participants)

        stepFreq = stepData.([Participants{iParticipant}]).(Conditions{1,1}).sRate;


        stabilityIndex(iParticipant,1)   = resultsEEG.([Participants{iParticipant}]).([Conditions{iCondition}]).stabilityIndex;
        IBIdeviation(iParticipant,1)     = resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).IBIDeviation;
        IBIdeviation(iParticipant,1)     = IBIdeviation(iParticipant,1) / (stepFreq/1000); %Convert to milliseconds
        meanAsynchrony(iParticipant, 1)  = mean(resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).Asynchrony);
        meanAsynchrony(iParticipant, 1)  = meanAsynchrony(iParticipant, 1) / (stepFreq/1000); %Convert to milliseconds
        %         meanAsynchrony(iParticipant,1)   = mean(resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).circularAsynchrony); % En faisant la moyenne classique
        %         asynchronyMean(iParticipant,1)   = (resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).asynchronyMean) * 180/pi;           % En utilisant la toolbox pour faire la moyenne
        %         circularVariance(iParticipant,1) = resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).circularVariance*180/pi;
        %         meanPhaseAngle(iParticipant,1)   = mean(resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).phaseAngle); % En faisant la moyenne classique
        phaseAngleMean(iParticipant,1)   = resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).phaseAngleMean;   % En utilisant la toolbox pour faire la moyenne
        resultantLength(iParticipant,1)  = resultsBehavioural.([Participants{iParticipant}]).([Conditions{iCondition}]).resultantLength;
    end
    [rhoIBI,pIBI]     = corr(stabilityIndex, IBIdeviation,     'Type', 'Spearman');
    [rhoAsync,pAsync] = corr(stabilityIndex, meanAsynchrony,   'Type', 'Spearman');
    %     [rhoCV,pCV]       = corr(stabilityIndex, circularVariance, 'Type', 'Spearman');
    [rhoPA,pPA]       = corr(stabilityIndex, phaseAngleMean,   'Type', 'Spearman');
    [rhoRVL,pRVl]     = corr(stabilityIndex, resultantLength,  'Type', 'Spearman');

    h = figure(1);
    ax = gca;
    ax.FontSize = 14;
    subplot(2,2,1); scatter(stabilityIndex, meanAsynchrony);
    txt = (['\rho = ' num2str(rhoAsync)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FitBoxToText','on', 'FontSize', 14);
    xlabel('Stability Index (Hz)', 'FontSize', 14); ylabel('Mean Asynchrony (ms)', 'FontSize', 14);
    title('Mean Asynchrony', 'FontSize', 14);
    %     subplot(2,2,2); scatter(stabilityIndex, circularVariance);
    %         txt = (['\rho = ' num2str(rhoCV)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FitBoxToText','on', 'FontSize', 14);
    %         xlabel('Stability Index (Hz)'); ylabel('Circular Variance (°)');
    %         title('Circular Variance');
    %     subplot(3,2,3); scatter(stabilityIndex, meanPhaseAngle);
    % %         txt = (['\rho = ' num2str(rhoPA)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FitBoxToText','on', 'FontSize', 14);
    %         xlabel('Stability Index (Hz)'); ylabel('Relative Phase Angle (°)');
    %         title('Relative Phase Angle -Regular Mean');
    subplot(2,2,2); scatter(stabilityIndex, IBIdeviation);
    txt = (['\rho = ' num2str(rhoIBI)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FitBoxToText','on', 'FontSize', 14);
    xlabel('Stability Index (Hz)', 'FontSize', 14); ylabel('IBI Deviation', 'FontSize', 14);
    title('Interbeat Interval Deviation', 'FontSize', 14);
    subplot(2,2,3); scatter(stabilityIndex, phaseAngleMean);
    txt = (['\rho = ' num2str(rhoPA)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FitBoxToText','on', 'FontSize', 14);
    xlabel('Stability Index (Hz)', 'FontSize', 14); ylabel('Relative Phase Angle (°)', 'FontSize', 14);
    title('Relative Phase Angle', 'FontSize', 14);
    subplot(2,2,4); scatter(stabilityIndex, resultantLength);
    txt = (['\rho = ' num2str(rhoRVL)]); dim = [.2 .5 .3 .3]; annotation('textbox',dim,'String',txt, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FitBoxToText','on', 'FontSize', 14);
    xlabel('Stability Index (Hz)', 'FontSize', 14); ylabel('Resultant Vector Length', 'FontSize', 14);
    title('Resultant Vector Length', 'FontSize', 14);
    sgtitle('Behavioural Measures As A Function Of The Stability Index', 'FontSize', 16)

end

% plotregression(stabilityIndex', asynchronyMean',   'Mean Circular Asynchrony -Toolbox Mean', ...
%     stabilityIndex', circularVariance', 'Circular Variance', ...
%     stabilityIndex', meanPhaseAngle',   'Relative Phase Angle -Regular Mean', ...
%     stabilityIndex', phaseAngleMean',   'Relative Phase Angle -Toolbox Mean', ...
%     stabilityIndex', resultantLength',  'Resultant Vector Length',  ...
%     stabilityIndex', IBIdeviation',     'IBI Deviation (ratio)');