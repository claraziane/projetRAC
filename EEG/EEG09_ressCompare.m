clear all;
close all;
clc;

[ret, Computer] = system('hostname');
if strcmp(Computer(1:5), 'BIMEC')
    pathResults = '\\10.89.24.15\e\Bureau\Clara\UdeM\S2M\projetRAC\dataAnalysis\Results\EEG\';
    addpath(genpath('C:\Users\p1208638\Documents\MATLAB\eeglab_current\eeglab2021.0'));
    addpath(genpath('\\10.89.24.15\e\Bureau\Clara\UdeM\S2M\projetSAB\Analyse_Donnees\drawingFigures\')); %To draw figures
    load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\Kinetics\All\stepData.mat')
    load('\\10.89.24.15\j\Projet_RAC\DATA\Preprocessed\RAC\All\RAC.mat')
else
    pathResults = '/Volumes/Seagate/Clara/UdeM/projetRAC/dataAnalysis/Results/EEG/';
    addpath(genpath('/Users/claraziane/Documents/Académique/Informatique/MATLAB/eeglab2021.1'));
    addpath '\\10.89.24.15\j\Projet_RAC\dataAnalysis\Scripts\Clara\EEG\Functions'
    addpath(genpath('/Volumes/Clara/UdeM/S2M/projetSAB/Analyse_Donnees/drawingFigures/')); %To draw figures
    load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/RAC/All/RAC.mat');
    load('/Volumes/10.89.24.15/Projet_RAC/DATA/Preprocessed/Kinetics/All/stepData.mat');
end

Conditions = {'RAC_preferredWalk',  'RAC_slowWalk',    'RAC_preferredRest'; ...
            'noRAC_preferredWalk', 'noRAC_slowWalk', 'noRAC_preferredRest'};
analysisType = 'flanking_noEpoch/';
load([pathResults 'All/' analysisType 'resultsEEG.mat'])
Color   =  [rgb('DarkOrange');  rgb('DarkOrange'); rgb('DodgerBlue')];
Colors  = [rgb('DarkOrange'); rgb('DarkOrange'); rgb('DodgerBlue')];
Titles  = {'Marche Préférentielle'; 'Marche Lente'; 'Repos'};
titleFigs  = {'preferredWalk'; 'slowWalk'; 'Rest'};
xLabels = {'SAR'; 'Controle'};
power                 = nan(11,3);
powerControl          = nan(11,3);
stabilityIndex        = nan(11,3);
stabilityIndexControl = nan(11,3);

eeglab;
for iCompare = 1%:length(Conditions)
    
    if iCompare == 1
        Participants = {       'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'};
    elseif iCompare == 2
        Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05';                      'P10'; 'P11'};
    else
        Participants = {'P01'; 'P02'; 'P03'; 'P04'; 'P05'; 'P06'; 'P08'; 'P09'; 'P10'; 'P11'; 'P12'};
    end
    Index = 1:2:length(Participants)*size(Conditions,2);

    for iParticipant = 1:length(Participants)
        
        %% Data extraction
        if strcmp(Conditions{1,iCompare}(end-3:end), 'Walk')
            Freq        = round(stepData.([Participants{iParticipant}]).([Conditions{1,iCompare}]).stepFreq,2);
            FreqControl = round(stepData.([Participants{iParticipant}]).([Conditions{2,iCompare}]).stepFreq,2);
        elseif strcmp(Conditions{1,iCompare}(end-3:end), 'Rest')
            Freq        = round(RAC.([Participants{iParticipant}]).([Conditions{1,iCompare}]).beatFrequency,2);
            FreqControl = round(RAC.([Participants{iParticipant}]).([Conditions{1,iCompare}]).beatFrequency,2);
        end

        EEGrate = 2500/3;
        ress            = resultsEEG.([Participants{iParticipant}]).([Conditions{1,iCompare}]).RESStime';
        ressControl     = resultsEEG.([Participants{iParticipant}]).([Conditions{2,iCompare}]).RESStime';
        ressSNR         = resultsEEG.([Participants{iParticipant}]).([Conditions{1,iCompare}]).ressSNR';
        controlSNR      = resultsEEG.([Participants{iParticipant}]).([Conditions{2,iCompare}]).ressSNR';      
        chanLocs        = resultsEEG.([Participants{iParticipant}]).([Conditions{1,iCompare}]).chanLocs;
        chanLocsControl = resultsEEG.([Participants{iParticipant}]).([Conditions{2,iCompare}]).chanLocs;
        if iCompare == size(Conditions,2)
            map             = (-1) * resultsEEG.([Participants{iParticipant}]).([Conditions{1,iCompare}]).Map;
            mapControl      = (-1) * resultsEEG.([Participants{iParticipant}]).([Conditions{2,iCompare}]).Map;
        else
            map             = resultsEEG.([Participants{iParticipant}]).([Conditions{1,iCompare}]).Map;
            mapControl      = resultsEEG.([Participants{iParticipant}]).([Conditions{2,iCompare}]).Map;
        end

        %% Topography
        figure(iCompare+1);
        subplot(2,length(Participants), iParticipant);...
            topoplot(map, chanLocs, 'maplimits', [-0.7 0.7], 'numcontour', 0, 'conv', 'off', 'electrodes', 'on', 'shading', 'interp'); hold on;
            title(Participants{iParticipant})
        subplot(2, length(Participants), iParticipant+length(Participants));...
            topoplot(mapControl, chanLocsControl, 'maplimits', [-0.7 0.7], 'numcontour', 0, 'conv', 'off', 'electrodes', 'on', 'shading', 'interp');

        %% Power

        nFFT = ceil(EEGrate/0.01);
        Hz    = linspace(0,EEGrate,nFFT);
        freqIndex        = dsearchn(Hz', Freq);
        freqIndexControl = dsearchn(Hz', FreqControl);
        
        % My method
%         ressx    = (abs(fft(ress, nFFT, 1)).^2)/length(ress);
%         controlx = (abs(fft(ressControl, nFFT, 1)).^2)/length(ressControl);
%         figure; plot(Hz,ressx); hold on; plot(Hz,controlx); xlim([0 25])
%         power(iParticipant)        = max(ressx(freqIndex-5:freqIndex+5));
%         powerControl(iParticipant) = max(controlx(freqIndexControl-5:freqIndexControl+5));
%        
%         % Matlab method
%         ressx = fft(ress,nFFT);
%         ressx = ressx(1:nFFT/2+1);
%         ressPower = (1/(EEGrate*nFFT)) * abs(ressx).^2;
%         ressPower(2:end-1) = 2*ressPower(2:end-1);
%         freq = 0:EEGrate/nFFT:EEGrate/2;
%         
%         controlx = fft(ressControl,nFFT);
%         controlx = controlx(1:nFFT/2+1);
%         controlPower = (1/(EEGrate*nFFT)) * abs(controlx).^2;
%         controlPower(2:end-1) = 2*controlPower(2:end-1);
% 
%         figure; plot(freq,ressPower); hold on; plot(freq,controlPower); xlim([0 25])
% 
%         % Etienne's method      
%         [ressFFT, ressF, Mean_Freq, Med_Freq, ressPSD]          = spec_fft(ress, EEGrate, 0);
%         [controlFFT, controlF, Mean_Freq, Med_Freq, controlPSD] = spec_fft(ressControl, EEGrate, 0);
%         freqIndex        = dsearchn(ressF', Freq);
%         freqIndexControl = dsearchn(controlF', FreqControl);
%         figure; plot(ressF,ressPSD); hold on; plot(controlF, controlPSD); xlim([0 25]);
% 
%         % Periodogram method
%         [ressPSD, ressF] = periodogram(ress, [], nFFT, EEGrate);
%         [controlPSD, controlF] = periodogram(ressControl, [], nFFT, EEGrate);
%         freqIndex        = dsearchn(ressF, Freq);
%         freqIndexControl = dsearchn(controlF, FreqControl);
%         figure; plot(ressF,ressPSD); hold on; plot(controlF, controlPSD); xlim([0 25]);
% 
%         power(iParticipant, iCompare)        = max(ressPSD(freqIndex-5:freqIndex+5));
%         powerControl(iParticipant, iCompare) = max(controlPSD(freqIndexControl-5:freqIndexControl+5));
           
        % RESS method
        power(iParticipant, iCompare)        = max(ressSNR(freqIndex-5:freqIndex+5));
        powerControl(iParticipant, iCompare) = max(controlSNR(freqIndexControl-5:freqIndexControl+5));
     
        %% Stability index
        stabilityIndex(iParticipant,iCompare)        = resultsEEG.([Participants{iParticipant}]).([Conditions{1,iCompare}]).stabilityIndex;
        stabilityIndexControl(iParticipant,iCompare) = resultsEEG.([Participants{iParticipant}]).([Conditions{2,iCompare}]).stabilityIndex;

        clear ress ressx ressControl controlx chanLocs map mapControl mapDiff
    end

    %% Stats
    powerDelta(:,iCompare) = power(:,iCompare) - powerControl(:,iCompare);
    siDelta(:,iCompare)    = stabilityIndex(:,iCompare) - stabilityIndexControl(:,iCompare);

    [powerH(iCompare),  powerP(iCompare),  powerCI(iCompare,:),   powerStats] = ttest(power(~isnan(power(:,iCompare)),iCompare), powerControl(~isnan(powerControl(:,iCompare)), iCompare));
    [powerDH(iCompare), powerDP(iCompare), powerDCI(iCompare,:), powerDStats] = ttest(powerDelta(~isnan(powerDelta(:,iCompare)), iCompare));
    [siH(iCompare),     siP(iCompare),     siCI(iCompare,:),         siStats] = ttest(stabilityIndex(~isnan(stabilityIndex(:,iCompare)), iCompare), stabilityIndexControl(~isnan(stabilityIndexControl(:,iCompare)), iCompare));
    [siDH(iCompare),    siDP(iCompare),    siDCI(iCompare,:),       siDStats] = ttest(siDelta(~isnan(siDelta(:,iCompare)), iCompare));
    
    %% Plotting

     % Power
     figure(5);
     subplot(1,3,iCompare); fig5 = UnivarScatter([power(:,iCompare), powerControl(:,iCompare)], 'PointSize', 200, 'LineWidth', 1.5, 'MarkerFaceColor', Color(iCompare,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
     plot([fig5(:,1) fig5(:,2)]', [power(:,iCompare), powerControl(:,iCompare)]', 'Color', Colors(iCompare,:), 'LineWidth', 1, 'LineStyle', '-'); hold on;
     ax1 = gca;
     set(ax1, 'xticklabel', xLabels);
     set(ax1, 'FontWeight', 'bold', 'FontSize', 16);
     if powerH(iCompare) == 1
         ylim = get(ax1, 'ylim');
         set(ax1, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
         ylim = get(ax1, 'ylim');
         plot(1.5, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
     end
     if iCompare == 1
         ylabel('SNR (Delta)')
     end
     title(Titles{iCompare})

     % Stability index
     figure(6);
     subplot(1,3,iCompare); fig6 = UnivarScatter([stabilityIndex(:,iCompare), stabilityIndexControl(:,iCompare)], 'PointSize', 200, 'LineWidth', 1.5, 'MarkerFaceColor', Color(iCompare,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines','Compression', 2); hold on;
     plot([fig6(:,1) fig6(:,2)]', [stabilityIndex(:,iCompare), stabilityIndexControl(:,iCompare)]', 'Color', Colors(iCompare,:), 'LineWidth', 1.5, 'LineStyle', '-'); hold on;
     ax2 = gca;
     set(ax2, 'xticklabel', xLabels);
     set(ax2, 'FontWeight', 'bold', 'FontSize', 16);
     if siH(iCompare) == 1
         ylim = get(ax2, 'ylim');
         set(ax2, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
         ylim = get(ax2, 'ylim');
         plot(1.5, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
     end
     if iCompare == 1
         ylabel('Index de Stabilité')
     end
     title(Titles{iCompare})

end

%% More plotting
saveas(figure(iCompare-1), [pathResults 'All/' analysisType '/topo_' titleFigs{1} '.png']);
saveas(figure(iCompare), [pathResults 'All/' analysisType '/topo_' titleFigs{2} '.png']);
saveas(figure(iCompare+1), [pathResults 'All/' analysisType '/topo_' titleFigs{3} '.png']);
saveas(figure(5), [pathResults 'All/' analysisType '/fig_Power.png']);
saveas(figure(6), [pathResults 'All/' analysisType '/fig_stabilityIndex.png']);

% Power
figure(iCompare+4);
UnivarScatter(powerDelta, 'PointStyle', '^', 'PointSize', 200, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax3 = gca;
xlim = get(ax3, 'xlim');
plot([xlim(1) xlim(2)], [0 0], 'color', [0.80,0.80,0.80]);
set(ax3, 'xticklabel', Titles);
ylabel('SNR (Delta)')
set(ax3, 'FontWeight', 'bold', 'FontSize', 16);
ylim = get(ax3, 'ylim');
set(ax3, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
ylim = get(ax3, 'ylim');
for iCompare = 1:3
    if powerDH(iCompare) == 1
        plot(iCompare, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
    end
end
title('Puissance')
saveas(figure(iCompare+4), [pathResults 'All/' analysisType '/fig_powerDelta.png']);

figure(iCompare+5);
UnivarScatter(siDelta, 'PointStyle', '^', 'PointSize', 200, 'LineWidth', 1.5, 'MarkerFaceColor', Colors(:,:),'MarkerEdgeColor', 'none', 'Whiskers', 'lines', 'Compression', 2); hold on;
ax4 = gca;
xlim = get(ax4, 'xlim');
plot([xlim(1) xlim(2)], [0 0], 'color', [0.80,0.80,0.80]);
set(ax4, 'xticklabel', Titles);
ylabel('Index de Stabilité (Delta)')
set(ax4, 'FontWeight', 'bold', 'FontSize', 16);
ylim = get(ax4, 'ylim');
set(ax4, 'ylim', [ylim(1) ylim(2)+((ylim(2)-ylim(1))*0.10)]);
ylim = get(ax4, 'ylim');
for iCompare = 1:3
    if siDH(iCompare) == 1
        plot(iCompare, ((ylim(2)-ylim(1))*0.93)+ylim(1), '*', 'color', 'black', 'LineWidth', 1.5,  'MarkerSize', 30);
    end
end
title('Index de Stabilité')
saveas(figure(iCompare+5), [pathResults 'All/' analysisType '/fig_stabilityIndexDelta.png']);
