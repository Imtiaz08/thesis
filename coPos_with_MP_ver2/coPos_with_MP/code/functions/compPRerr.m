function [prErr, SNR] = compPRerr(sampFactor, corr_spac, noChips_del, ...
    noChips_SNR, noise_phase, noise_amp, cb, sys, freqBands, comps, satID,...
    LOS_received, no_refls, delta_d, delta_P, const)

% sampFactor = 0.1 %24 % sampFactor/12 must be intenger for E5!
% corr_spac = 0.5 % spacing in EPL correlator, number of chips
% noChips_del = 2 % Number of chips in each direction to find correlation peak
% noChips_SNR = 1000 % Number of chips in each direction to decide noise level for SNR
% noise_phase = 0 % phase noise of receiver
% noise_amp = 0 % amplitude noise of receiver
% cb = 1e3; % clock bias of receiver (not yet used)
% 
% % General signal setting (extend to several frequencies)
% sys = 'GPS'
% freqBands = {'L1', 'L2'}
% comps = {'CA', 'CM'}
% satID = 30
% LOS_received = 0
% refl_present = 1
% delta_d = [26 26]
% delta_P = [0.7 0.7];


L1_corrmax = 157542;
L2_corrmax = 2455200;

for f = 1:1:length(freqBands)
    freqBand = freqBands{f};
    comp = comps{f};
    
    
    % Generate LOS time signal
    t_end = const.(sys).(freqBand).(comp).T_code;
    f_samp = sampFactor*const.(sys).(freqBand).(comp).f_carrier;
    t_vec = (0:1/f_samp:t_end-1/f_samp);
    code = codeGenerator(sys, freqBand, comp, satID, f_samp, const);
      [s_LOS, codeToCorr, ~]...
            = signalGenerator(const, sys, freqBand, comp, code,...
            t_vec, f_samp, noise_phase, noise_amp, 0);
    if LOS_received == 0
        s_LOS = zeros(1,length(t_vec));
    end
    
    % [corr, del] = xcorr(s_LOS.(sys).(freqBand).*exp(-2*1j*pi*const.(sys).(freqBand).(comp).f_carrier*t), codeToCorr.(sys).(freqBand), 1e5);
    % figure;
    % plot(del, corr)
    % hold on
    % plot(del, imag(corr), 'r')
    % plot(del, abs(corr), 'g')
    
    % Generate multipath time signal(s)
    s_MP_scaled = zeros(1,length(t_vec));
    for r = 1:1:no_refls
        time_delay_MP = delta_d(r,f)/const.c;
        [s_MP, ~, ~]...
            = signalGenerator(const, sys, freqBand, comp, code, t_vec, f_samp,...
            noise_phase, noise_amp, time_delay_MP);
        s_MP_scaled = s_MP_scaled + s_MP .* delta_P(r,f);
    end
    
    % Combine LOS and MP to complete time signal:
    
    % sqrt(mean(s_LOS.^2))
    % sqrt(mean(s_MP_scaled.^2))
    
    s_total = s_LOS + s_MP_scaled;
    
    % Number of samples in correlation form in each direction
    noSamplesInCorr_SNR = noChips_SNR * f_samp/const.GPS.L1.CA.f_chip;
    noSamplesInCorr_del = noChips_del * f_samp/const.GPS.L1.CA.f_chip;
    
    
    % Generate correlation form (large version, for SNR computation)
    [corrForm_SNR, delays_SNR] = xcorr(...
        s_total.*exp(-2*1j*pi*const.(sys).(freqBand).(comp).f_carrier*t_vec),...
        codeToCorr, noSamplesInCorr_SNR);
    
    chipVec_SNR = delays_SNR/f_samp;
    
    corrForm_SNR_abs = abs(corrForm_SNR);
    
    % Filter out compact correlation form, for delay tracking
    centerSample_SNR = floor(length(chipVec_SNR)/2)+1;
    chipVec_del = chipVec_SNR(centerSample_SNR - noSamplesInCorr_del:centerSample_SNR + noSamplesInCorr_del);
    corrForm_del_abs = corrForm_SNR_abs(centerSample_SNR - noSamplesInCorr_del:centerSample_SNR + noSamplesInCorr_del);
    
    % find reference early phase late samples
    timeInCorrSpac = corr_spac/const.(sys).(freqBand).(comp).f_chip;
    [~, esamp_ref] = min(abs(chipVec_del + timeInCorrSpac));
    [~, psamp_ref] = min(abs(chipVec_del));
    [~, lsamp_ref] = min(abs(chipVec_del - timeInCorrSpac));
    
    
    % Method to estimate peak:
    noSamplesInCorrSpac = timeInCorrSpac * f_samp;
    
    
    middleSamp = ceil(length(chipVec_del)/2);
    samplesForTracking = [noSamplesInCorrSpac + 1: 1 : length(chipVec_del) - noSamplesInCorrSpac];
    
    earlyLateDiff = NaN*ones(size(chipVec_del));
    phaseCorrPow = NaN*ones(size(chipVec_del));
    earlySamples = NaN*ones(size(chipVec_del));
    phaseSamples = NaN*ones(size(chipVec_del));
    lateSamples = NaN*ones(size(chipVec_del));
    
    for s = samplesForTracking
        
        earlySamples(s) = s - noSamplesInCorrSpac;
        phaseSamples(s) = s ;
        lateSamples(s) = s + noSamplesInCorrSpac;
        
        % Difference between early and late:
        earlyLateDiff(s) = abs(corrForm_del_abs(earlySamples(s)) - corrForm_del_abs(lateSamples(s)));
        
        % Correlation power of phase:
        phaseCorrPow(s) = corrForm_del_abs(phaseSamples(s));
    end
    trackingOutput = phaseCorrPow./earlyLateDiff;
    [~, bestSamp] = max(trackingOutput);
    
    
    % Decide SNR:
    nbefore = corrForm_SNR_abs(1:  floor(length(corrForm_SNR_abs)/2)+1-corr_spac/const.(sys).(freqBand).(comp).f_chip*f_samp);
    nafter = corrForm_SNR_abs( floor(length(corrForm_SNR_abs)/2)+1+corr_spac/const.(sys).(freqBand).(comp).f_chip*f_samp:end);
    
    
    % Save results to struct
    trackRes.(freqBand).maxPow = corrForm_del_abs(bestSamp);
    trackRes.(freqBand).bestSamp = bestSamp; % Decide delay error
    trackRes.(freqBand).MP_error = chipVec_del(bestSamp)*const.c; % Decide delay error
    trackRes.(freqBand).SNR = max(corrForm_SNR_abs)/ mean([nbefore nafter]);
    trackRes.(freqBand).trackingOutput = trackingOutput;
    trackRes.(freqBand).earlySamples = earlySamples;
    trackRes.(freqBand).phaseSamples = phaseSamples;
    trackRes.(freqBand).lateSamples = lateSamples;
    trackRes.(freqBand).earlyLateDiff = earlyLateDiff;
    trackRes.(freqBand).phaseCorrPow = phaseCorrPow;
    trackRes.(freqBand).chipVec_del = chipVec_del;
    trackRes.(freqBand).corrForm_del_abs = corrForm_del_abs;
    trackRes.(freqBand).esamp_ref = esamp_ref;
    trackRes.(freqBand).psamp_ref = psamp_ref;
    trackRes.(freqBand).lsamp_ref = lsamp_ref;
end
prErr = [trackRes.(freqBands{1}).MP_error, trackRes.(freqBands{2}).MP_error];
SNR =   [trackRes.(freqBands{1}).maxPow/L1_corrmax,  trackRes.(freqBands{2}).maxPow/L2_corrmax];
%SNR2  =[trackRes.(freqBands{1}).SNR  trackRes.(freqBands{2}).SNR ]

  
 
% %Visualization step
% close all
% figure;
% 
% for f = 1:1:length(freqBands)
% 
%     % Figure with correlation form and reference ELP:
%     
%     subplot(2,2,f)
%     % Received correlation form (synced)
%     plot(1e6*trackRes.(freqBands{f}).chipVec_del, trackRes.(freqBands{f}).corrForm_del_abs, '-ok')
%     hold on
%     % Correct locations of EPL
%     vline(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).esamp_ref), 'g')
%     vline(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).psamp_ref), 'g')
%     vline(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).lsamp_ref), 'g')
%     xlim(1e6*[trackRes.(freqBands{f}).chipVec_del(1) trackRes.(freqBands{f}).chipVec_del(end)])
%     ylabel([ freqBands{f} ' correlation power' ])
%    
%     % Figure with tracking result:
%     subplot(2,2,f+2)
%     hold on
%     
%     plot(1e6*trackRes.(freqBands{f}).chipVec_del, trackRes.(freqBands{f}).earlyLateDiff, '-or')
%     plot(1e6*trackRes.(freqBands{f}).chipVec_del, trackRes.(freqBands{f}).phaseCorrPow, '-og')
%     plot(1e6*trackRes.(freqBands{f}).chipVec_del, trackRes.(freqBands{f}).trackingOutput, '-ob')
%     
%     vline(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).earlySamples(trackRes.(freqBands{f}).bestSamp)), 'm')
%     vline(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).phaseSamples(trackRes.(freqBands{f}).bestSamp)), 'm')
%     vline(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).lateSamples(trackRes.(freqBands{f}).bestSamp)), 'm')
%     
%     scatter(1e6*trackRes.(freqBands{f}).chipVec_del(trackRes.(freqBands{f}).bestSamp), trackRes.(freqBands{f}).trackingOutput(trackRes.(freqBands{f}).bestSamp), 250, 'm')
%     
%     xlim(1e6*[trackRes.(freqBands{f}).chipVec_del(1) trackRes.(freqBands{f}).chipVec_del(end)])
%     text(-1.8, max( trackRes.(freqBands{f}).phaseCorrPow), ['MP error (m): ' num2str( trackRes.(freqBands{f}).MP_error)])
%     xlabel('Correlation delay (\mus)')
%     ylabel([ freqBands{f} ' correlation power' ])
% end
 end
%%

% 
% indToPlot = 1:1:tToPlot*f_samp;
% figure(1);
% hold on
% plot(1e6*t_vec(indToPlot),s(indToPlot),lineSpecs{sig})
% %title([sys ' ' freqBand ' ' comp  ' SV: '  num2str(satID)])
% xlabel('time [\mus]')
% ylabel('time domain signal')
% 
% 
% 
% end
% legend('GPS L1', 'GPS L5', 'GAL E1', 'GAL E5')
% 
% figure;
% histogram(abs(s))