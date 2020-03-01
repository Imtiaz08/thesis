function [s, codeToCorr, all_theta_MP] = signalGenerator(const, sys,...
    freqBand, comp, code, t, f_samp, w_phase, w_amp, time_delay_MP)
%SIGNALGENERATOR Summary of this function goes here
%   Detailed explanation goes here
samples_delay_MP = round(time_delay_MP*f_samp);

D = 1;
if strcmp(sys, 'GPS')
    if strcmp(freqBand, 'L1')
       theta_MP_carrier = 2*pi*mod(const.GPS.L1.CA.f_carrier*time_delay_MP, 1);
       %  theta_MP_carrier = 0;
        carrier = cos(2*pi*const.GPS.L1.CA.f_carrier*(t+0) + w_phase*randn(1,length(t)) + theta_MP_carrier);
        C = circshift(code.(sys).(freqBand).CA, samples_delay_MP);
        s = C.*D.*carrier + w_amp*randn(1, length(t));
        codeToCorr = code.(sys).(freqBand).CA;
        all_theta_MP = theta_MP_carrier;
    elseif strcmp(freqBand, 'L2') % Problem: L2CL duration is 1.5 seconds - very heavy to simulate!! 
         
        theta_MP_carrier = 2*pi*mod(const.GPS.L2.CL.f_carrier*time_delay_MP, 1);
        %theta_MP_carrier = 0;
        carrier = cos(2*pi*const.GPS.L2.CM.f_carrier*(t+0) + w_phase*randn(1,length(t)) + theta_MP_carrier);
        C = circshift(code.(sys).(freqBand).CM, samples_delay_MP);
        s = C.*D.*carrier + w_amp*randn(1, length(t));
        codeToCorr = code.(sys).(freqBand).CM;
        all_theta_MP = theta_MP_carrier;
    elseif strcmp(freqBand, 'L5')
        %theta_MP_carrier = 2*pi*mod(const.GPS.L5.SOL.f_carrier*time_delay_MP, 1);
        theta_MP_carrier = 0;
        carrier = cos(2*pi*const.GPS.L5.SOL.f_carrier*t + w_phase*randn(1,length(t)) + theta_MP_carrier);
        code_L5_I = circshift(code.(sys).(freqBand).I, samples_delay_MP);
        code_L5_Q = circshift(code.(sys).(freqBand).Q, samples_delay_MP);
        xL5 = (code_L5_I.*D + 1j*code_L5_Q);
        s = xL5.*carrier;
        codeToCorr = code.(sys).(freqBand).I;
        all_theta_MP = theta_MP_carrier;
    end
elseif strcmp(sys, 'GAL')
    if strcmp(freqBand, 'E1')
        f_carrier = const.GAL.E1.OS.f_carrier;
        f_subcarrier = const.GAL.E1.OS.f_subcarrier;
       % theta_MP_carrier = 2*pi*mod(f_carrier*time_delay_MP, 1);
       theta_MP_carrier = 0;
        code_E1B = circshift(code.(sys).(freqBand).B, samples_delay_MP);
        code_E1C = circshift(code.(sys).(freqBand).C, samples_delay_MP);
        alpha = const.GAL.E1.OS.alpha;
        beta = const.GAL.E1.OS.beta;
        m = const.GAL.E1.OS.m; % 6
        n = const.GAL.E1.OS.n; % 1       
        theta_MP_sc1 = 2*pi*mod(1*f_subcarrier*time_delay_MP, 1);
        theta_MP_sc6 = 2*pi*mod(6*f_subcarrier*time_delay_MP, 1);
        sc_E1_Ba = sign(sin(2*pi* n*f_subcarrier*t + w_phase*randn(1,length(t)) + theta_MP_sc1 ));
        sc_E1_Bb = sign(sin(2*pi* m*f_subcarrier*t + w_phase*randn(1,length(t)) + theta_MP_sc6));
        Bcomp = (code_E1B.*D).*(alpha*sc_E1_Ba + beta*sc_E1_Bb);
        Ccomp = code_E1C.*1.*(alpha*sc_E1_Ba - beta*sc_E1_Bb);
        s = ((Bcomp - Ccomp).*cos(2*pi*f_carrier*t + ...
            w_phase*randn(1, length(t)) + theta_MP_carrier)) +...
            w_amp*randn(1, length(t));
        
        sc = sqrt(10/11).*sign(sin(2*pi* 1*const.GAL.E1.OS.f_subcarrier*t)) + sqrt(1/11).*sign(sin(2*pi* 6*const.GAL.E1.OS.f_subcarrier*t));
        codeToCorr = code.(sys).(freqBand).B.*sc;

        all_theta_MP = [theta_MP_carrier, theta_MP_sc1, theta_MP_sc6];
        
    elseif strcmp(freqBand, 'E5')
        %theta_MP_carrier = 2*pi*mod(const.GAL.E5.OS.f_carrier*time_delay_MP, 1);
        theta_MP_carrier = 0;
        carrier = cos(2*pi*const.GAL.E5.OS.f_carrier*t + w_phase*randn(1,length(t)) + theta_MP_carrier);
        tbE5 = 1/const.GAL.E5.OS.f_subcarrier;
        dE5aI = 1; % Data bits simulated to always 1
        dE5bI = 1;
        eE5aI = circshift(code.(sys).(freqBand).aI, samples_delay_MP) .* dE5aI;
        eE5aQ = circshift(code.(sys).(freqBand).aQ, samples_delay_MP);
        eE5bI = circshift(code.(sys).(freqBand).bI, samples_delay_MP) .* dE5bI;
        eE5bQ = circshift(code.(sys).(freqBand).bQ, samples_delay_MP);
        
        econjE5aI = eE5aQ .* eE5bI .* eE5bQ;
        econjE5aQ = eE5aI .* eE5bI .* eE5bQ;
        econjE5bI = eE5bQ .* eE5aI .* eE5aQ;
        econjE5bQ = eE5bI .* eE5aI .* eE5aQ;

        AS = 1/2 * [sqrt(2)+1, 1, -1, -sqrt(2)-1, -sqrt(2)-1, -1, 1, sqrt(2)+1];
        AP = 1/2 * [-sqrt(2)+1, 1, -1, sqrt(2)-1, sqrt(2)-1, -1, 1, -sqrt(2)+1];
        
        scE5S_1per = upsampleCode(AS, tbE5*f_samp/8); % /8
        scE5P_1per = upsampleCode(AP, tbE5*f_samp/8); % /8
       
        delay_samples = f_samp*tbE5/4;


        
        
        scE5S_1per_del = circshift(scE5S_1per, delay_samples);
        scE5P_1per_del = circshift(scE5P_1per, delay_samples);
%         
%         figure;
%         plot((1:1:length(b1_base))/f_samp, b1_base)
%         hold on
%         plot((1:1:length(b2_base))/f_samp, b2_base)
%         plot((1:1:length(b1_base_del))/f_samp, b1_base_del)
%         plot((1:1:length(b2_base_del))/f_samp, b2_base_del)
        
        factor = floor(length(t)/length(scE5S_1per));
        rest = length(t) - factor*length(scE5S_1per);
        scE5S = [repmat(scE5S_1per, 1,factor) scE5S_1per(1:1:rest)];
        scE5S_del = [repmat(scE5S_1per_del, 1,factor) scE5S_1per_del(1:1:rest)];
        scE5P = [repmat(scE5P_1per, 1,factor) scE5P_1per(1:1:rest)];
        scE5P_del = [repmat(scE5P_1per_del, 1,factor) scE5P_1per_del(1:1:rest)];
         
        xE5 = 1/(2*sqrt(2)) * (eE5aI + 1j* eE5aQ) .* (scE5S - 1j*scE5S_del) + ...
            1/(2*sqrt(2)) * (eE5bI + 1j* eE5bQ) .* (scE5S + 1j*scE5S_del) + ...
            1/(2*sqrt(2)) * (econjE5aI + 1j* econjE5aQ) .* (scE5P - 1j*scE5P_del) + ...
            1/(2*sqrt(2)) * (econjE5bI + 1j* econjE5bQ) .* (scE5P + 1j* scE5P_del) ;
        
%           xE5 = 1/(2*sqrt(2)) * (eE5aI + 1j* eE5aQ) .* (scE5S - scE5S_del) + ...
%             1/(2*sqrt(2)) * (eE5bI + 1j* eE5bQ) .* (scE5S + scE5S_del) + ...
%             1/(2*sqrt(2)) * (econjE5aI + 1j* econjE5aQ) .* (scE5P - scE5P_del) + ...
%             1/(2*sqrt(2)) * (econjE5bI + 1j* econjE5bQ) .* (scE5P - scE5P_del) ;
        
        s = xE5.*carrier;
        codeToCorr = 1/(2*sqrt(2))*code.(sys).(freqBand).aI.*(scE5S - 1j*scE5S_del);
        % codeToCorr = 1/(2*sqrt(2))*code.(sys).(freqBand).aI.*(scE5S - scE5S_del);
        all_theta_MP = theta_MP_carrier;
      
    end
end

end

