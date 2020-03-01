function codes_upsampled = codeGenerator(sys, freqBand, comp, satID, f_samp, const)
%CODEGENERATOR Summary of this function goes here
%   Detailed explanation goes here

if strcmp(sys, 'GPS')
    if strcmp(freqBand, 'L1')
        code = cacode(satID);
        code(code==0) = -1;
        codes.(sys).(freqBand).CA = code;
    elseif strcmp(freqBand, 'L2')
        if strcmp(comp, 'CL')
            codes.(sys).(freqBand).CL = codes_L2CL(:,satID)';
            load('codes_L2CL.mat')
        elseif strcmp(comp, 'CM')
            load('codes_L2CM.mat')
            codes.(sys).(freqBand).CM = codes_L2CM(:,satID)';
        end

    elseif strcmp(freqBand, 'L5')  
        load('codes_L5I.mat')
        load('codes_L5Q.mat')
        codes.(sys).(freqBand).I = codes_L5I(:,satID)';
        codes.(sys).(freqBand).Q = codes_L5Q(:,satID)';
    end

elseif strcmp(sys, 'GAL')
    if strcmp(freqBand, 'E1')
        load('codes_E1B.mat')
        load('codes_E1C.mat')
        codes.(sys).(freqBand).B = codes_E1B(:,satID)';
        codes.(sys).(freqBand).C = codes_E1C(:,satID)';
    elseif strcmp(freqBand, 'E5')
        load('codes_E5aI.mat')
        load('codes_E5aQ.mat')
        load('codes_E5bI.mat')
        load('codes_E5bQ.mat')
        codes.(sys).(freqBand).aI = codes_E5aI(:,satID)';
        codes.(sys).(freqBand).aQ = codes_E5aQ(:,satID)';
        codes.(sys).(freqBand).bI = codes_E5bI(:,satID)';
        codes.(sys).(freqBand).bQ = codes_E5bQ(:,satID)';
    end
end

allFields = fields(codes.(sys).(freqBand));
for f = 1:1:length(allFields)
    currCode = codes.(sys).(freqBand).(allFields{f});
    if strcmp(sys, 'GAL')
        setLength = const.(sys).(freqBand).OS.T_code*f_samp;
    elseif strcmp(sys, 'GPS') & strcmp(freqBand, 'L5')
        const.(sys).(freqBand).SOL.T_code;     
        setLength = const.(sys).(freqBand).SOL.T_code*f_samp;
    else
        setLength = const.(sys).(freqBand).(allFields{f}).T_code*f_samp;
    end
    factor = setLength/size(currCode, 2);
    if factor > 1
        code_upsampled = upsampleCode(currCode , factor);
    else
        code_upsampled = currCode;
    end
    codes_upsampled.(sys).(freqBand).(allFields{f}) = code_upsampled;
    
end

end

