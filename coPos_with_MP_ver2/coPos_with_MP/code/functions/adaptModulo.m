function [allrho, allcp] = adaptModulo(rawrho, rawcp, nvals, mod_m, allFreq)
%ADAPTMODULO Summary of this function goes here
%   Detailed explanation goes here

%
% nvals = -5:1:5;
% mod_m = const.c * 20e-3;
% allFreq = {'f1', 'f2'}





for f = 1:1:length(allFreq)
    currF = allFreq{f};
    
    targetInts = mode(floor(rawrho.(currF) / mod_m), 2); % targetInt is the lower value, targetInt +1 the higher value
    intsToAdd = floor(rawrho.(currF)/ mod_m) ;
    allrho.(currF) = rawrho.(currF) + (( targetInts - intsToAdd )* mod_m);
    
    for prn_test = 1:1:32
        
        modjumps = find(abs(diff(allrho.(currF)( prn_test,:) -mod_m )) > mod_m /6); % indexes for where the signal jumps between two modulos
        direcs = sign(allrho.(currF)( prn_test,modjumps+1) - allrho.(currF) ( prn_test,modjumps) );
        
        for m = 1:1:length(modjumps)
            allrho.(currF) (prn_test,1:modjumps(m)) = allrho.(currF)(prn_test,1:modjumps(m)) + direcs(m)*mod_m;
        end
        
        % Match to best n:
        levDiffs = NaN*ones(1, length(nvals));
        for n = 1:1:length(nvals)
            if f == 1
                refMatch = rawrho.(currF)(prn_test,~isnan(rawrho.(currF)(prn_test,:))) ;
            else
                refMatch = rawrho.f1(prn_test,~isnan(rawrho.(currF)(prn_test,:))) ;
            end
            toMatch = (allrho.(currF)(prn_test,~isnan(allrho.(currF)(prn_test,:))) + nvals(n) * mod_m) ;
            levDiffs(n) = sum( abs(refMatch-toMatch) > mod_m/2) / length(refMatch);
        end
        [~, bestn] = min(levDiffs);
        allrho.(currF)(prn_test,:) = allrho.(currF) (prn_test,:) + nvals(bestn) * mod_m;
        %         figure(prn_test);
        %         hold on
        %         plot(allrho.(currF)(prn_test,:), 'marker', verMarkers{1}, 'color', freqCols{f}, 'linestyle',verLineStyles{1})
        %         plot(rawrho.(currF)(prn_test,:), 'marker', verMarkers{2}, 'color', freqCols{f}, 'linestyle',verLineStyles{2})
        %         hline( [-3:1:3]* mod_m, 'k')
        %         title(['PRN: ' num2str(prn_test)])
    end
end

for f = 1:1:length(allFreq)
    currF = allFreq{f};
    
    targetInts = mode(floor(rawcp.(currF) / mod_m), 2); % targetInt is the lower value, targetInt +1 the higher value
    intsToAdd = floor(rawcp.(currF)/ mod_m) ;
    allcp.(currF) = rawcp.(currF) + (( targetInts - intsToAdd )* mod_m);
    
    for prn_test = 1:1:32
        
        modjumps = find(abs(diff(allcp.(currF)( prn_test,:) -mod_m )) > mod_m /6); % indexes for where the signal jumps between two modulos
        direcs = sign(allcp.(currF)( prn_test,modjumps+1) - allcp.(currF) ( prn_test,modjumps) );
        
        for m = 1:1:length(modjumps)
            allcp.(currF) (prn_test,1:modjumps(m)) = allcp.(currF)(prn_test,1:modjumps(m)) + direcs(m)*mod_m;
        end
        
        % Match to best n:
        levDiffs = NaN*ones(1, length(nvals));
        for n = 1:1:length(nvals)
            if f == 1
                refMatch = rawcp.(currF)(prn_test,~isnan(rawcp.(currF)(prn_test,:))) ;
            else
                refMatch = rawcp.f1(prn_test,~isnan(rawcp.(currF)(prn_test,:))) ;
            end
            toMatch = (allcp.(currF)(prn_test,~isnan(allcp.(currF)(prn_test,:))) + nvals(n) * mod_m) ;
            levDiffs(n) = sum( abs(refMatch-toMatch) > mod_m/2) / length(refMatch);
        end
        [~, bestn] = min(levDiffs);
        allcp.(currF)(prn_test,:) = allcp.(currF) (prn_test,:) + nvals(bestn) * mod_m;
        %         figure(prn_test);
        %         hold on
        %         plot(allcp.(currF)(prn_test,:), 'marker', verMarkers{1}, 'color', freqCols{f}, 'linestyle',verLineStyles{1})
        %         plot(rawcp.(currF)(prn_test,:), 'marker', verMarkers{2}, 'color', freqCols{f}, 'linestyle',verLineStyles{2})
        %         hline( [-3:1:3]* mod_m, 'k')
        %         title(['PRN: ' num2str(prn_test)])
    end
end

end

