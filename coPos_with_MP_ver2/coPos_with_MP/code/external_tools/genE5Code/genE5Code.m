function [refE5a,refE5b] = genE5Code(sv,fs)
% function [refE5a,refE5b] = genE5Code(sv,fs)
% Generates 10230 length E5a and E5b codes for Galileo 
%
% 
% sv: a row or column vector of the SV's to be generated
% 	 valid entries are 1 to 50
% fs: optional number of samples per chip (defaults to 1), fractional samples allowed, must be 1 or greater.
%
% For multiple samples per chip, function is a zero order hold. 
%
%
% For example to generate the E5a and E5b codes for PRN 1:
% [refE5a,refE5b] = E5abcode(1),
% to generate the E5a and E5b codes for PRN 1 and PRN 7 use:
% [refE5a,refE5b] = E5abcode([1,7])
%
% Upsample made by Dan Boschen boschen@loglin.com 
%
% For more information refer to the "Galileo Open Service Signal In Space Interface Control Document (OS SIS ICD)"
% http://www.gsc-europa.eu/gnss-markets/segments-applications/os-sis-icd
%
% Konstantin Veprev 
% konstantin.veprev@ntlab.com

% Revision History
% rev 1.0 Konstantin Veprev 16-03-2016  Initial Release

%%

L=10230;

if nargin<2
	fs=1;
end


if (max(sv)>50) || (min(sv)<1) || (min(size(sv))~=1)
	error('sv must be a row or column vector with integers between 1 and 37\n')
end

if fs<1
	error('fs must be 1 or greater\n')
end	

% force integers
testint=round(sv)-sv;
if testint ~= 0 
	warning('non-integer value entered for sv, rounding to closest integer\n');
	sv = round(sv);
end

%%
% E5 Primary Codes Register 1 and Register 2 Values 
feedBackReg1E5aI = 40503;
feedBackReg2E5aI = 50661;
feedBackReg1E5aQ = 40503;
feedBackReg2E5aQ = 50661;

feedBackReg1E5bI = 64021;
feedBackReg2E5bI = 51445;
feedBackReg1E5bQ = 64021;
feedBackReg2E5bQ = 43143;

% Base Register 2 Start Values
[loadReg2E5aI,loadReg2E5aQ,loadReg2E5bI,loadReg2E5bQ] = startvalues(sv);

%%

refE5a = complex(genGoldCode(feedBackReg1E5aI,feedBackReg2E5aI,loadReg2E5aI),genGoldCode(feedBackReg1E5aQ,feedBackReg2E5aQ,loadReg2E5aQ)).';
refE5b = complex(genGoldCode(feedBackReg1E5bI,feedBackReg2E5bI,loadReg2E5bI),genGoldCode(feedBackReg1E5bQ,feedBackReg2E5bQ,loadReg2E5bQ)).';

%%

%upsample to desired rate
if fs~=1
	%fractional upsampling with zero order hold
	index=0;
	for cnt = 1/fs:1/fs:L
		index=index+1;
		if ceil(cnt) > L   %traps a floating point error in index
			gfsa(:,index)=refE5a(:,L);
            gfsb(:,index)=refE5b(:,L);

		else
			gfsa(:,index)=refE5a(:,ceil(cnt));
            gfsb(:,index)=refE5b(:,ceil(cnt));

		end
	end 
	refE5a=gfsa;
    refE5b=gfsb;

end
