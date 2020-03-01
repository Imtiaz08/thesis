function out = upsampleCode(in, factor)
%UPSAMPLE Summary of this function goes here
%   Detailed explanation goes here
%upsample to desired rate

L = length(in);
if factor~=1
	%fractional upsampling with zero order hold
	index=0;
	for cnt = 1/factor:1/factor:L
		index=index+1;
		if ceil(cnt) > L   %traps a floating point error in index
			gfs(:,index)=in(:,L);
		else
			gfs(:,index)=in(:,ceil(cnt));
		end
	end 
	out=gfs;
end

end

