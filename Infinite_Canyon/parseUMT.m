function tpos = parseUMT(path)
fid = fopen(path);
tline = fgetl(fid);

samp = 0;
clear tpos
while ischar(tline)
    samp = samp + 1;
    C = strsplit(tline,',');
    tline = fgetl(fid);
    tpos(1:1:4, samp) = [str2num(C{1}), str2num(C{4}), str2num(C{5}), str2num(C{6})];
end
fclose(fid);
end