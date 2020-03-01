function material = UCMP_getMaterial(type)
%GETMATERIAL Summary of this function goes here
%   Detailed explanation goes here
if strcmp(type, 'build')
    material = 'concrete';%'concrete'; % Glass could also be incorporated here, looks like reflections get slighly stronger
    %material = 'glass';
elseif strcmp(type, 'ground')
    material = 'dryGround'; 
elseif strcmp(type, 'car')
    material = 'steel';
end

end

