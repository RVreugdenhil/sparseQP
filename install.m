% Install function
% Add the src folder into the working directory of matlab


path1 = pwd;
path_src = strcat(path1,filesep,'src');
path(genpath(path_src),path); 
clear path1 path_src