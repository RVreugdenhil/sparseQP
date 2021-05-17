function [path_data] = get_path_data()
% return the path to the folder containing all datasets
    path1 = mfilename('fullpath');
    path2 = strsplit(path1,filesep);
    path3 = strjoin(path2(1:end-3),filesep);
    path_data = strcat(path3,filesep,'Data',filesep);
end

