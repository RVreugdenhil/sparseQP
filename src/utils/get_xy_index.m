function [xindex, yindex, khat] = get_xy_index(filename)
% return the indices of x and y from the data
% Input:
% filename: a string containing the file name of the dataset

% Output:
% xindex, yindex: index for x and y
% khat: targeted number of subspaces

    if strcmp(filename, 'Crime.csv')
        % column 28 in this file has NaN
        xindex = [1:27,29:101];
        yindex = 102;
        khat = 53;
    elseif strcmp(filename, 'Facebook.csv')
        % column 38 in this file has NaN
        xindex = [1:37, 39:53];
        yindex = 54;
        khat = 14;
    elseif strcmp(filename, 'OnlineNewsPopularity.csv')
        % column 1,2 are non predictive
        xindex = 3:60;
        yindex = 61;
        khat = 21;
    elseif strcmp(filename, 'UJIndoor.csv')
        % left out colums have constant value throughout every row
        xindex = [1:2, 5:91, 96:151, 153:157, 161:214, 216, 218:225, 228:237, 248:253, 255:292,294,295,297:300,302,305,306,308:332,334:348,350:352,354:359,361:364,366:415,417, 418, 420:422, 424:428,430:432,434:437,439,440,443, 446:450,452:457,459:481,483, 484, 486, 489, 490, 492:496,498:519];
        yindex = 521;
        khat = 49;
    elseif strcmp(filename, 'Superconductivity.csv')
        xindex = 1:81;
        yindex = 82;
        khat = 9;
    end
end

