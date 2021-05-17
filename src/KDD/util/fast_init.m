
try
    parpool;
catch
    disp('matlab pool already started');
end

try
    matlabpool;
catch
    disp('matlab pool already started');
end



%   getenv('NUMBER_OF_PROCESSORS')

