function [steps_candidates] = get_realval(x)
steps_candidates=[];
for i=1:length(x)
    if(isreal(x(i)))
        steps_candidates=[steps_candidates;x(i)];
    end
end

steps_candidates = sort(steps_candidates(:));