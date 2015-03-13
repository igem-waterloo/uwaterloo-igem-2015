function [final] = chris_salahub_get_min(x)
%CHRIS_SALAHUB_GET_MIN consumes an array of floats and outputs the
%   minimum and first index
ind = 1;
curr_max = x(1);
for j = 2:length(x)
    if x(j) < curr_max
        ind = j;
        curr_max = x(j);
    end
end
final = [curr_max, ind];
end

