
function [min_index] = minfcn(A)
%1 - return the minimum element of array
%2 - return index

N = numel(A);
min = A(1);

for i = 1:1:N
    if A(i) < min
        min = A(i);
        index = find(A == A(i));
        
    end
end

min_index = 1;
min_index(1) = min;
min_index(2) = index;


end

