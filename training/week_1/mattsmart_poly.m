function [y] = mattsmart_poly(x, coeffs)
    
sum_recurse = 0;
coeffs_recurse = coeffs;

function sum = poly_recurse(coeffs_recurse, sum_recurse)
    n = length(coeffs_recurse);
    if n == 1
        sum = sum_recurse + coeffs_recurse(n);
    else
        sum_recurse = sum_recurse + coeffs_recurse(n)*x.^(n-1);
        coeffs_recurse(n) = [];  % remove the last element
        sum = poly_recurse(coeffs_recurse, sum_recurse);
    end
    
end

y = poly_recurse(coeffs_recurse, sum_recurse);

end

