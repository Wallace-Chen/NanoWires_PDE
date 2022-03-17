% a function to construct a char array as the f coefficient for the poisson
% equations of the exchange potentials.

function f = f_ex_coefficient(N)

    f = strings(N*N, 1);
    for i = 1 : N*N
        f(i) = sprintf('n_e_xc_func(x, y, %d)', i);
    end
    f = char( convertStringsToChars(f) );
    
end