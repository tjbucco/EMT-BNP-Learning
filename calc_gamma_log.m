function gamma_expos = calc_gamma_log(x, a, b)
        gamma_expos = a.*log(b) - log(gamma(a)) + (a - 1).*log(x) - b*x;
end