function inverse_gamma_expos = calc_inverse_gamma_log(x, a, b)
        inverse_gamma_expos = sum(a.*log(b) - log(gamma(a)) - (a + 1).*log(x) - b./x);
end