function probs = calc_indicator_probs_seq(Q, D, Sigmax1, Zsqsum, no_of_targets, Ccs, Cwoi, hyper, T, no_of_measurements, Vsqsum, initial_ind, X)

if hyper.alphadp == inf
    probs = zeros(1,max(Ccs) + 1);
    probs(1) = 1;
return
end

alstarc = hyper.al0 + 1/2*no_of_measurements*[1; 1];
blstarc = hyper.bl0 + 1/2*Zsqsum;
lsqstarcinv = gamrnd(alstarc, 1./blstarc, 2, 1);
lsqstarc = 1./lsqstarcinv;
D = cat(3,diag(lsqstarc), D);

avstarc = hyper.av0 + 1/2*(T)*[2; 2];
Vsqsum_pos = sum(Vsqsum(1:2));
Vsqsum_velo = sum(Vsqsum(3:4));
bvstarc = hyper.bv0 + 1/2*[Vsqsum_pos; Vsqsum_velo];
ssqstarcinv = gamrnd(avstarc, 1./bvstarc, 2, 1);
ssqstarc = 1./ssqstarcinv;
Q = cat(3,diag(ssqstarc([1, 1, 2, 2])), Q);
% diag(Vsqsum/(T-1))
% Q(end + 1) = diag(Vsqsum/(T-1));

aexpo = NaN(1,max(Ccs) + 1);
for c = 1:max(Ccs) + 1
    initial_model = 0;
    motion_model = 0;
    %Ynksqsum = nansum((transpose(Rkt)*(Z(:, :, t) - X(1:2, k*ones(1, maxnm), t))).^2, 2);
    %measurement_model = measurement_model + log((lambda(c)^no_of_measurements(t))/factorials(min(no_of_measurements(t), 150))) - lambda(c) - 1/2*(1./transpose(diag(D(:,:,c)))*Yntsqsum);
    measurement_model = 1/2*(1./transpose(diag(D(:,:,c)))*Zsqsum);
    if initial_ind
        initial_model = - 1/2*(1./transpose(diag(Sigmax1))*X.^2);
    else
        motion_model = motion_model - 1/2 * 1./transpose(diag(Q(:,:,c))) * Vsqsum;
    end
    
    if c == 1
        alt = hyper.al0 + 1/2*no_of_measurements;
        blt = hyper.bl0 + 1/2*Zsqsum;
        extent_conditional = calc_inverse_gamma_log(diag(D(:,:,c)), alt, blt);
        
        avt = hyper.av0 + T;
        bvt = hyper.bv0 + 1/2*[Vsqsum_pos; Vsqsum_velo];
        driving_conditional = calc_inverse_gamma_log([Q(1,1,c);Q(3,3,c)], avt, bvt);
        
        conditional_probs = extent_conditional + driving_conditional;
        
        base_pdf = calc_inverse_gamma_log(diag(D(:,:,c)), hyper.al0, hyper.bl0) + calc_inverse_gamma_log([Q(1,1,c);Q(3,3,c)], hyper.av0, hyper.bv0);
        
        aexpo(1, 0 + 1) = log(hyper.alphadp/(hyper.alphadp + no_of_targets - 1)) + base_pdf - conditional_probs + initial_model + measurement_model + motion_model;
        %aexpo(1, 0 + 1) = log(class_hyperparameters.alphadp/(class_hyperparameters.alphadp + K - 1)) + base_pdf - conditional_probs + initial_model + measurement_model + motion_model;
        %aexpo(1, 0 + 1) = log(alphadp*(2*pi)^(2/3)*(hyper.bl0(1).^hyper.al0(1)*hyper.bl0(2).^hyper.al0(2)*bM0^aM0)/(gamma(hyper.al0(1))*gamma(hyper.al0(2))*gamma(aM0))) + aln(1)*log(aln(1)/bln(1)) + aln(2)*log(aln(2)/bln(2)) + aMn*log(aMn/bMn) - 1/2*log(aln(1)*aln(2)*aMn) - (aln(1) + aln(2) + aMn);
    else
        aexpo(c) = log(sum(Cwoi == (c - 1))) - log(hyper.alphadp + no_of_targets - 1) + initial_model + measurement_model + motion_model;
        %aexpo(c + 1) = log(sum(Cwon == c)) - log(class_hyperparameters.alphadp + no_of_targets - 1) + initial_model + measurement_model + motion_model;
        %1/(class_hyperparameters.alphadp + maxN - 1)*(sum(Cwon == c))*(T products of measurement model)*(target state prior)*(T - 1 products of state transition model);
        %aexpo(c + 1) = log(sum(Cwon == c)) + Msum(1, n)*log(lambdaMstarcs(1, c)) - (Te(n) - Ts(n) + 1)*lambdaMstarcs(1, c) - 1/2*Msum(1, n)*log(Dstarcs(1, 1, c)*Dstarcs(2, 2, c)) - 1/2*(Zsqsum(1)/Dstarcs(1, 1, c) + Zsqsum(2)/Dstarcs(2, 2, c));
    end
end
if max(aexpo) == Inf
    a = isinf(aexpo);
    a = a/sum(a);
else
    aexpo = aexpo - max(aexpo);
    a = exp(aexpo);
    a = a/sum(a);
end
probs = a;