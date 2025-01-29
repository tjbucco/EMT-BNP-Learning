function probs = calc_indicator_probs(Q, D, lambda, Sigmax1, Zsqsum, no_of_targets, Ccs, Cwoi, hyper, T, no_of_measurements, Vsqsum, denom_expo_factorial, Xinit, premeasurements, driving_noises, no_of_measurements_k)

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
driving_noises_1 = reshape(driving_noises(1:2,2:end),1, []);
driving_noises_2 = reshape(driving_noises(3:4,2:end),1, []);


aMstarc = hyper.aM0 + no_of_measurements;
bMstarc = hyper.bM0 + T + 1;
lambda = [gamrnd(aMstarc, 1./bMstarc); lambda];

[no_of_meas_uniq, ~, i_unq] = unique(no_of_measurements_k(1:end-1));

initial_model = mvnpdfln(Xinit.', 0, diag(Sigmax1).');

aexpo = NaN(1,max(Ccs) + 1);
for c = 1:max(Ccs) + 1    
    if c == 1
        
        evid_A_1 = mvtpdfln_fast(driving_noises_1, hyper.bv0(1)/hyper.av0(1) * eye(size(driving_noises_1,2)), 2*hyper.av0(1));
        evid_A_2 = mvtpdfln_fast(driving_noises_2, hyper.bv0(2)/hyper.av0(2) * eye(size(driving_noises_2,2)), 2*hyper.av0(2));
        evid_B_1 = mvtpdfln_fast(premeasurements(1,:), hyper.bl0(1)/hyper.al0(1) * eye(size(premeasurements,2)), 2*hyper.al0(1));
        evid_B_2 = mvtpdfln_fast(premeasurements(2,:), hyper.bl0(2)/hyper.al0(2) * eye(size(premeasurements,2)), 2*hyper.al0(2));
        evid_C_1 = negmnpdfln(no_of_measurements_k.', 1/(size(no_of_measurements_k,1) + hyper.bM0) * ones(1, size(no_of_measurements_k,1) + 1), hyper.aM0);

        aexpo(1, 0 + 1) = log(hyper.alphadp/(hyper.alphadp + no_of_targets - 1)) + evid_A_1 +  evid_A_2 + evid_B_1 + evid_B_2 + evid_C_1 + initial_model;
    else
        
        motion_model_A1 = mvnpdfln(driving_noises_1, 0, Q(1,1,c)* ones(1,size(driving_noises_1,2)));
        motion_model_A2 = mvnpdfln(driving_noises_2, 0, Q(3,3,c)* ones(1,size(driving_noises_2,2)));

        measurement_model_B1 = mvnpdfln(premeasurements(1,:), 0, D(1,1,c)*ones(1, size(premeasurements,2)));
        measurement_model_B2 = mvnpdfln(premeasurements(2,:), 0, D(2,2,c)*ones(1, size(premeasurements,2)));

        poisspdf_uniq = poisspdfln(no_of_meas_uniq, lambda(c));
        measurement_model_C1 = sum(poisspdf_uniq(i_unq));

        motion_model = motion_model_A1 + motion_model_A2;
        measurement_model = measurement_model_B1 + measurement_model_B2 + measurement_model_C1;

        aexpo(c) = log(sum(Cwoi == (c - 1))) - log(hyper.alphadp + no_of_targets - 1) + initial_model + measurement_model + motion_model;
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