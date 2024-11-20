function [object_class_est, C, Dstar, Qstar, lambdaMstar, Ccs, Dstarcs, Qstarcs, lambdaMstarcs] = match_objects_w_knownparams(init_tracking, data, Xcs, object_class_est_prev, hyper, factorials, mcmc, ell, C, Dstar, Qstar, lambdaMstar)

assignments = associate_esttracks_w_truetracks_notmunkres(data.X, Xcs); %omega_tracks, column index corresponds to true track, row to estimated track
for i = 1:size(assignments,2)
    iprime = assignments(i);
    Ccs(i) = data.C(iprime);
end
C(:,ell) = Ccs;
Dstar(:,:, 1:max(data.C), ell) = data.Dstar; Dstarcs = Dstar(:,:,:, ell);
Qstar(:,:,1:max(data.C),ell) = data.Sigmavstar; Qstarcs = Qstar(:,:,:, ell);
lambdaMstar(1:max(data.C),ell) = data.lambdaMstar; lambdaMstarcs = lambdaMstar(:,ell);

object_class_est = object_class_est_prev;
existsteps = zeros(1, init_tracking.maxNactive);
Msum = zeros(1, init_tracking.maxNactive);
denom_factorial = zeros(1, init_tracking.maxNactive);
Zsqsum = zeros(2, init_tracking.maxNactive);
Vsqsum = zeros(4, init_tracking.maxNactive);

for i = 1:init_tracking.maxNactive
    denom_factorial_exp = 0;
    for t = init_tracking.Ts(i):init_tracking.Te(i)
        rot = calcrotmat(Xcs(:,t, i));
        premeasurements_i = rot.'* (init_tracking.Z(:, :, t, i) - Xcs(1:2, t, i));
        existsteps(i) = existsteps(i) + 1;
        Msum(i) = Msum(i) + init_tracking.M(t,i);
        denom_factorial_exp = denom_factorial_exp + log(factorials(max(min(init_tracking.M(t,i), 150),1)));
        Zsqsum(:, i) = Zsqsum(:, i) + nansum(premeasurements_i.^2, 2);
        if t == init_tracking.Ts(i)
            driving_noise_est = [0;0;0;0];
        else
            driving_noise_est = Xcs(:,t,i) - hyper.F*Xcs(:,t-1,i);
        end
        Vsqsum(:,i) = Vsqsum(:,i) + nansum(driving_noise_est.^2, 2);
    end
    denom_factorial(i) = denom_factorial_exp;
end
if init_tracking.maxNactive > 0
    object_class_est.existsteps = existsteps; object_class_est.Msum = Msum; object_class_est.Zsqsum = Zsqsum; object_class_est.Vsqsum = Vsqsum;
else
    object_class_est = object_class_est_prev;
    return
end

end