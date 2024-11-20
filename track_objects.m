function [state_est_up] = track_objects(daoed_est, object_class_est, state_est_prev, measurements, data, hyper, script)

X = state_est_prev.Xest;
Sigmax = state_est_prev.Sigmax;
rotangle = state_est_prev.rotangle;
Sigmarot = state_est_prev.Sigmarot;
driving_noise_est = state_est_prev.driving_noise_est;

C = object_class_est.Cest;

for iprime = find(daoed_est.exist_ind_all)
    %if daoed_est.exist_indt(i)
    %i = iprime - find(daoed_est.exist_ind_all, 1, "first") + 1;
    i = iprime;
    measurements_i = measurements(:, daoed_est.DA == i);
    if i > daoed_est.total_no_of_objects - daoed_est.no_of_newobjectst
        [X(:,i), Sigmax(:, :, i), rotangle(:,i), Sigmarot(:,:,i)] = ekf_update(state_est_prev.Xest(:,i), measurements_i, state_est_prev.Sigmax(:, :, i), object_class_est.Dest(:,:,C(iprime)), state_est_prev.rotangle(:, i), state_est_prev.Sigmarot(:, :, i), object_class_est.Qest(:,:,C(iprime)), hyper.F, hyper.H);
    else
        [X(:,i), Sigmax(:, :, i), rotangle(:,i), Sigmarot(:,:,i)] = ekf_update(state_est_prev.Xest(:,i), measurements_i, state_est_prev.Sigmax(:, :, i), object_class_est.Dest(:,:,C(iprime)), state_est_prev.rotangle(:, i), state_est_prev.Sigmarot(:, :, i), object_class_est.Qest(:,:,C(iprime)), hyper.F, hyper.H);
    end
    %     else
    %         X(:,i) = NaN(); Sigmax(:, :, i) = NaN(); rotangle(:,i) = NaN(); Sigmarot(:,:,i) = NaN();
    %end
    driving_noise_est(:, i) = X(:,i) - hyper.F*state_est_prev.Xest(:,i);
end
% for i = find(not(daoed_est.exist_indt))
%     X(:,i) = NaN; Sigmax(:, :, i) = NaN; rotangle(:,i) = NaN; Sigmarot(:,:,i) = NaN; driving_noise_est(:, i) = NaN;
% end

state_est_up = struct('Xest', X, 'Sigmax', Sigmax, 'rotangle', rotangle, 'Sigmarot', Sigmarot, 'driving_noise_est', driving_noise_est);
end