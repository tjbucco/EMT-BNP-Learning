function [data, hyperparams] = generatemeasurements(data, hyperparams, script)

maxT = script.maxT;
% generate target measurements
M = NaN(1, maxT, data.maxNactive);
for n = 1:data.maxNactive
    M(1, find(data.nactive(n, :), 1, 'first'):find(data.nactive(n, :), 1, 'last'), n) = poissrnd(data.lambdaMstar(data.C(n)), 1, find(data.nactive(n, :), 1, 'last') - find(data.nactive(n, :), 1, 'first') + 1);
    
    if M(1, find(data.nactive(n, :), 1, 'first'), n) == 0 && script.MTsnz
        M(1, find(data.nactive(n, :), 1, 'first'), n) = 1;
    end
end
M(M > script.maxnm) = script.maxnm;
maxm = max(nansum(M, 3));

if script.clutterflag
Z = NaN(2, maxm, maxT, data.maxNactive + 1);
else
Z = NaN(2, maxm, maxT, data.maxNactive);
end
for n = 1:data.maxNactive
    for t = find(data.nactive(n, :), 1, 'first'):find(data.nactive(n, :), 1, 'last')
            aij = vecangle(data.X(3:4, t, n));
            Rij = [cos(aij), -sin(aij); sin(aij), cos(aij)];
            Z(:, 1:M(1, t, n), t, n) = transpose(mvnrnd(data.X(1:2, t, n), Rij*data.Dstar(:, :, data.C(n))*transpose(Rij), M(1, t, n)));
    end
end

if script.clutterflag
    % generate clutter measurements
    max_xy = [max(max(max(Z(1,:,:,1:data.maxNactive)))); max(max(max(Z(2,:,:,1:data.maxNactive))))];
    min_xy = [min(min(min(Z(1,:,:,1:data.maxNactive)))); min(min(min(Z(2,:,:,1:data.maxNactive))))]; % min and max on both axes according to target measurements
    V = norm(max_xy - min_xy); % volume of area
    lambdaclutter = hyperparams.clutter_density*V; %expectation of number of false alarms
    lambdaclutter = hyperparams.lambdaclutter;
    m_clutter = poissrnd(lambdaclutter, maxT, 1); % number of false alarms
    Z = [Z NaN(2, max(m_clutter) + maxm, maxT, data.maxNactive + 1)];

    for t = 1:maxT
        Z(1:2,1:m_clutter(t),t,data.maxNactive + 1) = min_xy + (max_xy - min_xy).*rand(2,m_clutter(t)); % generate false alarms
    end

    max_xy = [max(max(max(Z(1,:,:,1:data.maxNactive + 1))));max(max(max(Z(2,:,:,1:data.maxNactive + 1))))];
    min_xy = [min(min(min(Z(1,:,:,1:data.maxNactive + 1))));min(min(min(Z(2,:,:,1:data.maxNactive + 1))))];
else
    hyperparams.lambdaclutter = 0;
    max_xy = [max(max(max(Z(1,:,:,1:data.maxNactive)))); max(max(max(Z(2,:,:,1:data.maxNactive))))];
    min_xy = [min(min(min(Z(1,:,:,1:data.maxNactive)))); min(min(min(Z(2,:,:,1:data.maxNactive))))]; % min and max on both axes according to target measurements
    V = norm(max_xy - min_xy); % volume of area
end

data.Z = Z;
data.V = V;
data.max_xy = max_xy;
data.min_xy = min_xy;

hyperparams.clutter_gauss_approx_cov = [1/12*(max_xy(1) - min_xy(1))^2, 0; 0, 1/12*(max_xy(2) - min_xy(2))^2];

% determine max indices, number of measurements per time

M = NaN(maxT, data.maxNactive);
for n = 1:data.maxNactive
    for t = 1:maxT
        if all(isnan(Z(1, :, t, n)))
            M(t, n) = 0;
        else
            M(t, n) = find(~isnan(Z(1, :, t, n)), 1, 'last');
        end
    end
end
Msum = nansum(M, 1);

measurements_all = cell(1,maxT);
no_of_measurements = NaN(1, maxT);
for t = 1:maxT
    measurements = Z(:,:,t,:);
    measurements_all{t} = reshape(measurements(~isnan(measurements)),2,[]);
    no_of_measurements(t) = size(measurements_all{t},2);
end
data.measurements_all = measurements_all;
data.max_no_of_measurements = max(no_of_measurements);
data.M = M;
data.Msum = Msum;
data.Mclutter = m_clutter;

end
