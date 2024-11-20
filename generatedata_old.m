 function [data, hyper] = generatedata(hyper, script)

%if script.fixed_dataparams_flag
    alphadp = hyper.alphadp_data;
    al0 = hyper.al0_data;
    bl0 = hyper.bl0_data;
    D0 = hyper.D0_data;
    av0 = hyper.av0_data;
    bv0 = hyper.bv0_data;
    Sigmav0 = hyper.Sigmav0_data;
% else
%     alphadp = hyper.alphadp;
%     al0 = hyper.al0;
%     bl0 = hyper.bl0;
%     D0 = hyper.D0;
%     av0 = hyper.av0;
%     bv0 = hyper.bv0;
%     Sigmav0 = hyper.Sigmav0;
% end

maxN = script.maxN;
maxT = script.maxT;
nactive = zeros(maxN, maxT);

% set active targets
n0 = min([binornd(script.binoN, script.binop), maxN]);
if script.fixed_birth_flag
    n0 = ceil(script.binoN * script.binop);
end
if n0 == 0
    n0 = 1;
end
nactive(1:n0, 1) = 1;
maxNactive = n0;

for t = 2:maxT
%     if t <= script.birth_heuristic_memory
%         survprop = 1;
%     else
        survprop = script.ps;
%     end
    nactive(:, t) = nactive(:, t - 1).*(rand(maxN, 1) < survprop);
    nactive(:, t) = nactive(:, t - 1).*(t < 1/(1 - survprop));
    nactive(:, t) = nactive(:, t - 1).*(sum(nactive(:, 1:t - 1),2) < script.minT | (rand(maxN, 1) < survprop));
 %|| (sum(nactive(:, (t - 1 - script.birth_heuristic_memory):(t - 1))) < script.birth_heuristic_memory))
    if script.fixed_birth_flag
        nbirth = (mod(t - 1, 1/(1 - survprop))  == 0) * ceil(script.binoN * script.binop);
        nbirth = (mod(t - 1, script.birthT)  == 0) * (t < script.stopbirthT) * ceil(1 * (script.binoN));
        if t == maxT
            nbirth = 0;
        end
    else
        nbirth = min([binornd(script.binoN, script.binop), maxN - maxNactive]);
    end
    nactive((maxNactive + 1):(maxNactive + nbirth), t) = 1;
    maxNactive = maxNactive + nbirth;
end
maxTactive = find(sum(nactive, 1), 1, 'last');

data = struct('nactive', nactive);

% generate clusters
U = rand(1, maxNactive);
C = zeros(1, maxNactive);
C(1) = 1;

for n = 2:maxNactive
    Ccounts = histcounts(C);

    for t = 2:length(Ccounts)
        if U(n) <= sum(Ccounts(2:t))/(alphadp + (n - 1))
            C(n) = t - 1;
            break
        end
    end

    if C(n) == 0
        C(n) = max(C) + 1;
    end
end
if script.fixed_dataparams_flag
    C = ones(1, maxNactive);
end
data.C = C;
    
% generate cluster parameters
Dstar = NaN(2, 2, max(C));
lsqinv = gamrnd(al0(:, ones(1, max(C))), 1./bl0(:, ones(1, max(C))), 2, max(C));
lsq = 1./lsqinv;
for c = 1:max(C)
    Dstar(:, :, c) = diag(lsq(:, c));
end
if script.driving_noise_only || script.fixed_dataparams_flag
    Dstar = repmat(D0, 1, 1, max(C));
end

Sigmavstar = NaN(4, 4, max(C));
ssqinv = gamrnd(av0(:, ones(1, max(C))), 1./bv0(:, ones(1, max(C))), 2, max(C));
ssq = 1./ssqinv;
for c = 1:max(C)
    Sigmavstar(:, :, c) = diag(ssq([1, 1, 2, 2], c));
end
if script.driving_noise_only || script.fixed_dataparams_flag
    Sigmavstar = repmat(Sigmav0, 1, 1, max(C));
end

lambdaMstar = min(max(gamrnd(hyper.aM0(:, ones(1, max(C))), 1./hyper.bM0(:, ones(1, max(C))), 1, max(C)), script.minlamb), script.maxlamb);
if script.driving_noise_only
    lambdaMstar = repmat(hyper.ElambdaM, 1, max(C));
end

data.Dstar = Dstar;
data.Sigmavstar = Sigmavstar;
data.lambdaMstar = lambdaMstar;

% generate kinematic state
X = NaN(4, maxT, maxNactive);
v = NaN(4, maxT, maxNactive);
if script.traj_mode == 1
    for n = 1:maxNactive
        X(:, find(nactive(n, :), 1, 'first'), n) = transpose(mvnrnd([0, 0, 0, 0], hyper.Sigmax1));
        v(:, :, n) = transpose(mvnrnd([0, 0, 0, 0], Sigmavstar(:, :, C(n)), maxT));

        for t = (find(nactive(n, :), 1, 'first') + 1):find(nactive(n, :), 1, 'last')
            X(:, t, n) = hyper.F*X(:, t - 1, n) + v(:, t - 1, n);
        end
    end
elseif script.traj_mode==2
    %% targets are born on circumference of a circle, targets move towards center of circle
    shift = 0;
    for n = 1:maxNactive
        tnborn = find(nactive(n, :), 1, 'first');
        angles = linspace(0 + mod(tnborn,  script.birthT - 1) * pi/3, 2*pi + mod(tnborn,  script.birthT - 1) * pi/3, script.binoN + 1);
        angles = linspace(0 + shift * 2*pi/5, 2*pi + shift * 2*pi/5, script.binoN + 1);
        shift = mod(shift + 1, 10);
        %angles = linspace(0, 2*pi, script.binoN + 1);
        %angles = shift * pi/3; %* mod(shift,2);
        x = cos(angles);
        y = sin(angles);
        X(:, tnborn, n) = transpose([hyper.radius * x(mod(n - 1, script.binoN) + 1), hyper.radius * y(mod(n - 1, script.binoN) + 1), hyper.velo * - x(mod(n - 1, script.binoN) + 1), hyper.velo * - y(mod(n - 1, script.binoN) + 1)]);
        v(:, :, n) = transpose(mvnrnd([0, 0, 0, 0], Sigmavstar(:, :, C(n)), maxT));
        for t = (find(nactive(n, :), 1, 'first') + 1):find(nactive(n, :), 1, 'last')
            X(:, t, n) = hyper.F*X(:, t - 1, n) + v(:, t - 1, n);
        end
    end
else
    %% targets are born along 2 straight lines well separated, targets move perpendicular to the straight lines towards one another
    for n = 1:maxNactive
        tnborn = find(nactive(n, :), 1, 'first');
        birth_loc = linspace(-hyper.radius/4, hyper.radius/4, script.binoN);
        x = birth_loc;
        y = (-1).^(0:script.binoN); %alternate between radius and -radius
        X(:, tnborn, n) = transpose([x(mod(n - 1, script.binoN) + 1), hyper.radius * y(mod(n - 1, script.binoN) + 1), 0, hyper.velo * - y(mod(n - 1, script.binoN) + 1)]);
        v(:, :, n) = transpose(mvnrnd([0, 0, 0, 0], Sigmavstar(:, :, C(n)), maxT));
        for t = (find(nactive(n, :), 1, 'first') + 1):find(nactive(n, :), 1, 'last')
            X(:, t, n) = hyper.F*X(:, t - 1, n) + v(:, t - 1, n);
        end
    end
end
data.X = X;

% generate target measurements
M = NaN(1, maxT, maxNactive);
for n = 1:maxNactive
    M(1, find(nactive(n, :), 1, 'first'):find(nactive(n, :), 1, 'last'), n) = poissrnd(lambdaMstar(C(n)), 1, find(nactive(n, :), 1, 'last') - find(nactive(n, :), 1, 'first') + 1);

    if M(1, find(nactive(n, :), 1, 'first'), n) == 0 && script.MTsnz
        M(1, find(nactive(n, :), 1, 'first'), n) = 1;
    end
end
M(M > script.maxnm) = script.maxnm;
maxm = max(nansum(M, 3));

if script.clutterflag
    Z = NaN(2, maxm, maxT, maxNactive + 1);
else
    Z = NaN(2, maxm, maxT, maxNactive);
end
for n = 1:maxNactive
    for t = find(nactive(n, :), 1, 'first'):find(nactive(n, :), 1, 'last')
        aij = vecangle(X(3:4, t, n));
        Rij = [cos(aij), -sin(aij); sin(aij), cos(aij)];
        Z(:, 1:M(1, t, n), t, n) = transpose(mvnrnd(X(1:2, t, n), Rij*Dstar(:, :, C(n))*transpose(Rij), M(1, t, n)));
    end
end

if script.clutterflag
    % generate clutter measurements
    max_xy = [hyper.radius;hyper.radius];
    min_xy = [-hyper.radius;-hyper.radius]; % min and max on both axes according to target measurements
    V = norm(max_xy - min_xy); % volume of area
    lambdaclutter = hyper.clutter_density*V; %expectation of number of false alarms
    %hyper.lambdaclutter = lambdaclutter;
    lambdaclutter = hyper.lambdaclutter;
    m_clutter = poissrnd(lambdaclutter, maxT, 1); % number of false alarms
    Z = [Z NaN(2, max(m_clutter) + maxm, maxT, maxNactive + 1)];

    for t = 1:maxT
        Z(1:2,1:m_clutter(t),t,maxNactive + 1) = min_xy + (max_xy - min_xy).*rand(2,m_clutter(t)); % generate false alarms
    end

    max_xy = [max(max(max(Z(1,:,:,1:maxNactive + 1))));max(max(max(Z(2,:,:,1:maxNactive + 1))))];
    min_xy = [min(min(min(Z(1,:,:,1:maxNactive + 1))));min(min(min(Z(2,:,:,1:maxNactive + 1))))];
else
    hyper.lambdaclutter = 0;
    max_xy = [max(max(max(Z(1,:,:,1:maxNactive)))); max(max(max(Z(2,:,:,1:maxNactive))))];
    min_xy = [min(min(min(Z(1,:,:,1:maxNactive)))); min(min(min(Z(2,:,:,1:maxNactive))))]; % min and max on both axes according to target measurements
    V = norm(max_xy - min_xy); % volume of area
end

data.Z = Z;
data.V = V;
data.max_xy = max_xy;
data.min_xy = min_xy;

hyper.clutter_gauss_approx_cov = [1/12*(max_xy(1) - min_xy(1))^2, 0; 0, 1/12*(max_xy(2) - min_xy(2))^2];

% determine max indices, number of measurements per time, and existence
% tracks
nactive = zeros(maxNactive, maxT);
Ts = NaN(1, maxNactive);
Te = NaN(1, maxNactive);
for n = 1:maxNactive
    nactive(n, :) = ~isnan(X(1, :, n));
    Ts(1, n) = find(nactive(n, :), 1, 'first');
    Te(1, n) = find(nactive(n, :), 1, 'last');
end

M = NaN(maxT, maxNactive);
for n = 1:maxNactive
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
data.maxTactive = maxTactive;
data.maxNactive = maxNactive;
data.nactive = nactive;
data.Ts = Ts;
data.Te = Te;
data.M = M;
data.Msum = Msum;
data.Mclutter = m_clutter;
end
