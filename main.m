loadseed_flag = 1;
if loadseed_flag
    load(strcat("data\srng.mat"));
else
   srng = rng;
end
rng(srng);

%close all
clearvars -except srng

loaddata_flag = 1;
if loaddata_flag
    matfiles = dir("data\*.mat");
    file_inc = 0;
    loadsi = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for data model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fixed_dataparams_flag = 0; %make one class with mean parameters

alphadp_data = 0.5; % DP concentration parameter, set high for no clustering

%driving noise covariance stuff
Essq_data = [100; 50];
Sigmav0_data = diag(Essq_data([1, 1, 2, 2])); % prior expectation of driving noise covariance matrix, increase for MORE variance
alpha_IG_Sig = 6;
av0_data = [alpha_IG_Sig; alpha_IG_Sig]; % driving noise covariance matrix parameters, increase for LESS variance
bv0_data = Essq_data.*(av0_data - 1); % driving noise covariance matrix parameters
var_Sigmav0 = Essq_data.^2./((av0_data -1).^2 .* (av0_data - 2));

%object extent stuff
Elsq_data = [675; 400];
alpha_IG_D = 12;
D0_data = diag(Elsq_data); % prior expectation of extent matrix, increase for MORE variance
al0_data = [alpha_IG_D; alpha_IG_D];% object extent covariance matrix parameters, increase for LESS variance
bl0_data = Elsq_data.*(al0_data - 1); % object extent covariance matrix parameters
var_D0 = Elsq_data.^2./((al0_data -1).^2 .* (al0_data - 2));

%mean number of measurements stuff
ElambdaM_data = 6; % prior expectation of number of measurements, increase for MORE variance
aM0_data = 4; % number of measurements parameters, increase for LESS variance
bM0_data = aM0_data/ElambdaM_data; % number of measurements parameters
var_lambdaM = aM0_data/(bM0_data^2);
minlamb = 0;
maxlamb = 35;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for Bayesian model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%class structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphadp = 0.5; %alphadp_data;%0.5; % DP concentration parameter, set high for no clustering

%driving noise covariance stuff
Essq = Essq_data .* 1; %[100; 10];
Sigmav0 = diag(Essq([1, 1, 2, 2])); % prior expectation of driving noise covariance matrix, increase for MORE variance
av0 = av0_data;%./av0_data*2; %[alpha_IG_Sig; alpha_IG_Sig]; % driving noise covariance matrix parameters, increase for LESS variance
bv0 = Essq.*(av0 - 1); % driving noise covariance matrix parameters

%object extent stuff
Elsq = Elsq_data .* 1; %[1000;500];
D0 = diag(Elsq); % prior expectation of extent matrix, increase for MORE variance
al0 = al0_data;%./al0_data*2; %[alpha_IG_D; alpha_IG_D];% object extent covariance matrix parameters, increase for LESS variance
bl0 = Elsq.*(al0 - 1); % object extent covariance matrix parameters

%mean number of measurements stuff
ElambdaM = ElambdaM_data .* 1; % 30; % prior expectation of number of measurements, increase for MORE variance
aM0 = aM0_data; %2; % number of measurements parameters, increase for LESS variance
bM0 = aM0/ElambdaM; % number of measurements parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%motion model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaT = 1;
F = [1, 0, deltaT*1, 0; 0, 1, 0, deltaT*1; 0, 0, 1, 0; 0, 0 , 0, 1]; % state evolution matrix

%ROI
ROI_x = 4000;
ROI_y = 4000;

%initial object states
init_radius = 100;
init_velo = 10;
Ex1 = [0; 0; 0; 0]; % mean for the initial object state (location and velocity)
Sigmax1 = diag([init_radius^2, init_radius^2, init_velo^2, init_velo^2]); % object generation covariance matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%measurement model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%
H = [eye(2) zeros(2,2)];

clutter_density = 5*10^-3; % clutter density
lambdaclutter = 1; % clutter density

clutter_gauss_approx_cov = [1/12*(2*ROI_x)^2, 0; 0, 1/12*(2*ROI_y)^2];


hyperparams = struct("alphadp", alphadp, "Sigmav0", Sigmav0, "av0", av0, "bv0", bv0, "D0", D0, "al0", al0, "bl0", bl0, "ElambdaM", ElambdaM, "aM0", aM0, "bM0", bM0, "alphadp_data", alphadp_data, "Sigmav0_data", Sigmav0_data, "av0_data", av0_data, "bv0_data", bv0_data, "D0_data", D0_data, "al0_data", al0_data, "bl0_data", bl0_data, "ElambdaM_data", ElambdaM_data, "aM0_data", aM0_data, "bM0_data", bM0_data, "F", F, "Ex1", Ex1, "Sigmax1", Sigmax1, "H", H, "clutter_density", clutter_density, 'lambdaclutter', lambdaclutter,"radius", init_radius, "velo", init_velo, "ROI_x", ROI_x, "ROI_y", ROI_y, "clutter_gauss_approx_cov", clutter_gauss_approx_cov, "Essq_data", Essq_data, "Elsq_data", Elsq_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other parameters for method and script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sims = 400;

J = 2; %iterations between mcmc and sequential tracker

factorials = [1 cumprod(1:150) Inf];

statustextoutput = 1;
plotflag = 0;
traj_mode = 2; %trajectory modes: 1=gaussian+random, 2=starts on circle moves towards center, 3=starts on parallel lines
if traj_mode ~= 1
    Sigmax1 = diag([init_radius^2/2, init_radius^2/2, init_velo^2/2, init_velo^2/2]); % object generation covariance matrix
end
clutterflag = 1;
driving_noise_only = 0;

%parameters for OSPA metric
p = 1; cutoff = 120; cutoff_GWD = 1000; cutoff_mse = cutoff*100.*[Elsq_data.^2; Essq_data.^2; ElambdaM_data]; alpha = 2; %OSPA metric parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DA and object existence detection parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxT = 100; % maximum number of observations per object
MTsnz = 1; % set initial number of measurements for each object >= 1
maxnm = 100; % maximum number of measurements per object and timestep

%parameters for gating
gate_parameter = 1-10^-16;
thresh = 10*chi2inv(gate_parameter, 2);

%parameters for gating in DA during initialization
gate_parameter_initial = 1-10^-16;
initial_thresh = 10*chi2inv(gate_parameter_initial, 2);
%initial_thresh = 10*chi2inv(gate_parameter_initial, 2);

%parameters for gating clusters
gate_parameter_cluster = 1-10^-16;
thresh_cluster = 100000*chi2inv(gate_parameter_cluster, 2);

%measurement em cov compensation
cov_comp_init = 0.04;
cov_comp = 0.04;

%birth and death stuff
binoN = 1; % parameter for object birth, number of POs born per time step, cf. I_{birth} in manuscript
Nbirth = 1; % number of POs born per time step, cf. I_{birth} in manuscript
birthT = 1;
binop = 1; % parameter for object birth
ps = 1/10; %= 2/3; survival probability of objects
minT = 80;
stopbirthT = 16;

maxN = binoN*maxT; % number of objects
use_binoN_init = 0;

birth_heuristic_memory = 2; %heuristic parameter, new PO generated when unmatched measurements persist for over x consectutive time steps
death_heuristic_memory = 2; %heuristic parameter, PO track is ended after x consecutive missed detections

fixed_birth_flag = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%object state estimation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sigmarot0 = diag([0.6, 0.005, 0.005]);
SigmaXinit = diag([900, 900, 16, 16]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MCMC parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mult = 1;
L = 100; % number of MCMC iterations, adjust for better results
bis = 20; % number of burn in samples
R = 10; % number of particles per target
Rinit = 10; % number of particles per target during initialization
fCi = round(1*L);%101;%9900; % fix C to final result after fCi iterations, set about 10% below L
Xsi = 1; % save current X every Xsi iterations
Xrsi = 1; % resample X every Xrsi iterations
DPsi = 1; % save current DP data every DPsi iterations
DPrsi = 1; % resample DP data every DPrsi iterations
Dfs = 1; % set D to D0 for first Dfs iterations
Qfs = 1; % set Q to Sigmav0 for first Sigmavfs iterations


mcmc_params = struct("L",L, "bis",bis, "R", R, "Rinit", Rinit, "fCi", fCi, "Xsi", Xsi, "Xrsi", Xrsi, "DPsi",DPsi,"DPrsi",DPrsi,"Dfs",Dfs,"Qfs",Qfs);
script_params = struct("sims", sims, "J", J, "statustextoutput", statustextoutput, "traj_mode", traj_mode, "fixed_dataparams_flag", fixed_dataparams_flag, "fixed_birth_flag", fixed_birth_flag, "maxT", maxT, "Nbirth", Nbirth, "maxN", maxN, "MTsnz", MTsnz, "maxnm", maxnm, "binoN", binoN, "binop", binop, "ps", ps, "birth_heuristic_memory", birth_heuristic_memory, "death_heuristic_memory", death_heuristic_memory, "clutterflag", clutterflag, "initial_thresh", initial_thresh, "thresh", thresh, "thresh_cluster", thresh_cluster, "driving_noise_only", driving_noise_only, "factorials", factorials, "known_params", 1, "known_labels", 1, 'Sigmarot0', Sigmarot0, 'SigmaXinit', SigmaXinit, 'birthT',birthT, 'stopbirthT', stopbirthT, 'cov_comp', cov_comp, "minlamb", minlamb, "maxlamb", maxlamb, "minT", minT, "use_binoN_init", use_binoN_init);
script_params_init = struct("sims", sims, "J", J, "statustextoutput", statustextoutput, "traj_mode", traj_mode, "fixed_dataparams_flag", fixed_dataparams_flag, "fixed_birth_flag", fixed_birth_flag, "maxT", maxT, "Nbirth", Nbirth, "maxN", maxN, "MTsnz", MTsnz, "maxnm", maxnm, "binoN", binoN, "binop", binop, "ps", ps, "birth_heuristic_memory", birth_heuristic_memory, "death_heuristic_memory", death_heuristic_memory, "clutterflag", clutterflag, "initial_thresh", initial_thresh, "thresh", thresh, "thresh_cluster", thresh_cluster, "driving_noise_only", driving_noise_only, "factorials", factorials, "known_params", 1, "known_labels", 1, 'Sigmarot0', Sigmarot0, 'SigmaXinit', SigmaXinit, 'birthT',birthT, 'stopbirthT', stopbirthT, 'cov_comp', cov_comp_init, "use_binoN_init", use_binoN_init);


generatedataflag = 1 && ~loaddata_flag;
only_measurements_flag = 0;


sim = 1; Cmax_DPavg=0; Cmax_EMavg=0; totalcard = [];
while sim <= sims
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate dataset
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if generatedataflag
        if sim == 1 || ~only_measurements_flag
            [data] = generatedata(hyperparams, script_params);
        else
            [data, hyperparams] = generatemeasurements(data, hyperparams, script_params);
        end
    elseif loaddata_flag && ~mod(sim, loadsi) || sim == 1
        file_inc = file_inc + 1;
        load(strcat('data\', matfiles(file_inc).name), 'data');
    elseif loaddata_flag
        [data] =generatemeasurements(data, hyperparams, script_params);
    end
    hyperparams_noclustering = hyperparams;
    hyperparams_noclustering.alphadp = inf;
    hyperparams_fusion = hyperparams;
    al0_fus = 70.*[1; 1]; bl0_fus = Elsq.*(al0_fus - 1); % object extent covariance matrix parameters
    hyperparams_fusion.al0 = al0_fus; hyperparams_fusion.bl0 = bl0_fus; %hyperparams_fusion.av0 = av0_data; hyperparams_fusion.bv0 = bv0_data;

    script_nc = script_params; %script_nc.cov_comp = 1;
    data.q = [squeeze(data.Sigmavstar(1,1,data.C)).' ; squeeze(data.Sigmavstar(3,3,data.C)).'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize tracking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    object_class_est = struct('Dest', repmat(hyperparams.D0, 1, 1, binoN), 'Qest', repmat(hyperparams.Sigmav0, 1, 1, binoN), 'lambdaMest', hyperparams.ElambdaM*ones(binoN,1), 'Cest',(1:binoN).', 'existsteps', zeros(1, binoN), 'Msum', zeros(1, binoN), 'Zsqsum', zeros(2, binoN), 'Vsqsum', zeros(4, binoN));
    init_tracking = sequentialTracking(data, hyperparams, object_class_est, script_params_init, mcmc_params);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% perform tracking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Perfresults = tracking_w_knownclass(data, hyperparams, script_params, mcmc_params, init_tracking, sim);
    [DPresults, ~, ~] = joint_tracking_and_classlearning(data, hyperparams, script_params, mcmc_params, init_tracking, sim);
    [NCresults, ~, ~] = joint_tracking_and_classlearning(data, hyperparams_noclustering, script_nc, mcmc_params, init_tracking, sim);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% metrics and evaluation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Known parameters metrics
    [OSPA_Perf_nolabels(sim), OSPA_Perf_nolabels_pert(:, sim), ~, ~] = calculate_OSPA_notracklabels(data.X, Perfresults.X, p, cutoff, alpha);
    [OSPA_Perf_GWD_nolabels(sim), OSPA_Perf_GWD_nolabels_pert(:, sim), ~, ~, MSE_Perf_nolabels(:,sim), ~, ~]  = calculate_OSPA_notracklabels(data.X, Perfresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), Perfresults.D, data.q, Perfresults.q, data.lambdaMstar(data.C), Perfresults.lambdaM, data.C, Perfresults.object_class_est.Cest);
    [OSPA_Perf(sim), ~, ~, ~] = calculate_OSPA(data.X, Perfresults.X, data.Ts, data.Te, Perfresults.Ts, Perfresults.Te, p, cutoff, alpha);
        
    %% DP clustering metrics
    [OSPA_DP_nolabels(sim), OSPA_DP_nolabels_pert(:, sim), ~, ~] = calculate_OSPA_notracklabels(data.X, DPresults.X, p, cutoff, alpha);
    [OSPA_DP_GWD_nolabels(sim), OSPA_DP_GWD_nolabels_pert(:, sim), ~, ~, MSE_DP_nolabels(:,sim), ~, ~]  = calculate_OSPA_notracklabels(data.X, DPresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), DPresults.D, data.q, DPresults.q, data.lambdaMstar(data.C), DPresults.lambdaM, data.C, DPresults.object_class_est.Cest);
    [OSPA_DP(sim), ~, ~, ~] = calculate_OSPA(data.X, DPresults.X, data.Ts, data.Te, DPresults.Ts, DPresults.Te, p, cutoff, alpha);
    
    %% No clustering metrics
    [OSPA_NC_nolabels(sim), OSPA_NC_nolabels_pert(:, sim), ~, ~] = calculate_OSPA_notracklabels(data.X, NCresults.X, p, cutoff, alpha);
    [OSPA_NC_GWD_nolabels(sim), OSPA_NC_GWD_nolabels_pert(:, sim), ~, ~, MSE_NC_nolabels(:,sim), ~, ~]  = calculate_OSPA_notracklabels(data.X, NCresults.X, p, cutoff_GWD, alpha, cutoff_mse, data.Dstar(:,:, data.C), NCresults.D, data.q, NCresults.q, data.lambdaMstar(data.C), NCresults.lambdaM, data.C, NCresults.object_class_est.Cest);
    [OSPA_NC(sim), ~, ~, ~] = calculate_OSPA(data.X, NCresults.X, data.Ts, data.Te, NCresults.Ts, NCresults.Te, p, cutoff, alpha);
    
    %% Save data
    current_dir = pwd;
    filename_results = strcat(pwd, '/results/', 'results', '.mat');
    save(filename_results, 'OSPA_Perf_GWD_nolabels_pert', 'OSPA_Perf_nolabels', 'OSPA_Perf', 'OSPA_NC_GWD_nolabels_pert', 'OSPA_NC_nolabels', 'OSPA_NC', 'OSPA_DP_GWD_nolabels_pert', 'OSPA_DP_nolabels', 'OSPA_DP', 'maxT')
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% iterate sim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sim = sim + 1
end

plotOSPAGWD_pert