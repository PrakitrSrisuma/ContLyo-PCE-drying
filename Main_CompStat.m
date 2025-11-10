% ==============================================================================
% This is a top-level routine for PCE modeling.
% Measuring computational performance.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
clc, clear, close all
addpath('PoCETfunctions','auxiliary','Calculations','Plotting', ...
    'Exporting Graphics','Input Data','Figures','Outputs','Computational Performance');

%% Setup
% Define mode of calculation
PCE_A1 = 'on';
MC_A1 = 'on';
PCE_A2 = 'on';
MC_A2 = 'on';
PCE_B1 = 'on';
MC_B1 = 'on';
PCE_B2 = 'on';
MC_B2 = 'on';
savedata = 1;


%% PCE for Case A1
switch PCE_A1
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_1Dry(ip);
m = ip.nz2;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.tspan = [0 7*3600];
simoptions.dt = ip.dt2;
simoptions.setup = odeset('RelTol',1e-9,'AbsTol',1e-9);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_1Dry.m');
writeMCOUTfile(sys,'OUTPUT_1Dry.m');

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);

% Apply PCE
for i = 0:nrun
    rng(127)
    tic;
    col  = PoCETsample(sys,'basis',col_samples);  % prepare sampling
    results  = PoCETsimCollocation(sys,'ODE_1Dry','OUTPUT_1Dry',col,simoptions);  % run PCE
    MomMats = PoCETmomentMatrices(sys,4);  % prepare moments
    S_PCE = PoCETcalcMoments(sys,MomMats,results.S.pcvals);  % cal moments
    S_PCE = filterS(S_PCE,ip);
    Savg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.S.pcvals);  % sample using PCE coefficients
    Tavg_PCE = PoCETcalcMoments(sys,MomMats,results.T_avg.pcvals);  % cal moments
    Tavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.T_avg.pcvals);  % sample using PCE coefficients
    
    if i > 0  % neglect the first run
        simtime(i) = toc;
        check(i) = Tavg_PCE(1,end);
    end

end

% Save data
if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_PCE_A1.mat');
    filename2 = fullfile('Computational Performance','check_PCE_A1.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_PCE_A1 = mean(simtime);
tSD_PCE_A1 = std(simtime);

end


%% MC for Case A1
switch MC_A1
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_1Dry(ip);
m = ip.nz2;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.tspan = [0 7*3600];
simoptions.dt = ip.dt2;
simoptions.setup = odeset('RelTol',1e-9,'AbsTol',1e-9);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_1Dry.m');
writeMCOUTfile(sys,'OUTPUT_1Dry.m');

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);

% Apply MC
for i = 0:nrun
    rng(127)
    tic;
    mc = PoCETsample(sys,'variables',mc_samples);
    mcresults = PoCETsimMonteCarlo(sys,'ODE_1Dry','OUTPUT_1Dry',mc,simoptions,'method','complete'); 
    
    if i > 0  % neglect the first run
        simtime(i) = toc;
        check(i) = max(mcresults.T_avg.mcvals,[],'all');
    end

end

% Save data
if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_MC_A1.mat');
    filename2 = fullfile('Computational Performance','check_MC_A1.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_MC_A1 = mean(simtime);
tSD_MC_A1 = std(simtime);

end


%% PCE for Case A2
switch PCE_A2
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.tspan = [0 7*3600];
simoptions.dt = ip.dt3;
simoptions.setup = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_2Dry.m');
writeMCOUTfile(sys,'OUTPUT_2Dry.m');

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);

% Apply PCE
for i = 0:nrun
    rng(127)
    tic;
    col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
    results  = PoCETsimCollocation(sys,'ODE_2Dry','OUTPUT_2Dry',col,simoptions);
    MomMats = PoCETmomentMatrices(sys,4);  % prepare moments
    cavg_PCE = PoCETcalcMoments(sys,MomMats,results.c_avg.pcvals);  % cal moments
    cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
    Tavg_PCE = PoCETcalcMoments(sys,MomMats,results.T_avg.pcvals);  % cal moments
    Tavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.T_avg.pcvals);  % sample using PCE coefficients
    
    if i > 0  % neglect the first run
        simtime(i) = toc;
        check(i) = Tavg_PCE(1,end);
    end

end

if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_PCE_A2.mat');
    filename2 = fullfile('Computational Performance','check_PCE_A2.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_PCE_A2 = mean(simtime);
tSD_PCE_A2 = std(simtime);

end


%% MC for Case A2
switch MC_A2
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.tspan = [0 7*3600];
simoptions.dt = ip.dt3;
simoptions.setup = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_2Dry.m');
writeMCOUTfile(sys,'OUTPUT_2Dry.m')

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);

% Apply MC
for i = 0:nrun
    rng(127)
    tic;
    mc = PoCETsample(sys,'variables',mc_samples);
    mcresults = PoCETsimMonteCarlo(sys,'ODE_2Dry','OUTPUT_2Dry',mc,simoptions,'method','complete'); 
    
    if i > 0  % neglect the first run
        simtime(i) = toc;
        check(i) = max(mcresults.T_avg.mcvals,[],'all');
    end

end

if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_MC_A2.mat');
    filename2 = fullfile('Computational Performance','check_MC2_A2.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_MC_A2 = mean(simtime);
tSD_MC_A2 = std(simtime);

end


%% PCE for Case B1
switch PCE_B1
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.dt = ip.dt3;
simoptions.setup = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_2Dry_opt.m');
writeMCOUTfile(sys,'OUTPUT_2Dry_opt.m');

% Probability
P0 = .95;
Ctarget = 0.01;

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);
col  = PoCETsample(sys,'basis',col_samples);  % get collocation points


% Apply PCE
for k = 0:nrun
    rng(127)
    tic;
    simoptions.tspan = [0 7*3600];
    Tb = linspace(270,330,15);
    nsim = length(Tb);
    P = zeros(nsim,1);
    
    tic;
    for i = 1:nsim
        col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
        results  = PoCETsimCollocation(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',col,simoptions,'Tb_profile', Tb(i));
        cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals(:,end));  % sample using PCE coefficients
    end

    if k > 0  % neglect the first run
        simtime(k) = toc;
        check(k) = cavg_PCEsamp(end);
    end
end

if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_PCE_B1.mat');
    filename2 = fullfile('Computational Performance','check_PCE_B1.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_PCE_B1 = mean(simtime);
tSD_PCE_B1 = std(simtime);

end


%% MC for Case B1
switch MC_B1
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.dt = ip.dt3;
simoptions.setup = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_2Dry_opt.m');
writeMCOUTfile(sys,'OUTPUT_2Dry_opt.m');

% Probability
P0 = .95;
Ctarget = 0.01;

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);
col  = PoCETsample(sys,'basis',col_samples);  % get collocation points


% MC Simulation
for k = 0:nrun
    rng(127)
    tic;
    simoptions.tspan = [0 7*3600];
    Tb = linspace(270,330,15);
    nsim = length(Tb);
    P = zeros(nsim,1);
    
    tic;
    for i = 1:nsim
        mc = PoCETsample(sys,'variables',mc_samples);
        mcresults = PoCETsimMonteCarlo(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',mc,simoptions,'method','complete','Tb_profile', Tb(i));
        cavg_MCsamp = mcresults.c_avg.mcvals(:,end);
    end
    if k > 0  % neglect the first run
        simtime(k) = toc;
        check(k) = cavg_MCsamp(end);
    end
end

if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_MC_B1.mat');
    filename2 = fullfile('Computational Performance','check_MC_B1.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_MC_B1 = mean(simtime);
tSD_MC_B1 = std(simtime);

end


%% PCE for Case B2
switch PCE_B2
case 'on'

% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.dt = ip.dt3;
simoptions.setup = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_2Dry_opt.m');
writeMCOUTfile(sys,'OUTPUT_2Dry_opt.m');

% Probability
P0 = .95;
Ctarget = 0.01;

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);

% Run fmincon
tf = 10*3600;
simoptions.tspan = [0 tf];
tspan = (0:ip.dt3:tf)';
col  = PoCETsample(sys,'basis',col_samples); 

% Run fmincon
options = optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','Central','FiniteDifferenceStepSize',1e-3,'OptimalityTolerance',1e-4);
for i = 0:nrun
    rng(127)
    tic;[u,fval] = fmincon(@(u) obj_OCP(u,sys,col,simoptions,P0,Ctarget), 290,[],[],[],[],273,295,[],options); 
    
    if i > 0  % neglect the first run
        simtime(i) = toc;
        check(i) = u;
    end

end

% Save data
if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_PCE_B2.mat');
    filename2 = fullfile('Computational Performance','check_PCE_B2.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_PCE_B2 = mean(simtime);
tSD_PCE_B2 = std(simtime);

end


%% MC for Case B2
switch MC_B2
case 'on'
    
% Setup
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;

% ODE solvers
simoptions.dt = ip.dt3;
simoptions.setup = odeset('RelTol',ip.tol,'AbsTol',ip.tol);
simoptions.solver = 'ode15s';

% Number of sampling points
mc_samples = 2000;
col_samples = 50;
pce_order = 2;  % if changed; create a new system (newprob = 1)
N_samples = factorial(2+pce_order)/(factorial(2)*factorial(pce_order));  % check number of samples

% Create a problem structure
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
writeMCRHSfile(sys,'ODE_2Dry_opt.m');
writeMCOUTfile(sys,'OUTPUT_2Dry_opt.m');

% Probability
P0 = .95;
Ctarget = 0.01;

% Start
nrun = 10;
simtime = zeros(nrun,1);
check = zeros(nrun,1);

% Run fmincon
tf = 10*3600;
simoptions.tspan = [0 tf];
tspan = (0:ip.dt3:tf)';
col  = PoCETsample(sys,'basis',col_samples); 

% Run fmincon
options = optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','Central','FiniteDifferenceStepSize',1e-3,'OptimalityTolerance',1e-4);
for i = 0:nrun
    rng(127)
    tic;[u,fval] = fmincon(@(u) obj_OCP_MC(u,sys,simoptions,P0,Ctarget), 290,[],[],[],[],273,295,[],options); 
    
    if i > 0  % neglect the first run
        simtime(i) = toc;
        check(i) = u;
    end

end

% Save data
if savedata == 1
    filename1 = fullfile('Computational Performance','simtime_MC_B2.mat');
    filename2 = fullfile('Computational Performance','check_MC_B2.mat');
    save(filename1,'simtime'); save(filename2,'check')
end

tmean_MC_B2 = mean(simtime);
tSD_MC_B2 = std(simtime);

end


%% Objective functions
function outputs = obj_OCP(u,sys,col,simoptions,P0,Ctarget)
    
    mc_samples = 2000;
    simoptions.tspan = [0 10*3600];
    results  = PoCETsimCollocation(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',col,simoptions,'Tb_profile', u);
    cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
    P = zeros(length(results.time),1);
    for i = 1:length(results.time)
        c_tmp = find(cavg_PCEsamp(:,i) <= Ctarget);
        P(i) = length(c_tmp)/mc_samples;
    end

    i_opt = find(P>P0,1);
    if isempty(i_opt)
        outputs = 1e7*(1-P(end));
    else
        outputs = results.time(i_opt);
    end

end

function outputs = obj_OCP_MC(u,sys,simoptions,P0,Ctarget)
    
    mc_samples = 2000;
    simoptions.tspan = [0 10*3600];
    mc = PoCETsample(sys,'variables',mc_samples);
    mcresults = PoCETsimMonteCarlo(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',mc,simoptions,'method','complete','Tb_profile', u);
    cavg_MCsamp = mcresults.c_avg.mcvals;
    P = zeros(length(mcresults.time),1);
    for i = 1:length(mcresults.time)
        c_tmp = find(cavg_MCsamp(:,i) <= Ctarget);
        P(i) = length(c_tmp)/mc_samples;
    end

    i_opt = find(P>P0,1);
    if isempty(i_opt)
        outputs = 1e7*(1-P(end));
    else
        outputs = mcresults.time(i_opt);
    end

end