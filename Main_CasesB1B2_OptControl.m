% ==============================================================================
% This is a top-level routine for Cases B1 and B2.
% Stochastic optmization and control for countinuous lyophilization.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group, MIT.
% ==============================================================================
clc, clear, close all
addpath('PoCETfunctions','auxiliary','Calculations','Plotting', ...
    'Exporting Graphics','Input Data','Figures','Outputs');
rng(127)


%% Exract all input data
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_2Dry(ip);
m = ip.nz3;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;


%% Setup
% Choose the figure(s)
OptPCE_Fig6 = 'off';
OptPCE_Fig7 = 'on';

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


%% Multiple simulations to find the shelf temperature
switch OptPCE_Fig6
case 'on'

Nhis = 15;
simoptions.tspan = [0 7*3600];
Tb = linspace(270,330,12);
nsim = length(Tb);
P = zeros(nsim,1);

% Apply PCE
tic;
for i = 1:nsim
    col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
    results  = PoCETsimCollocation(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',col,simoptions,'Tb_profile', Tb(i));
    cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals(:,end));  % sample using PCE coefficients
    
    c_tmp = find(cavg_PCEsamp <= Ctarget);
    P(i) = length(c_tmp)/mc_samples;

end
toc;

% Plotting
SecDrying_Opt2A = figure;
tiledlayout(1,2,"TileSpacing","loose","Padding","compact")
nexttile
plot(Tb,P,'-^b','markerfacecolor','b','linewidth',1,'markersize',3); xlabel('Shelf temperature (K)'); ylabel('Probability');
text(-.25,1,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
xlim([270 330])
xticks(270:20:330)
graphics_setup('1by2')

nexttile
col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
Tb_opt = 310;
results  = PoCETsimCollocation(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',col,simoptions,'Tb_profile', Tb_opt);
MomMats = PoCETmomentMatrices(sys,4);  % prepare moments
cavg_PCE = PoCETcalcMoments(sys,MomMats,results.c_avg.pcvals(:,end));  % cal moments
cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals(:,end));  % sample using PCE coefficients
c_tmp = find(cavg_PCEsamp <= Ctarget);
P_opt = length(c_tmp)/mc_samples;

cavg_opt = cavg_PCEsamp; 
cmin = min(cavg_opt); cmax = max(cavg_opt);
plot_cw_cdf(cavg_opt,linspace(cmin,cmax,Nhis),'1b',P0,Ctarget); hold on
xlim([cmin cmax])
xticks(linspace(cmin,cmax,4))
graphics_setup('1by2')
text(-.27,1,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');

% Save plot (uncomment the below code)
% export_figures(SecDrying_Opt2A,'SecDrying_Opt2A'); 

end


%% Solve a stochastic optimal control problem with PCE
switch OptPCE_Fig7
case 'on'

tf = 10*3600;
simoptions.tspan = [0 tf];
tspan = (0:ip.dt3:tf)';
col  = PoCETsample(sys,'basis',col_samples); 

% Run fmincon
options = optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','Central','FiniteDifferenceStepSize',1e-3,'OptimalityTolerance',1e-4);
tic;[u,fval] = fmincon(@(u) obj_OCP(u,sys,col,simoptions,P0,Ctarget), 290,[],[],[],[],273,295,[],options); toc;
T_opt = u; t_opt = fval;

% Plotting
col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
Tb = linspace(273,295,10);
P = zeros(length(Tb),length(tspan));
for i = 1:length(Tb)
    results  = PoCETsimCollocation(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',col,simoptions,'Tb_profile',Tb(i));
    cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
    for j = 1:length(results.time)
        c_tmp = find(cavg_PCEsamp(:,j) <= Ctarget);
        P(i,j) = length(c_tmp)/mc_samples;
    end
end

SecDrying_Opt2B = figure;
tiledlayout(1,2,"TileSpacing","loose","Padding","compact")

[X,Y] = meshgrid(Tb,tspan);
nexttile
surf(X,Y/3600,P','FaceAlpha',0.85,'FaceColor','interp','EdgeColor','none')
zlabel('Probability');
xlh = xlabel('Shelf temperature (K)'); xlh.Position(2) = xlh.Position(2) - 2;
ylh = ylabel('Time (h)'); ylh.Position(1) = ylh.Position(1) + 5; ylh.Position(2) = ylh.Position(2) - 1;
xticks(linspace(273,295,3)); xlim([273 295]); 
text(-.27,1,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
graphics_setup('1by2')

nexttile
simoptions.tspan = [0 t_opt];
results  = PoCETsimCollocation(sys,'ODE_2Dry_opt','OUTPUT_2Dry_opt',col,simoptions,'Tb_profile',T_opt);
cavg_PCEsamp = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals(:,end)); 
cavg_opt = cavg_PCEsamp; 
cmin = min(cavg_opt); cmax = max(cavg_opt);
plot_cw_cdf(cavg_opt,linspace(cmin,cmax,15),'1b',P0,Ctarget); hold on
xlim([cmin cmax])
xticks(linspace(cmin,cmax,4))
graphics_setup('1by2')
text(-.27,1,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');

% Save plot (uncomment the below code)
export_figures(SecDrying_Opt2B,'SecDrying_Opt2B');

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