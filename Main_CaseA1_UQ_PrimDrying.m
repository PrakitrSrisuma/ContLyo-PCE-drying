% ==============================================================================
% This is a top-level routine for Case A1.
% Uncertainty quantification for primary drying in countinuous lyophilization.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group, MIT.
% ==============================================================================
clc, clear, close all
addpath('PoCETfunctions','auxiliary','Calculations','Plotting', ...
    'Exporting Graphics','Input Data','Figures','Outputs');
rng(10)  % make the results deterministic


%% Exract all input data
input = get_inputdata;  % default inputs
ip = input_processing(input);
ipPCE = get_inputPCE_1Dry(ip);
m = ip.nz2;
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;


%% Setup
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


%% Generating Fig. 3 in the ACC paper
% Apply PCE
tic;
col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
results  = PoCETsimCollocation(sys,'ODE_1Dry','OUTPUT_1Dry',col,simoptions);
S_PCE = PoCETsample(sys, 'PCE', mc_samples, results.S.pcvals);  % sample using PCE coefficients
S_PCE(S_PCE>ip.H2) = ip.H2;
Tavg_PCE = PoCETsample(sys, 'PCE', mc_samples, results.T_avg.pcvals);  % sample using PCE coefficients
S_PCEsamp = S_PCE(:,end);
Tavg_PCEsamp = Tavg_PCE(:,end);
[T_PCE, Tmin_PCE, Tmax_PCE] = cal_CI(Tavg_PCE,.95);
[S_PCE, Smin_PCE, Smax_PCE] = cal_CI(S_PCE,.95);
toc;

% Apply MC
mc = PoCETsample(sys,'variables',mc_samples);
tic; mcresults = PoCETsimMonteCarlo(sys,'ODE_1Dry','OUTPUT_1Dry',mc,simoptions,'method','complete'); toc;
S_MCsamp = mcresults.S.mcvals(:,end);
S_MCsamp(S_MCsamp>ip.H2) = ip.H2;
S_MC(2,:) = var(mcresults.S.mcvals,0,1);
S_MC(1,:) = mean(mcresults.S.mcvals,1);
Tavg_MCsamp = mcresults.T_avg.mcvals(:,end);
Tavg_MC(2,:) = var(mcresults.T_avg.mcvals,0,1);
Tavg_MC(1,:) = mean(mcresults.T_avg.mcvals,1);

% Plotting
Nhis = 15;
PrimDrying_PCE_MC = figure;
tiledlayout(1,4,"TileSpacing","compact","Padding","compact")
linePCE = {'Color', [1 0 0],'Linewidth',2};
lineMC = {'Linestyle',':','Color', [0 0 1],'Linewidth',2.5};
Smin = (floor(min([S_MCsamp;S_PCEsamp])*1e4)/1e4)*100;
Smax = (ceil(max([S_MCsamp;S_PCEsamp])*1e4)/1e4)*100;
Tmin = floor(min([Tavg_MCsamp;Tavg_PCEsamp]));
Tmax = ceil(max([Tavg_MCsamp;Tavg_PCEsamp]));

nexttile(1)
plot_T_MeanCI(results.time/3600,T_PCE,Tmin_PCE,Tmax_PCE,{'Color', [0 100/256 0],'Linewidth',1.5},{[0 1 .5], 'FaceAlpha',.3,'LineStyle','none'})
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by4')
text(-.33,1,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
xlim([0 7])
xticks(0:1:7)

nexttile(2)
dT = 0.3*(Tmax-Tmin)/Nhis;  % shift the plot to make both diagrams visible
plot_T_distrb(Tavg_PCEsamp,linspace(Tmin-dT,Tmax,Nhis),'1a'); hold on
plot_T_distrb(Tavg_MCsamp,linspace(Tmin,Tmax+dT,Nhis),'1b');
lg = legend({'PCE','MC'},'location','northwest'); lg.ItemTokenSize(1) = 15;
xlim([Tmin-dT Tmax+dT])
xticks(linspace(Tmin-dT,Tmax+dT,4))
graphics_setup('1by4')
text(-.33,1,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');

nexttile(3)
plot_S_MeanCI(results.time/3600,S_PCE*100,Smin_PCE*100,Smax_PCE*100,{'Color', [0 100/256 0],'Linewidth',1.5},{[0 1 .5], 'FaceAlpha',.3,'LineStyle','none'})
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by4')
xlim([0 7])
xticks(0:1:7)
text(-.33,1,'(C)','Units','normalized','FontSize', 10,'fontweight', 'bold');

nexttile(4)
dS = .3*(Smax-Smin)/Nhis;  % shift the plot to make both diagrams visible
plot_S_distrb(S_PCEsamp*100,linspace(Smin-dS,Smax,Nhis),'1a'); hold on
plot_S_distrb(S_MCsamp*100,linspace(Smin,Smax+dS,Nhis),'1b');
xlim([Smin-dS Smax+dS])
xticks(linspace(Smin-dS,Smax+dS,4))
lg = legend({'PCE','MC'},'location','west'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by4')
text(-.33,1,'(D)','Units','normalized','FontSize', 10,'fontweight', 'bold');

% Save plot (uncomment the below code)
% export_figures(PrimDrying_PCE_MC,'PrimDrying_PCE_MC');
