% ==============================================================================
% This is a top-level routine for Case A2.
% Uncertainty quantification for secondary drying in countinuous lyophilization.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group, MIT.
% ==============================================================================
clc, clear, close all
addpath('PoCETfunctions','auxiliary','Calculations','Plotting', ...
    'Exporting Graphics','Input Data','Figures','Outputs');
rng(127)  % make the results deterministic


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
PCE_Fig4 = 'on';
PCE_Fig5 = 'off';

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


%%  PCE + MC (Figure used in the ACC paper)
switch PCE_Fig4 
case 'on'

% Apply PCE
tic;
col  = PoCETsample(sys,'basis',col_samples);  % get collocation points
results  = PoCETsimCollocation(sys,'ODE_2Dry','OUTPUT_2Dry',col,simoptions);
cavg_PCE = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
Tavg_PCE = PoCETsample(sys, 'PCE', mc_samples, results.T_avg.pcvals);  % sample using PCE coefficients
cavg_PCEsamp = cavg_PCE(:,end);
Tavg_PCEsamp = Tavg_PCE(:,end);
[T_PCE, Tmin_PCE, Tmax_PCE] = cal_CI(Tavg_PCE,.95);
[c_PCE, cmin_PCE, cmax_PCE] = cal_CI(cavg_PCE,.95);
toc;

% Apply MC
mc = PoCETsample(sys,'variables',mc_samples);
tic; mcresults = PoCETsimMonteCarlo(sys,'ODE_2Dry','OUTPUT_2Dry',mc,simoptions,'method','complete'); toc;
cavg_MCsamp = mcresults.c_avg.mcvals(:,end);
cavg_MC(2,:) = var(mcresults.c_avg.mcvals,0,1);
cavg_MC(1,:) = mean(mcresults.c_avg.mcvals,1);
Tavg_MCsamp = mcresults.T_avg.mcvals(:,end);
Tavg_MC(2,:) = var(mcresults.T_avg.mcvals,0,1);
Tavg_MC(1,:) = mean(mcresults.T_avg.mcvals,1);

% Plotting
Nhis = 15;
SecDrying_PCE_MC = figure;
tiledlayout(1,4,"TileSpacing","compact","Padding","compact")
linePCE = {'Color', [1 0 0],'Linewidth',2};
lineMC = {'Linestyle',':','Color', [0 0 1],'Linewidth',2.5};
cmin = floor(min([cavg_MCsamp;cavg_PCEsamp])*1e3)/1e3;
cmax = ceil(max([cavg_MCsamp;cavg_PCEsamp])*1e3)/1e3;
Tmin = floor(min([Tavg_MCsamp;Tavg_PCEsamp])*1e2)/1e2;
Tmax = ceil(max([Tavg_MCsamp;Tavg_PCEsamp])*1e2)/1e2;

nexttile(1)
plot_T_MeanCI(results.time/3600,T_PCE,Tmin_PCE,Tmax_PCE,{'Color', [0 100/256 0],'Linewidth',1.5},{[0 1 .5], 'FaceAlpha',.3,'LineStyle','none'}); hold on
text(-.33,1,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
xlim([0 7])
xticks(0:1:7)
graphics_setup('1by4')
% ylim([270 300])

nexttile(2)
dT = 0.3*(Tmax-Tmin)/Nhis;
plot_T_distrb(Tavg_PCEsamp,linspace(Tmin-dT,Tmax,Nhis),'1a'); hold on
plot_T_distrb(Tavg_MCsamp,linspace(Tmin,Tmax+dT,Nhis),'1b');
text(-.33,1,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');
lg = legend({'PCE','MC'},'location','west'); lg.ItemTokenSize(1) = 15;
xlim([Tmin-dT Tmax+dT])
xticks(linspace(Tmin-dT,Tmax+dT,4))
graphics_setup('1by4')

nexttile(3)
plot_cw_MeanCI(results.time/3600,c_PCE,cmin_PCE,cmax_PCE,{'Color', [0 100/256 0],'Linewidth',1.5},{[0 1 .5], 'FaceAlpha',.3,'LineStyle','none'})
text(-.33,1,'(C)','Units','normalized','FontSize', 10,'fontweight', 'bold');
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
xlim([0 7])
xticks(0:1:7)
ylim([0 .12])
yticks(0:.02:.12)
graphics_setup('1by4')

nexttile(4)
dc = 0.3*(cmax-cmin)/Nhis;
plot_cw_distrb(cavg_PCEsamp,linspace(cmin-dc,cmax,Nhis),'1a'); hold on
plot_cw_distrb(cavg_MCsamp,linspace(cmin,cmax+dc,Nhis),'1b');
xlim([cmin-dc cmax+dc])
xticks(linspace(cmin-dc,cmax+dc,4))
text(-.33,1,'(D)','Units','normalized','FontSize', 10,'fontweight', 'bold');
lg = legend({'PCE','MC'},'location','best'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by4')

% Save plot (uncomment the below code)
% export_figures(SecDrying_PCE_MC,'SecDrying_PCE_MC'); 

end


%%  PCE for 3 cases (Figure used in the ACC paper)
switch PCE_Fig5
case 'on'

SecDrying_PCE_All = figure;
tiledlayout(1,3,"TileSpacing","loose","Padding","compact")

% Concentration only
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;
parameters(8).dist = 'none'; 
parameters(8).data = ip.fa;
parameters(9).dist = 'none'; 
parameters(9).data = ip.hb3;
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
col  = PoCETsample(sys,'basis',col_samples); 

tic;
results  = PoCETsimCollocation(sys,'ODE_2Dry','OUTPUT_2Dry',col,simoptions);
cavg_PCE = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
[c_PCE, cmin_PCE, cmax_PCE] = cal_CI(cavg_PCE,.95);
toc;

nexttile(1); plot_cw_MeanCI(results.time/3600,c_PCE,cmin_PCE,cmax_PCE,{'Color', [0 0 1],'Linewidth',2},{[0 0.5 1], 'FaceAlpha',.3,'LineStyle','none'})
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by3')
text(.02,.95,'(A)','Units','normalized','FontSize', 10,'fontweight', 'bold');
xlim([0 7])
xticks(0:1:7)
ylim([0 .12])
yticks(0:.02:.12)

% Kinetics only
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;
for i = 1:m
    states(i+m).dist = 'none';
    states(i+m).data = ip.cw0;
end
parameters(9).dist = 'none'; 
parameters(9).data = ip.hb3;
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
col  = PoCETsample(sys,'basis',col_samples); 
tic;
results  = PoCETsimCollocation(sys,'ODE_2Dry','OUTPUT_2Dry',col,simoptions);
cavg_PCE = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
[c_PCE, cmin_PCE, cmax_PCE] = cal_CI(cavg_PCE,.95);
toc;

nexttile(2)
plot_cw_MeanCI(results.time/3600,c_PCE,cmin_PCE,cmax_PCE,{'Color', [1 0 .0],'Linewidth',2},{[1 0.2 .2], 'FaceAlpha',.3,'LineStyle','none'})
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by3')
text(.02,.95,'(B)','Units','normalized','FontSize', 10,'fontweight', 'bold');
xlim([0 7])
xticks(0:1:7)


% HTC only
states = ipPCE.states;
parameters = ipPCE.parameters;
outputs = ipPCE.outputs;
inputs = ipPCE.inputs;
for i = 1:m
    states(i+m).dist = 'none';
    states(i+m).data = ip.cw0;
end
parameters(8).dist = 'none'; 
parameters(8).data = ip.fa;
sys = PoCETcompose(states,parameters,inputs,outputs,pce_order);
col  = PoCETsample(sys,'basis',col_samples); 
tic;
results  = PoCETsimCollocation(sys,'ODE_2Dry','OUTPUT_2Dry',col,simoptions);
cavg_PCE = PoCETsample(sys, 'PCE', mc_samples, results.c_avg.pcvals);  % sample using PCE coefficients
[c_PCE, cmin_PCE, cmax_PCE] = cal_CI(cavg_PCE,.95);
toc;

nexttile(3)
plot_cw_MeanCI(results.time/3600,c_PCE,cmin_PCE,cmax_PCE,{'Color', [.5 0 .7],'Linewidth',2},{[.7 0 .5], 'FaceAlpha',.3,'LineStyle','none'})
lg = legend({'Mean','95% CI'},'location','best'); lg.ItemTokenSize(1) = 15;
graphics_setup('1by3')
text(.02,.95,'(C)','Units','normalized','FontSize', 10,'fontweight', 'bold');
xlim([0 7])
xticks(0:1:7)

% Save plot (uncomment the below code)
% export_figures(SecDrying_PCE_All,'SecDrying_PCE_All');

end