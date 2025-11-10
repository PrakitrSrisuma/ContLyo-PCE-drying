function outputsPCE = get_inputPCE_2Dry(ip)

% Extract input data
m = ip.nz3;

% Define all states
states = struct([]);
for i = 1:m
    states(i).name = ['T_' num2str(i)];
    states(i).dist = 'none';
    states(i).data = ip.T03; % lower and upper bound
    states(i).rhs  = ['((k/rho*Cp)/dz^2)*(T_',num2str(i-1),'-2*T_' num2str(i) '+T_' num2str(i+1) ') ' ...
        '+ (rhod*dHdes/(rho*Cp))*c_' num2str(i) '*(-fa*exp(-Ea/(R*T_' num2str(i) ')))'... 
        '- F1*SB*A3*(T_' num2str(i) '^4-Tc^4)/(V*rho*Cp)'];
end
states(1).rhs  = ['2*((k/rho*Cp)/dz^2)*(T_2-T_1) + (rhod*dHdes/(rho*Cp))*c_1*(-fa*exp(-Ea/(R*T_1)))' ...
    '- F1*SB*A3*(T_1^4-Tc^4)/(V*rho*Cp) - 2*(SB*F2/(rho*Cp*dz))*(T_1^4-Tc^4)'];
states(m).rhs  = ['2*((k/rho*Cp)/dz^2)*(T_' num2str(m-1) '-T_' num2str(m) ') ' ...
    '+ (rhod*dHdes/(rho*Cp))*c_' num2str(i) '*(-fa*exp(-Ea/(R*T_' num2str(m) ')))' ...
    '-2*(h/(rho*Cp*dz))*(T_' num2str(m) '-Tb)- F1*SB*A3*(T_' num2str(m) '^4-Tc^4)/(V*rho*Cp)'];


for i = 1:m
    states(i+m).name = ['c_' num2str(i)];
    % states(i+m).dist = 'none';
    % states(i+m).data = ip.cw0;
    states(i+m).dist = 'normal';
    states(i+m).data = [ip.cw0 ip.cw0*0.2];
    states(i+m).rhs  = ['c_' num2str(i) '*(-fa*exp(-Ea/(R*T_' num2str(i) ')))'];
end

% Define parameters and distribution
parameters(1).name = 'dz';
parameters(1).dist = 'none'; 
parameters(1).data = ip.dz3;

parameters(2).name = 'rho';
parameters(2).dist = 'none'; 
parameters(2).data = ip.rhoe; 

parameters(3).name = 'rhod';
parameters(3).dist = 'none'; 
parameters(3).data = ip.rhod; 

parameters(4).name = 'Cp';
parameters(4).dist = 'none'; 
parameters(4).data = ip.Cpe; 

parameters(5).name = 'k';
parameters(5).dist = 'none'; 
parameters(5).data = ip.ke; 

parameters(6).name = 'Ea';
parameters(6).dist = 'none'; 
parameters(6).data = ip.Ea;

parameters(7).name = 'R';
parameters(7).dist = 'none'; 
parameters(7).data = ip.R; 

parameters(8).name = 'fa';
% parameters(8).dist = 'none'; 
% parameters(8).data = ip.fa; 
parameters(8).dist = 'uniform'; 
parameters(8).data = [.3 .5]; % lower and upper bound

parameters(9).name = 'h';
% parameters(9).dist = 'none'; 
% parameters(9).data = ip.hb3;
parameters(9).dist = 'normal'; 
parameters(9).data = [ip.hb3 3]; % mean and SD

parameters(10).name = 'dHdes';
parameters(10).dist = 'none'; 
parameters(10).data = ip.dHdes; 

parameters(11).name = 'F1';
parameters(11).dist = 'none'; 
parameters(11).data = ip.F3*ip.eps1;

parameters(12).name = 'A3';
parameters(12).dist = 'none'; 
parameters(12).data = ip.A3;

parameters(13).name = 'V';
parameters(13).dist = 'none'; 
parameters(13).data = ip.Ac*ip.H3;

parameters(14).name = 'F2';
parameters(14).dist = 'none'; 
parameters(14).data = ip.eps1;

parameters(15).name = 'SB';
parameters(15).dist = 'none'; 
parameters(15).data = ip.SB;

% Define outputs
outputs(1).name = 'c_avg';
outputs(1).rhs = '(c_1+';
for i = 2:m-1
    outputs(1).rhs  = [outputs(1).rhs 'c_' num2str(i) '+'];
end
outputs(1).rhs  = [outputs(1).rhs 'c_' num2str(i) ')*(1/' num2str(m) ')'];

outputs(2).name = 'T_avg';
outputs(2).rhs = '(T_1+';
for i = 2:m-1
    outputs(2).rhs  = [outputs(2).rhs 'T_' num2str(i) '+'];
end
outputs(2).rhs  = [outputs(2).rhs 'T_' num2str(i) ')*(1/' num2str(m) ')'];

% Define inputs
inputs(1).name = 'Tb';
inputs(1).rhs = 'cal_Tb(t,Tb_profile)';
inputs(1).Tb_profile  = ip.Tb3;

inputs(2).name = 'Tc';
inputs(2).rhs = 'cal_Tc(t,Tc_profile)';
inputs(2).Tc_profile  = ip.Tc3;

% Export
outputsPCE.inputs = inputs;
outputsPCE.parameters = parameters;
outputsPCE.outputs = outputs;
outputsPCE.states = states;

return