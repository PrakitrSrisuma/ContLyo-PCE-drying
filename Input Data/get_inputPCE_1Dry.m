function outputsPCE = get_inputPCE_1Dry(ip)

% Extract input data
m = ip.nz2;

% Define all states
states = struct([]);
for i = 1:m
    states(i).name = ['T_' num2str(i)];
    states(i).dist = 'none';
    states(i).data = ip.T02; % lower and upper bound
    states(i).rhs  = ['((k/rho*Cp)/dz^2)*(1/(H-S)^2)*(T_',num2str(i-1),'-2*T_' num2str(i) '+T_' num2str(i+1) ') ' ...
        '- ((' num2str(i) '-1)*dz-1)*(Nw/(rho-rhoe))/(H-S)*(T_',num2str(i+1) '-T_' num2str(i-1) ')/(2*dz)' ...
        '- F1*SB*A2*(T_' num2str(i) '^4-Tc^4)/(Ac*(H-S)*rho*Cp)'];
    states(i).rhs = ['sw*(' states(i).rhs ')'];
end
states(1).rhs  = ['2*((k/rho*Cp)/dz^2)*(1/(H-S)^2)*(T_2-T_1 - Nw*dz*dHsub*(H-S)/k- (F2*SB*dz)*(H-S)*(T_1^4-Tc^4)/k)' ...
    '-(0-1)*(Nw/(rho-rhoe))/(H-S)*((H-S)*Nw*dHsub/k + F2*SB*(T_1^4-Tc^4)*(H-S)/k)  - F1*SB*A2*(T_1^4-Tc^4)/(Ac*(H-S)*rho*Cp)'];
states(m).rhs  = ['2*((k/rho*Cp)/dz^2)*(1/(H-S)^2)*(T_' num2str(m-1) '-T_' num2str(m) '+ (S-H)*h*dz*(T_' num2str(m) '-Tb)/k)' ...
    '-((' num2str(m) '-1)*dz-1)*(Nw/(rho-rhoe)/(H-S))*((S-H)*h*(T_' num2str(m) '-Tb)/k) ' ...
    '- F1*SB*A2*(T_' num2str(m) '^4-Tc^4)/(Ac*(H-S)*rho*Cp)'];
states(1).rhs = ['sw*(' states(1).rhs ')'];
states(m).rhs = ['sw*(' states(m).rhs ')'];

states(m+1).name = 'S';
states(m+1).dist = 'none';
states(m+1).data = ip.S0; % lower and upper bound
states(m+1).rhs  = 'sw*Nw/(rho-rhoe)';

% Define parameters and distribution
parameters(1).name = 'dz';
parameters(1).dist = 'none'; 
parameters(1).data = ip.dpsi;

parameters(2).name = 'rho';
parameters(2).dist = 'none'; 
parameters(2).data = ip.rhof; 

parameters(3).name = 'rhoe';
parameters(3).dist = 'none'; 
parameters(3).data = ip.rhoe; 

parameters(4).name = 'Cp';
parameters(4).dist = 'none'; 
parameters(4).data = ip.Cpf; 

parameters(5).name = 'k';
parameters(5).dist = 'none'; 
parameters(5).data = ip.kf; 

parameters(6).name = 'H';
parameters(6).dist = 'none'; 
parameters(6).data = ip.H2;

parameters(7).name = 'h';
% parameters(7).dist = 'none'; 
% parameters(7).data = ip.hb2;
parameters(7).dist = 'normal'; 
parameters(7).data = [ip.hb2 3]; % mean and SD

parameters(8).name = 'dHsub';
parameters(8).dist = 'none'; 
parameters(8).data = ip.dHsub; 

parameters(9).name = 'SB';
parameters(9).dist = 'none'; 
parameters(9).data = ip.SB; 

parameters(10).name = 'F1';
parameters(10).dist = 'none'; 
parameters(10).data = ip.F2*ip.eps1;

parameters(11).name = 'Ac';
parameters(11).dist = 'none'; 
parameters(11).data = ip.Ac;

parameters(12).name = 'A2';
parameters(12).dist = 'none'; 
parameters(12).data = ip.A2;

parameters(13).name = 'F2';
parameters(13).dist = 'none'; 
parameters(13).data = ip.eps1;

parameters(14).name = 'Rp0';
% parameters(14).dist = 'none'; 
% parameters(14).data = ip.Rp0;
parameters(14).dist = 'uniform'; 
parameters(14).data = [1e4 2e4];

parameters(15).name = 'Rp1';
% parameters(15).dist = 'none'; 
% parameters(15).data = ip.Rp1;
parameters(15).dist = 'uniform'; 
parameters(15).data = [1e7 3e7];

parameters(16).name = 'Rp2';
parameters(16).dist = 'none'; 
parameters(16).data = ip.Rp2;

% Define outputs
outputs(1).name = 'T_avg';
outputs(1).rhs = '(T_1+';
for i = 2:m-1
    outputs(1).rhs  = [outputs(1).rhs 'T_' num2str(i) '+'];
end
outputs(1).rhs  = [outputs(1).rhs 'T_' num2str(m) ')*(1/' num2str(m) ')'];

% Define inputs
inputs(1).name = 'Tb';
inputs(1).rhs = 'cal_Tb(t,Tb_profile)';
inputs(1).Tb_profile  = ip.Tb2;

inputs(2).name = 'Tc';
inputs(2).rhs = 'cal_Tc(t,Tc_profile)';
inputs(2).Tc_profile  = ip.Tc2;

inputs(3).name = 'Rp';
inputs(3).rhs  = 'Rp0 + Rp1*S/(1+S*Rp2)';

inputs(4).name = 'PwT';
inputs(4).rhs  = 'exp(-6139.9./T_1 + 28.8912)';

inputs(5).name = 'Pwc';
inputs(5).rhs  = num2str(ip.Pwc);

inputs(6).name = 'sw';
inputs(6).rhs = 'switching_1Dry(H,S)';

inputs(7).name = 'Nw';
inputs(7).rhs  = '(PwT-Pwc)/Rp';

% Export
outputsPCE.inputs = inputs;
outputsPCE.parameters = parameters;
outputsPCE.outputs = outputs;
outputsPCE.states = states;

return