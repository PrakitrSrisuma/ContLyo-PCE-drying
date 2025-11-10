function dXdt = ODE_1Dry(t,X,PAR,Tb_profile,Tc_profile)

 T_1 = X(1);
 T_2 = X(2);
 T_3 = X(3);
 S = X(4);
 
 Tb = cal_Tb(t,Tb_profile);
 Tc = cal_Tc(t,Tc_profile);
 Rp = PAR.Rp0 + PAR.Rp1*S/(1+S*PAR.Rp2);
 PwT = exp(-6139.9./T_1 + 28.8912);
 Pwc = 3;
 sw = switching_1Dry(PAR.H,S);
 Nw = (PwT-Pwc)/Rp;
 
 ddt_T_1 = sw*(2*((PAR.k/PAR.rho*PAR.Cp)/PAR.dz^2)*(1/(PAR.H-S)^2)*(T_2-T_1 - Nw*PAR.dz*PAR.dHsub*(PAR.H-S)/PAR.k- (PAR.F2*PAR.SB*PAR.dz)*(PAR.H-S)*(T_1^4-Tc^4)/PAR.k)-(0-1)*(Nw/(PAR.rho-PAR.rhoe))/(PAR.H-S)*((PAR.H-S)*Nw*PAR.dHsub/PAR.k + PAR.F2*PAR.SB*(T_1^4-Tc^4)*(PAR.H-S)/PAR.k)  - PAR.F1*PAR.SB*PAR.A2*(T_1^4-Tc^4)/(PAR.Ac*(PAR.H-S)*PAR.rho*PAR.Cp));
 ddt_T_2 = sw*(((PAR.k/PAR.rho*PAR.Cp)/PAR.dz^2)*(1/(PAR.H-S)^2)*(T_1-2*T_2+T_3) - ((2-1)*PAR.dz-1)*(Nw/(PAR.rho-PAR.rhoe))/(PAR.H-S)*(T_3-T_1)/(2*PAR.dz)- PAR.F1*PAR.SB*PAR.A2*(T_2^4-Tc^4)/(PAR.Ac*(PAR.H-S)*PAR.rho*PAR.Cp));
 ddt_T_3 = sw*(2*((PAR.k/PAR.rho*PAR.Cp)/PAR.dz^2)*(1/(PAR.H-S)^2)*(T_2-T_3+ (S-PAR.H)*PAR.h*PAR.dz*(T_3-Tb)/PAR.k)-((3-1)*PAR.dz-1)*(Nw/(PAR.rho-PAR.rhoe)/(PAR.H-S))*((S-PAR.H)*PAR.h*(T_3-Tb)/PAR.k) - PAR.F1*PAR.SB*PAR.A2*(T_3^4-Tc^4)/(PAR.Ac*(PAR.H-S)*PAR.rho*PAR.Cp));
 ddt_S = sw*Nw/(PAR.rho-PAR.rhoe);

 dXdt = [ddt_T_1; ddt_T_2; ddt_T_3; ddt_S];
end