function dXdt = ODE_2Dry_opt(t,X,PAR,Tb_profile,Tc_profile)

 T_1 = X(1);
 T_2 = X(2);
 T_3 = X(3);
 c_1 = X(4);
 c_2 = X(5);
 c_3 = X(6);
 
 Tb = cal_Tb(t,Tb_profile);
 Tc = cal_Tc(t,Tc_profile);
 
 ddt_T_1 = 2*((PAR.k/PAR.rho*PAR.Cp)/PAR.dz^2)*(T_2-T_1) + (PAR.rhod*PAR.dHdes/(PAR.rho*PAR.Cp))*c_1*(-PAR.fa*exp(-PAR.Ea/(PAR.R*T_1)))- PAR.F*PAR.hrad*(T_1^4-Tc^4)/(PAR.V*PAR.rho*PAR.Cp) - 2*PAR.ftop*(PAR.hrad_top/(PAR.rho*PAR.Cp*PAR.dz))*(T_1^4-Tc^4);
 ddt_T_2 = ((PAR.k/PAR.rho*PAR.Cp)/PAR.dz^2)*(T_1-2*T_2+T_3) + (PAR.rhod*PAR.dHdes/(PAR.rho*PAR.Cp))*c_2*(-PAR.fa*exp(-PAR.Ea/(PAR.R*T_2)))- PAR.F*PAR.hrad*(T_2^4-Tc^4)/(PAR.V*PAR.rho*PAR.Cp);
 ddt_T_3 = 2*((PAR.k/PAR.rho*PAR.Cp)/PAR.dz^2)*(T_2-T_3) + (PAR.rhod*PAR.dHdes/(PAR.rho*PAR.Cp))*c_3*(-PAR.fa*exp(-PAR.Ea/(PAR.R*T_3)))-2*(PAR.h/(PAR.rho*PAR.Cp*PAR.dz))*(T_3-Tb)- PAR.F*PAR.hrad*(T_3^4-Tc^4)/(PAR.V*PAR.rho*PAR.Cp);
 ddt_c_1 = c_1*(-PAR.fa*exp(-PAR.Ea/(PAR.R*T_1)));
 ddt_c_2 = c_2*(-PAR.fa*exp(-PAR.Ea/(PAR.R*T_2)));
 ddt_c_3 = c_3*(-PAR.fa*exp(-PAR.Ea/(PAR.R*T_3)));

 dXdt = [ddt_T_1; ddt_T_2; ddt_T_3; ddt_c_1; ddt_c_2; ddt_c_3];
end