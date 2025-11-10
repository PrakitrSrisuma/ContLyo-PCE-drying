function fout = OUTPUT_2Dry(t,X_in,PAR,Tb_profile,Tc_profile)

 T_1 = X_in(1,:);
 T_2 = X_in(2,:);
 T_3 = X_in(3,:);
 c_1 = X_in(4,:);
 c_2 = X_in(5,:);
 c_3 = X_in(6,:);
 
 Tb = cal_Tb(t,Tb_profile);
 Tc = cal_Tc(t,Tc_profile);
 
 c_avg = (c_1+c_2+c_2)*(1/3);
 T_avg = (T_1+T_2+T_2)*(1/3);
 
 fout = [c_avg; T_avg];
end