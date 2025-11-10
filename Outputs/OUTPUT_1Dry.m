function fout = OUTPUT_1Dry(t,X_in,PAR,Tb_profile,Tc_profile)

 T_1 = X_in(1,:);
 T_2 = X_in(2,:);
 T_3 = X_in(3,:);
 S = X_in(4,:);
 
 Tb = cal_Tb(t,Tb_profile);
 Tc = cal_Tc(t,Tc_profile);
 Rp = PAR.Rp0 + PAR.Rp1*S/(1+S*PAR.Rp2);
 PwT = exp(-6139.9./T_1 + 28.8912);
 Pwc = 3;
 sw = switching_1Dry(PAR.H,S);
 Nw = (PwT-Pwc)/Rp;
 
 T_avg = (T_1+T_2+T_3)*(1/3);
 
 fout = [T_avg];
end