function outputs = cal_Tc(t,Tc_profile)

nt = length(t);
Tc = zeros(nt,1);

if length(Tc_profile) == 1 && ~isa(Tc_profile,'function_handle')
    for i = 1:nt
        Tc(i) = Tc_profile;
    end
elseif isa(Tc_profile,'function_handle')
    for i = 1:nt
        Tc(i) = Tc_profile(t(i));
    end
else
    for i = 1:nt
        Tc(i) = interp1(Tc_profile(:,2),Tc_profile(:,1),t(i));
    end
end

outputs = Tc;

return