function outputs = filterS(S,ip)

mean = S(1,:);
mean(mean>ip.H2) = ip.H2;
S(1,:) = mean;

outputs = S;

end
    