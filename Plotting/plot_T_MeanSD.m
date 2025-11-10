function plot_T_MeanSD(t,T,linestyle1,linestyle2)

mean = T(1,:);
SD = sqrt(T(2,:));

curve1 = mean + SD;
curve2 = mean - SD;
a = [t, fliplr(t)];
area = [curve1, fliplr(curve2)];
fill(a, area, linestyle2{:});
hold on;
plot(t, mean, linestyle1{:});
ylabel('Temperature (K)')
xlabel('Time (h)')

end