function plot_S_MeanSD(t,S,linestyle1,linestyle2)

mean = S(1,:)*100;
SD = sqrt(S(2,:))*100;

curve1 = mean + SD;
curve2 = mean - SD;
a = [t, fliplr(t)];
area = [curve1, fliplr(curve2)];
plot(t, mean, linestyle1{:});
hold on
fill(a, area, linestyle2{:});
ylabel({'Sublimation front (cm)'})
xlabel('Time (h)')

end