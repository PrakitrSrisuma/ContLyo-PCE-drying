function plot_cw_cdf(c,Nhis,color,P0,Ctarget)

switch color
case '1a'
cc1 = [.2 1 .3];
cc2 = [0 0.5 0];

case '1b'
cc1 = [0 .8 .75];
cc2 = [0 0 1];

end

% h = cdfplot(c); 
% h.linewidth = 2;
% hold on
h = histogram(c,Nhis,'Normalization','cdf','FaceColor',cc1,'FaceAlpha',.3,'linewidth',.5); hold on
% [values, edges] = histcounts(c, Nhis,'Normalization', 'cdf');
% centers = (edges(1:end-1)+edges(2:end))/2;
% plot(centers, values, 'linewidth',2.5,'color',cc2)

ylabel('CDF')
xlabel('Concentration (wt/wt)'); hold on

plot([Ctarget Ctarget],[0 P0],':r','linewidth',1.5); hold on
cmin = min(c);
plot([cmin Ctarget],[P0 P0],':r','linewidth',1.5)
text(.04,.9,num2str(P0),'Units','normalized','FontSize', 7,'fontweight', 'bold','Color','r');

ax = ancestor(h, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.3f')

end