function graphics_setup(plot_type)

switch plot_type
case '1by1'
set(gcf, 'units', 'centimeters', 'Position',  [10, 6, 6, 4.5]);
set(gca,'fontsize',7,'XMinorTick','on','YMinorTick','on')

case '1by2'
set(gcf, 'units', 'centimeters', 'Position',  [10, 6, 10, 4.5]);
set(gca,'fontsize',7,'XMinorTick','on','YMinorTick','on')

case '1by3'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 15, 4.5]);
set(gca,'fontsize',7,'XMinorTick','on','YMinorTick','on')

case '2by3'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 17, 11]);
set(gca,'fontsize',7,'XMinorTick','on','YMinorTick','on')

case '1by4'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 20, 4.5]);
set(gca,'fontsize',7,'XMinorTick','on','YMinorTick','on')

end

return