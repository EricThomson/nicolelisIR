function gridfix(grid_color)
%gridfix(grid_color) %grid_color is 1x3 numeric array
%for older versions of matlab that have ugly-ass grid lines, make them pretty


set(gca,'gridlinestyle','-')
set(gca,'Xcolor',grid_color,'Ycolor',grid_color,'Zcolor',grid_color)
c=copyobj(gca,gcf);
set(c,'color','none','xcolor','k','xgrid','off', ...
    'ycolor','k','ygrid','off','zcolor','k','zgrid','off','Box','off');