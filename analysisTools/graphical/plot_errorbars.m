function plot_errorbars(x,y,err,style,bar_col,width_line)
%plot_errorbars(x,y,err,style,bar_col,width_line)
%x,y points with lines of length err either up, down, or both
%style 0: up and down, 1: just up, 2: just down
%bar_col is bar color (e.g., for black either 'b' or [0 0 0])

num_points=length(x);

switch style
    case 0 %up/down
          plot_err_up(x,y,err,bar_col, num_points);hold on
        plot_err_down(x,y,err,bar_col, num_points);
    case 1 %just up
        plot_err_up(x,y,err,bar_col,num_points);hold on
    case 2  %just down
        plot_err_down(x,y,err,bar_col, num_points);hold on
    otherwise
        warning('plot_errorbars: need style') %#ok<WNTAG>
        return
end

%Subfunctions called
function plot_err_up(x,y,err,bar_col, num_points)
    for ii=1:num_points
        plot([x(ii) x(ii)],[y(ii) y(ii)+err(ii)],'Color',bar_col, 'LineWidth',width_line);hold on;
    end %for i
end %plot_up function



function plot_err_down(x,y,err,bar_col, num_points)
    for ii=1:num_points
        plot([x(ii) x(ii)],[y(ii) y(ii)-err(ii)],'Color', bar_col,'LineWidth',width_line);hold on;
    end %for i
end %plot_up function

end %main function
