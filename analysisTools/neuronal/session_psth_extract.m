function [psth, bin_edges]=session_psth_extract(event_times, time_window, bin_width, plot_on)
%[psth_all, psths_sum, {bin_centers}] = ...
% session_psth_extract(event_times, ref_times, time_window, bin_width, plot_on)
%bin_centers, spike_times are optional outputs
%
%Pull a single session-long psth for a set of event times.
%
%Inputs:
%event_times: numerical array of times of events for psth (usually spike times)
%ref_times: numerical array of reference events (e.g., beam break). 
%time_window: window to use to calculate the psth (e.g., [0 1000])
%bin_width: bin width for psth
%plot_on: 1 if you want to plot the psth
%
%Outputs:
%psth: psth for trial
%bin_edges: 1 x numBinswhat bin centers were used 


bin_edges=[time_window(1): bin_width  :time_window(2)];
psth = histc(event_times,bin_edges);
   
if plot_on
    h=bar(bin_edges, psth, 1);axis tight;hold on
    set(h,'FaceColor','k', 'EdgeColor','k')
    box on
    xlabel('Time');
    ylabel('Number of Events');shg
end

