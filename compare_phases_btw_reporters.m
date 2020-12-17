function compare_phases_btw_reporters()

    fluorescence_extraction_method = 'ilastik'; % 'whole_z_stack';
    thresh = 0.5;

    addpath DABEST/utilities/panel_2.6/demo/;
    addpath DABEST/Plotting/;
    addpath DABEST/Utilities/;    
    addpath export_fig/
    addpath my_xticklabels/

    
    [FaceColor, EdgeColor] = generate_cmh_colors();
    

    dirnames = repmat({'data/'},[1,4]);
    strainnames = {'HML620','GR1395','HML699','HML692'};
    listnames = {'HML620_list.txt', 'GR1395_list.txt','HML699_list.txt','HML692_list.txt'};
    reporter_names = { 'zk180.5_{pro}\newline::GFP-pest','mlt-10_{pro}\newline::GFP-pest', 'let-7_{pro}\newline::GFP-pest','dpy-6_{pro}\newline::GFP-pest'};

    no_strains = length(dirnames);

    % this string_array is used for defining categories in the plotSpread.m
    % function
    string_array = {};
    for ii = 1:no_strains;string_array = {string_array{:},num2str(ii)};  end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read the peak stats for each WT strain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:no_strains
        peak_stats{ii} = get_peak_stats(dirnames{ii},strainnames{ii}, listnames{ii}, ...
            fluorescence_extraction_method, thresh);
    end

    do_printing = 1;

    fig_position = [0.1 0.1 0.2 0.40];
    axis_position = [0.35 0.2 0.6 0.7];

    ms = 14; fs = 19; lw = 1.3;
    tl  = [0.02 0.02];
    
    fig = figure(1); clf;
    set(fig,'color', 'w', 'units', 'normalized', 'position', fig_position);    
    
    for ii = 1:no_strains
        
        if strfind(reporter_names{ii}, 'zk180.5')
            color_index = 4;
        elseif strfind(reporter_names{ii}, 'mlt-10')
            color_index = 2;
        elseif strfind(reporter_names{ii}, 'let-7')
            color_index = 3;
        elseif strfind(reporter_names{ii}, 'dpy-6')
            color_index = 1;
        end
        hold on;

        ind = ~isnan(~isnan(peak_stats{ii}.L1_molttimes)) & ...
            ~isnan(peak_stats{ii}.L2_peaktimes) & ...
            ~isnan(~isnan(peak_stats{ii}.L2_molttimes)); 
        
        data2plot = mod(2*pi*(peak_stats{ii}.L2_peaktimes(ind) - peak_stats{ii}.L1_molttimes(ind))...
            ./(peak_stats{ii}.L2_molttimes(ind) - peak_stats{ii}.L1_molttimes(ind)), 2*pi);

        catIdx = repmat({num2str(ii)}, [length(data2plot),1]);    

        plotSpread([NaN(no_strains,1);data2plot'],'categoryIdx',...
            {string_array{:}, catIdx{:}},'xValues', ii, 'categoryMarkers',...
            repmat({'o'}, [1,no_strains]),'categoryColors',...
            FaceColor([4,2,3,1],:));

            replace_marker_color(fig, FaceColor(color_index,:),...
                EdgeColor(color_index,:), FaceColor(color_index,:),lw,ms);
    end
  
    

    ylabel('Peak phase in L2');
    set(gca,'position', axis_position);
    set(gca, 'box','off', 'tickdir', 'out', 'linewidth', lw, 'fontsize', fs, 'ticklength', tl);
    set(gca, 'Ylim', [0, 2*pi], 'ytick', [0 pi/2 pi 3*pi/2 2*pi], 'yticklabel', {'0','\pi/2','\pi','3\pi/2','2\pi'});
    set(gca, 'xtick', 1:length(reporter_names), 'xticklabel', reporter_names);
    
    ax = gca;
    ax.XTickLabelRotation = 30;
    
    set(gca,'xlim', [0.7 4.5]);
    
    % This flips X and Y Axis! 
    view([90 -90]);
    set(gca, 'xdir', 'reverse');
    
    if do_printing
        set(gca,'color', 'none');
        set(gcf,'color', 'none');
        export_fig(gcf, ['pics/'  'WT_reporter_phase_comparison.pdf'],'-pdf');        
        set(gca,'color', 'w');
        set(gcf,'color', 'w');
    end
    

end


%%%%%%%%%%%%%%%%%%% USED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function replace_marker_color(fig_handle, orig_color,  Edgecolor, Facecolor,lw,ms)
    set(findobj(fig_handle,'color',orig_color),'Color',Edgecolor ,'MarkerFaceColor',Facecolor,'Linewidth', lw,'Markersize', ms);
    set(findobj(fig_handle,'MarkerEdgeColor',orig_color),'MarkerEdgeColor',Edgecolor ,'MarkerFaceColor',Facecolor,'Linewidth', lw);
end