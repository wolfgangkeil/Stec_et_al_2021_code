function compare_molting_btw_genotypes()
    addpath DABEST/utilities/panel_2.6/demo/;
    addpath DABEST/Plotting/;
    addpath DABEST/Utilities/;    
    addpath export_fig/

    close all;
    
    
    do_printing = 1;
    no_genotypes = 4;
    
    stage = 'L2';

    
    fig_position = [0.1 0.3 0.3 0.3];

    ms = 10; fs = 17; lw = 1.1;
    tl  = [0.015 0.015];

    
    %Design of the figures is set here
    Colors = define_marker_colors(no_genotypes, 'hammell');
    % Marker and line properties
%     ms = 12;
%     fs = 26;
%     lw = 2;
%     tl  = [0.025 0.025];

    % load WT experiments 
    % mlt-10 reporter
    WT_stage_lengths = get_average_lethargi_length_new('data/', 'GR1395_list.txt');
    % dpy-6 reporter
    WT_stage_lengths = [WT_stage_lengths, get_average_lethargi_length_new('data/', 'HML692_list.txt')];
    % zk180.5 reporter
    WT_stage_lengths = [WT_stage_lengths, get_average_lethargi_length_new('data/', 'HML620_list.txt')];
    WT_stage_lengths = [WT_stage_lengths, get_average_lethargi_length_new('data/', 'HML696_list.txt')];
    % lin-4 reporter
    WT_stage_lengths = [WT_stage_lengths, get_average_lethargi_length_new('data/', 'HML661_list_cont_develop_from_L1.txt')];
    % let 7 reporter
    WT_stage_lengths = [WT_stage_lengths, get_average_lethargi_length_new('data/', 'HML699_list.txt')];
    % mir-47 reporter
    WT_stage_lengths = [WT_stage_lengths, get_average_lethargi_length_new('data/', 'HML619_list.txt')];

    
    % load blmp-1 experiments 
    % mlt-10 reporter
    blmp1_stage_lengths = get_average_lethargi_length_new('data/', 'HML274_list.txt');
    % dpy-6 reporter
    blmp1_stage_lengths = [blmp1_stage_lengths, get_average_lethargi_length_new('data/', 'HML693_list.txt')];
    % zk180.5 reporter
    blmp1_stage_lengths = [blmp1_stage_lengths, get_average_lethargi_length_new('data/', 'HML697_list.txt')];
    % mir-47 reporter
    blmp1_stage_lengths = [blmp1_stage_lengths, get_average_lethargi_length_new('data/', 'HML625_list.txt')];

    % load elt-3 single mutant experiments 
    % mlt-10 reporter
    elt3_stage_lengths = get_average_lethargi_length_new('data/', 'HML475_list.txt');
    % dpy-6 reporter
    elt3_stage_lengths = [elt3_stage_lengths, get_average_lethargi_length_new('data/', 'HML694_list.txt')];
    
    % load double mutant experiments 
    % mlt-10 reporter
    double_stage_lengths = get_average_lethargi_length_new('data/', 'HML474_list.txt');
    % dpy-6 reporter
    double_stage_lengths = [double_stage_lengths, get_average_lethargi_length_new('data/', 'HML695_list.txt')];
    % zk180.5 reporter
    double_stage_lengths = [double_stage_lengths, get_average_lethargi_length_new('data/', 'HML698_list.txt')];
    %let-7 reporter
    double_stage_lengths = [double_stage_lengths, get_average_lethargi_length_new('data/', 'HML701_list.txt')];
    % mir-47 reporter
    double_stage_lengths = [double_stage_lengths, get_average_lethargi_length_new('data/', 'HML635_list.txt')];
    
    
    % Generate identifiers for plotting routine
    
    WT_identifiers = repmat({'WT'}, [1, length(WT_stage_lengths)]);
    blmp1_identifiers = repmat({'blmp-1(0)'}, [1, length(blmp1_stage_lengths)]);
    elt3_identifiers = repmat({'elt-3(0)'}, [1, length(elt3_stage_lengths)]);
    double_identifiers = repmat({'blmp-1(0); elt-3(0)'}, [1, length(double_stage_lengths)]);

    identifiers = {WT_identifiers{:}, blmp1_identifiers{:}, elt3_identifiers{:}, double_identifiers{:}}';
    
    if strcmpi(stage, 'L2')
        data = [WT_stage_lengths(2,:), blmp1_stage_lengths(2,:), elt3_stage_lengths(2,:), double_stage_lengths(2,:)]';
    elseif strcmpi(stage, 'L3')
        data = [WT_stage_lengths(3,:), blmp1_stage_lengths(3,:), elt3_stage_lengths(3,:), double_stage_lengths(3,:)]';
    elseif strcmpi(stage, 'L4')
        data = [WT_stage_lengths(4,:), blmp1_stage_lengths(4,:), elt3_stage_lengths(4,:), double_stage_lengths(4,:)]';
    end
    
    [ss] = FscatJit2(identifiers, data);
    
    disp('==============================')
    disp('Number of animals per genotype')
    disp(['WT L2 stage: ' num2str(sum(~isnan((WT_stage_lengths(2,:)))))]);
    disp(['WT L3 stage: ' num2str(sum(~isnan((WT_stage_lengths(3,:)))))]);
    disp(['WT L4 stage: ' num2str(sum(~isnan((WT_stage_lengths(4,:)))))]);
    disp('------')
    
    disp(['blmp1 mut L2 stage: ' num2str(sum(~isnan((blmp1_stage_lengths(2,:)))))]);
    disp(['blmp1 mut L3 stage: ' num2str(sum(~isnan((blmp1_stage_lengths(3,:)))))]);
    disp(['blmp1 mut L4 stage: ' num2str(sum(~isnan((blmp1_stage_lengths(4,:)))))]);
    disp('------')

    
    disp(['elt3 mut L2 stage: ' num2str(sum(~isnan((elt3_stage_lengths(2,:)))))]);
    disp(['elt3 mut L3 stage: ' num2str(sum(~isnan((elt3_stage_lengths(3,:)))))]);
    disp(['elt3 mut L4 stage: ' num2str(sum(~isnan((elt3_stage_lengths(4,:)))))]);
    disp('------')

    disp(['double mut L2 stage: ' num2str(sum(~isnan((double_stage_lengths(2,:)))))]);
    disp(['double mut L3 stage: ' num2str(sum(~isnan((double_stage_lengths(3,:)))))]);
    disp(['double mut L4 stage: ' num2str(sum(~isnan((double_stage_lengths(4,:)))))]);

    % Modify figure properties
    set(figure(1),'color', 'w', 'units', 'normalized', 'position', fig_position);
    ax = get(figure(1), 'Children'); % Gets all axes, ax(1) is the lower axis for the differences

    set(ax([1,4]), 'box','off', 'tickdir', 'out', 'linewidth', lw, 'fontsize', fs);
    set(ax([1,4]), 'XLim', [0.5 4.7]);    
    set(ax(4), 'YLim', [0 30]);
    
    
    set(ax(1), 'position', [0.12 0.05 0.88 0.3]);
    set(ax(4), 'position', [0.12 0.5 0.88 0.45]);
    
    ax(4).YLabel.String = [stage ' duration [h]']; 
    ax(1).YLabel.String = '\Delta value [h]'; 
    ax(1).YLabel.FontSize = fs; 
    ax(4).YLabel.FontSize = fs;
    ax(4).XTickLabel = {''};

    set(ax(1), 'xcolor','w');
        
    % Find objects in upper axis (objects are scatter plots whose colors
    % can be changed
    MEColor = [0.4 0.4 0.4];
    MFColor = [0.7 0.7 0.7];
    CircleSize = 120;
    tt = ax(4).findobj;
    
    
    % Tweak the circles and colors
    for ii = 4:7
        tt(ii).MarkerFaceColor = MFColor;
        tt(ii).MarkerEdgeColor = MEColor;
        tt(ii).LineWidth = lw;
        tt(ii).SizeData = CircleSize;
        
    end
    % Change appearance of dotted line
    a = get(ax(1),'children');
    a(1).Color = [0.7 0.7 0.7];
    a(1).LineWidth = lw;
    a(1).LineStyle = '--';
    
    %change the appearance of the errorbar-plot
    a(2).MarkerSize = ms;
    a(2).MarkerFaceColor = [0.7 0.7 0.7];
    a(2).MarkerEdgeColor = [0.4 0.4 0.4];
    
    if do_printing
        export_fig(['pics/molting_comparison_btw_genotypes' stage '.pdf']);
    end

%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%% Compare L1 lengths
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 1. WT
%     fig = figure(1);
%     clf;
%     set(gcf,'units','normalized', 'position', fig_position);
%     set(gca, 'position', axis_position);
%     
%     tmp = WT_stage_lengths(1,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'WT'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 1, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     hold on;
%     
%     % 2. blmp-1(0)
%     tmp = blmp1_stage_lengths(1,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'blmp-1(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 2, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     % 3. elt-3
%     tmp = elt3_stage_lengths(1,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'elt-3(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 3, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     % 4. double
%     tmp = double_stage_lengths(1,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'double'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 4, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     hold off;
%     
%     % Make markers look nice
%     for ii = 1:no_genotypes
%         replace_marker_color(fig, Colors.Face{ii},Colors.Edge{ii}, Colors.Face{ii},lw,ms);
%     end
%     
%     axis([0.5 4.5 0 35]);
%     set(gca, 'tickdir', 'out', 'fontsize', fs, 'ticklength', tl , 'linewidth', lw,...
%         'xtick', [1,2,3,4], 'xticklabel', {''});
%     ylabel('L1 stage length [h]');
%     %set(gca,'XTickLabelRotation',20);
% 
% 
%     if do_printing
%         set(gcf,'color', 'w');
%         export_fig(fig, ['pics/' 'stage_lengths_comparison_L1.pdf']);        
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%% Compare L2 lengths
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 1. WT
%     fig = figure(2);
%     clf;
%     set(gcf,'units','normalized', 'position', fig_position);
%     set(gca, 'position', axis_position);
% 
%     tmp = WT_stage_lengths(2,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'WT'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 1, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     hold on;
%     
%     % 2. blmp-1(0)
%     tmp = blmp1_stage_lengths(2,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'blmp-1(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 2, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     % 3. elt-3(0)
%     tmp = elt3_stage_lengths(2,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'elt-3(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 3, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     % 4. double
%     tmp = double_stage_lengths(2,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'double'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 4, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     hold off;
%     
%     % Make markers look nice
%     for ii = 1:no_genotypes
%         replace_marker_color(fig, Colors.Face{ii},  Colors.Edge{ii}, Colors.Face{ii},lw,ms);
%     end
%         
%     axis([0.5 4.5 0 25]);
%     set(gca, 'tickdir', 'out', 'fontsize', fs, 'ticklength', tl , 'linewidth', lw,...
%         'xtick', [1,2,3,4], 'xticklabel', {''});
%     ylabel('L2 stage length [h]');
%     %set(gca,'XTickLabelRotation',20);
%     
%     if do_printing
%         set(gcf,'color', 'none');
%         export_fig(fig, ['pics/' 'stage_lengths_comparison_L2.pdf']);        
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%% Compare L3 lengths
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 1. WT
%     fig = figure(3);
%     clf;
%     set(gcf,'units','normalized', 'position', fig_position);
%     set(gca, 'position', axis_position);
% 
%     tmp = WT_stage_lengths(3,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'WT'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 1, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     hold on;
%     
%     % 2. blmp-1(0)
%     tmp = blmp1_stage_lengths(3,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'blmp-1(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 2, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     % 3. elt-3(0)
%     tmp = elt3_stage_lengths(3,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'elt-3(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 3, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     % 4. double
%     tmp = double_stage_lengths(3,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'double'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 4, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     hold off;
%     
%     % Make markers look nice
%     for ii = 1:no_genotypes
%         replace_marker_color(fig, Colors.Face{ii},  Colors.Edge{ii}, Colors.Face{ii},lw,ms);
%     end
%         
%     axis([0.5 4.5 0 25]);
%     set(gca, 'tickdir', 'out', 'fontsize', fs, 'ticklength', tl , 'linewidth', lw,...
%         'xtick', [1,2,3,4], 'xticklabel', {''});
%     ylabel('L3 stage length [h]');
%     %set(gca,'XTickLabelRotation',20);
%     
%     if do_printing
%         set(gcf,'color', 'w');
%         export_fig(fig, ['pics/' 'stage_lengths_comparison_L3.pdf']);        
%     end
%             
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%% Compare L4 lengths
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 1. WT
%     fig = figure(4);
%     clf;
%     set(gcf,'units','normalized', 'position', fig_position);
%     set(gca, 'position', axis_position);
% 
%     tmp = WT_stage_lengths(4,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'WT'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 1, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     hold on;
%     
%     % 2. blmp-1(0)
%     tmp = blmp1_stage_lengths(4,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'blmp-1(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 2, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     % 3. elt-3(0)
%     tmp = elt3_stage_lengths(4,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'elt-3(0)'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 3, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
% 
%     % 4. double
%     tmp = double_stage_lengths(4,:);
%     tmp = tmp(~isnan(tmp));
%     catIdx = repmat({'double'}, [length(tmp),1]);    
%     
%     plotSpread([NaN;NaN;NaN;NaN;tmp'],'categoryIdx',{'WT', 'blmp-1(0)', 'elt-3(0)', 'double', catIdx{:}},'xValues', 4, 'categoryMarkers',...
%         {'o', 'o', 'o','o'},'categoryColors',reshape([Colors.Face{1:no_genotypes}], [3,no_genotypes])');
%     
%     hold off;
%     
%     % Make markers look nice
%     for ii = 1:no_genotypes
%         replace_marker_color(fig, Colors.Face{ii},  Colors.Edge{ii}, Colors.Face{ii},lw,ms);
%     end
%         
%     axis([0.5 4.5 0 25]);
%     set(gca, 'tickdir', 'out', 'fontsize', fs, 'ticklength', tl , 'linewidth', lw,...
%         'xtick', [1,2,3,4], 'xticklabel', {''});
%     ylabel('L4 stage length [h]');
%     %set(gca,'XTickLabelRotation',20);
%     
%     if do_printing
%         set(gcf,'color', 'w');
%         export_fig(fig, ['pics/' 'stage_lengths_comparison_L4.pdf']);        
%     end
end

%%%%%%%%%%%%%%%%%%% USED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function replace_marker_color(fig_handle, orig_color,  Edgecolor, Facecolor,lw,ms)
    set(findobj(fig_handle,'color',orig_color),'Color',Edgecolor ,...
        'MarkerFaceColor',Facecolor,'Linewidth', lw,'Markersize', ms);
end


function axesoffwithlabels(h)
%AXESOFFWITHLABELS Make axes invisible but not the xlabel and ylabel.
%
%   AXESOFFWITHLABELS(H) makes axes invisible, keeping the x- and ylabel
%   with handle H.
%
%  Sample Usage
%    plot(rand(1,10))
%    h(1) = xlabel('x');
%    h(2) = ylabel('x');
%    axesoffwithlabels(h)
%
%   Thorsten.Hansen@psychol.uni-giessen.de  2018-08-08
set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(h, 'Color', 'k')
% get rid of the white ticks and tick labels, moving the labels closer to
% the axes
set(gca, 'XTick', []);
set(gca, 'YTick', []);
end