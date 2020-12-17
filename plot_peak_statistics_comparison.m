function [] = plot_peak_statistics_comparison(dirnames, strainnames,listnames,...
    genotypes, comparison_name)
%
% DESCRIPTION                      
% this function makes plots that compare different aspects of the peak
% statistics for two strains 
% Used for Figure 5 of Stec et al., Current Biology 2020
% 
% 
% EXAMPLE: 
% plot_peak_statistics_comparison({'data/', 'data/'}, {'GR1395', 'HML274'},{'GR1395_list.txt','HML274_list.txt'},{'WT','blmp-1(0)'}, 'mlt-10');
%
% plot_peak_statistics_comparison({'data/', 'data/'}, {'HML620', 'HML698'},{'HML620_list.txt','HML698_list.txt'},{'WT','blmp-1(0)'}, 'zk180.5');
%
% 
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    close all;
    % Check if input is correct
    if ~iscell(dirnames) || ~iscell(strainnames) || ~iscell(genotypes) || ~iscell(listnames)
        error('Input is wrong format. Must be cell arrays');
    end
    
    if length(dirnames) ~= length(genotypes)  ... 
            || length(dirnames) ~= length(strainnames) ...
            || length(genotypes) ~= length(strainnames) ...
            || length(genotypes) ~= length(listnames) ...
            || length(listnames) ~= length(strainnames)
        error('Input arrays must have the same lengths.');
    end
    
    % Load plotting routines for estimation plot
    addpath DABEST/utilities/panel_2.6/demo/;
    addpath DABEST/Plotting/;
    addpath DABEST/Utilities/;    
    
    
    
    %%%%%%%%%%%%%%%%% Design of the figures is set here
    close all;
    %thresh = 0.5;
        
    
    
    disp('Getting average stage lengths for WT to scale ...');
    % all WT animals
    WT_stage_lengths = get_average_larval_stage_length('data/', 'all_WT_list.txt');
    for ii = 1:size(WT_stage_lengths,1)
      avg_WT_stage_lengths(ii) = mean(WT_stage_lengths(ii,~isnan(WT_stage_lengths(ii,:))));
    end
    

    [FaceColor, EdgeColor] = generate_cmh_colors();
    % Marker and line properties    
    if length(dirnames) == 2 % means we are plotting WT and one other mutant only (e.g. for figure 7)
        
        if contains(comparison_name, 'dpy-6') || contains(comparison_name, 'dpy6')
            FaceColor = FaceColor([1,5],:);
            EdgeColor = EdgeColor([1,5],:);
            division_factor = 3; % this is the factor to make all axis the same in the figure
            lp_cutoff = 1.3;

        elseif contains(comparison_name, 'mlt-10') || contains(comparison_name, 'mlt10')
            FaceColor = FaceColor([2,6],:);
            EdgeColor = EdgeColor([2,6],:);
            division_factor = 120; % this is the factor to make all axis the same in the figure
            lp_cutoff = 0.3;
            
        elseif contains(comparison_name, 'let-7') || contains(comparison_name, 'let7')
            FaceColor = FaceColor([3,7],:);
            EdgeColor = EdgeColor([3,7],:);
            division_factor = 4; % this is the factor to make all axis the same in the figure
            lp_cutoff = 1;
                     
        elseif contains(comparison_name, 'zk180.5')
            FaceColor = FaceColor([4,8],:);
            EdgeColor = EdgeColor([4,8],:);
            division_factor = 80; % this is the factor to make all axis the same in the figure
            lp_cutoff = 0.3;         
        end
   else
        FaceColor = FaceColor([2,6],:);
        EdgeColor = EdgeColor([2,6],:);
        division_factor = 1;
   end
        
    % plots for individual stages
    fig_position = [0.2 0.2 0.9*[0.28 0.5]];
    axis_position = [0.2 0.2 0.67 0.7];
        
    ms = 12;
    fs = 28;
    lw = 2;
    tl  = [0.025 0.025];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    no_strains = length(dirnames);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read the peak stats for each strain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:no_strains
        peak_stats{ii} = get_peak_stats(dirnames{ii},strainnames{ii}, listnames{ii});
    end
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%% Compare L2 peaks heights
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. L2 peak height
    
    fig = figure(5); clf;
    % Modify figure properties
    set(fig, 'units', 'normalized', 'position', fig_position);

    % First gather all the data and identifiers
    % Generate identifiers for plotting routine
    data2plot = [];
    identifiers= {};
    for ii = 1:no_strains
        if isempty(data2plot)
            data2plot = peak_stats{ii}.L2_peakheights(~isnan(peak_stats{ii}.L2_peakheights))'/division_factor;
            identifiers = repmat({genotypes{ii}}, [1, length(data2plot)]);
        else
            data2plot = [data2plot; peak_stats{ii}.L2_peakheights(~isnan(peak_stats{ii}.L2_peakheights))'/division_factor];
            tmp = repmat({genotypes{ii}}, [1, length(peak_stats{ii}.L2_peakheights(~isnan(peak_stats{ii}.L2_peakheights)))]);
            identifiers = {identifiers{:}, tmp{:}};

        end
    end
    
    [ss] = FscatJit2(identifiers', data2plot, [0 35]);
    modify_figure_properties(fig, FaceColor, EdgeColor, lw, ms, fs,tl);

    % 
    ax = get(fig, 'Children'); % Gets all axes, ax(1) is the lower axis for the differences
    if strcmpi(comparison_name,'mlt-10')
        ax(1).YTick = [-20 -10 0];
    end
    
    set(gcf,'color', 'w');
    set(ax(:),'color', 'w');
    

%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%% Compare L2 peak phase 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig = figure(6); clf;
    % Modify figure properties
    set(fig, 'units', 'normalized', 'position', fig_position);

    % First gather all the data and identifiers
    % Generate identifiers for plotting routine 
        
    data2plot = [];
    identifiers= {};
    for ii = 1:no_strains
        ind = ~isnan(~isnan(peak_stats{ii}.L1_molttimes)) & ...
            ~isnan(peak_stats{ii}.L2_peaktimes) & ...
            ~isnan(~isnan(peak_stats{ii}.L2_molttimes)); 
        
        if isempty(data2plot)
%             data2plot = mod(2*pi*(peak_stats{ii}.L2_peaktimes(ind) - peak_stats{ii}.L1_molttimes(ind))...
%                 ./(peak_stats{ii}.L2_molttimes(ind) - peak_stats{ii}.L1_molttimes(ind)), 2*pi)';
            data2plot = (2*pi*(peak_stats{ii}.L2_peaktimes(ind) - peak_stats{ii}.L1_molttimes(ind))...
                ./(peak_stats{ii}.L2_molttimes(ind) - peak_stats{ii}.L1_molttimes(ind)))';

            identifiers = repmat({genotypes{ii}}, [1, length(data2plot)]);
        else
%             data2plot = [data2plot; mod(2*pi*(peak_stats{ii}.L2_peaktimes(ind) - peak_stats{ii}.L1_molttimes(ind))...
%                 ./(peak_stats{ii}.L2_molttimes(ind) - peak_stats{ii}.L1_molttimes(ind)), 2*pi)'];
            
            data2plot = [data2plot; (2*pi*(peak_stats{ii}.L2_peaktimes(ind) - peak_stats{ii}.L1_molttimes(ind))...
                ./(peak_stats{ii}.L2_molttimes(ind) - peak_stats{ii}.L1_molttimes(ind)))'];
            
            tmp = repmat({genotypes{ii}}, [1, length(mod(2*pi*(peak_stats{ii}.L2_peaktimes(ind) - peak_stats{ii}.L1_molttimes(ind))...
                    ./(peak_stats{ii}.L2_molttimes(ind) - peak_stats{ii}.L1_molttimes(ind)), 2*pi))]);
            identifiers = {identifiers{:}, tmp{:}};

        end
    end
    
    [ss] = FscatJit2(identifiers(data2plot>0)', data2plot(data2plot>0), [0 2*pi]);
    modify_figure_properties(fig, FaceColor, EdgeColor, lw, ms, fs,tl);

    % 
    ax = get(fig, 'Children'); % Gets all axes, ax(1) is the lower axis for the differences    
    ax(3).YTick = [0 pi  2*pi];
    ax(3).YTickLabel = {'0', '\pi', '2\pi'};
    
    ax(1).YTick = [-pi/4 0 pi/4];
    if strcmpi(comparison_name, 'dpy-6');
        ax(1).YTick = [-pi/4 0 pi/4];
        ax(1).YTickLabel = {'-\pi/4', '0', '\pi/4'};
    elseif strcmpi(comparison_name, 'mlt-10');
        ax(1).YTick = [-pi/2 0];
        ax(1).YTickLabel = {'-\pi/2', '0'};
    else
        ax(1).YTick = [-pi/4 0 pi/4];
        ax(1).YTickLabel = {'-\pi/4', '0', '\pi/4'};
    end

    set(gcf,'color', 'w');
    set(ax(:),'color', 'w');
    

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Compare L2 transcription onset phase
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(7); clf;
    % Modify figure properties
    set(fig, 'units', 'normalized', 'position', fig_position);

    % First gather all the data
    data2plot = [];
    identifiers = {};
    for ii = 1:no_strains
        onset_times = [];

        fid = fopen([dirnames{ii} listnames{ii}]);
        C1 = textscan(fid, '%s');
        fclose(fid);
        for jj = 1:length(C1{1,1})
            if ~strcmpi(C1{1,1}{jj}(1), '%')  % this allows to comment out animals in the list
                
                % Get the normalization of the duration 
                disp(['Reading worm: ' C1{1,1}{jj}]);
                worm = read_single_worm_molting_data(C1{1,1}{jj});
                % read molting times
                molts = get_molt_times(worm);
                
                stage_lengths = molts(2:4) - molts(1:3);
                stage_lengths = [NaN stage_lengths]; % We have no data on the L1 stage lengths
                
                onsets = get_peak_onset_offsets(dirnames{ii}, strainnames{ii}, listnames{ii}, comparison_name, jj, lp_cutoff);
                
                if ~isnan(onsets(2)) && ~isnan(stage_lengths(2)) 
                    onset_times = [onset_times ((onsets(2)-molts(1))/stage_lengths(2)*2*pi)];
                end
            end
        end
        if isempty(identifiers)
            identifiers = repmat({genotypes{ii}}, [1, length(onset_times)]);
            data2plot = onset_times;
        else
            data2plot = [data2plot, onset_times];
            tmp = repmat({genotypes{ii}}, [1, length(onset_times)]);
            identifiers = {identifiers{:}, tmp{:}};

        end
       
    end
    [ss] = FscatJit2(identifiers', data2plot', [0 2*pi]);
    modify_figure_properties(fig, FaceColor, EdgeColor, lw, ms, fs,tl);
    
    ax = get(fig, 'Children'); % Gets all axes, ax(1) is the lower axis for the differences    
    ax(3).YTick = [0 pi  2*pi];
    ax(3).YTickLabel = {'0', '\pi', '2\pi'};
    
    ax(1).YTick = [-pi/4 0 pi/4];
    if strcmpi(comparison_name, 'dpy-6');
        ax(1).YTick = [-pi/4 0 pi/4];
        ax(1).YTickLabel = {'-\pi/4', '0', '\pi/4'};
    elseif strcmpi(comparison_name, 'mlt-10');
        ax(1).YTick = [-pi/2 0];
        ax(1).YTickLabel = {'-\pi/2', '0'};
    else
        ax(1).YTick = [-pi/4 0 pi/4];
        ax(1).YTickLabel = {'-\pi/4', '0', '\pi/4'};
    end

       
    set(gcf,'color', 'w');
    set(ax(:),'color', 'w');
        



%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Compare L2 peak duration 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(8); clf;
    % Modify figure properties
    set(fig, 'units', 'normalized', 'position', fig_position);

    % First gather all the data
    data2plot = [];
    identifiers = {};
    for ii = 1:no_strains
        durations = [];

        fid = fopen([dirnames{ii} listnames{ii}]);
        C1 = textscan(fid, '%s');
        fclose(fid);
        for jj = 1:length(C1{1,1})
            if ~strcmpi(C1{1,1}{jj}(1), '%')  % this allows to comment out animals in the list
                
                % Get the normalization of the duration 
                disp(['Reading worm: ' C1{1,1}{jj}]);
                worm = read_single_worm_molting_data(C1{1,1}{jj});
                % read molting times
                molts = get_molt_times(worm);
                stage_lengths = molts(2:4) - molts(1:3);
                stage_lengths = [NaN stage_lengths];
                normalization_factor = stage_lengths(2)/avg_WT_stage_lengths(2);
                             
                [onsets, offsets] = get_peak_onset_offsets(dirnames{ii}, strainnames{ii}, listnames{ii}, comparison_name, jj, lp_cutoff);
                if ~isnan(offsets(2) - onsets(2)) && ~isnan(normalization_factor)
                    durations = [durations (offsets(2) - onsets(2))/normalization_factor];
                end
            end
        end
        if isempty(identifiers)
            identifiers = repmat({genotypes{ii}}, [1, length(durations)]);
            data2plot = durations;
        else
            data2plot = [data2plot, durations];
            tmp = repmat({genotypes{ii}}, [1, length(durations)]);
            identifiers = {identifiers{:}, tmp{:}};

        end
       
    end
    if contains(comparison_name, 'zk180.5')
        [ss] = FscatJit2(identifiers', data2plot', [0 20]);
    elseif contains(comparison_name, 'mlt-10')
        [ss] = FscatJit2(identifiers', data2plot', [0 15]);
    elseif contains(comparison_name, 'let-7')
        [ss] = FscatJit2(identifiers', data2plot', [0 18]);
    elseif contains(comparison_name, 'dpy-6')
        [ss] = FscatJit2(identifiers', data2plot', [0 10]);
    end
    modify_figure_properties(fig, FaceColor, EdgeColor, lw, ms, fs,tl);
    
    
    set(gcf,'color', 'w');
    set(ax(:),'color', 'w');
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Compare L2 cumulative output (area under the curve
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(8); clf;
    % Modify figure properties
    set(fig, 'units', 'normalized', 'position', fig_position);

    % First gather all the data
    data2plot = [];
    identifiers = {};
    for ii = 1:no_strains
        AuCs = [];

        fid = fopen([dirnames{ii} listnames{ii}]);
        C1 = textscan(fid, '%s');
        fclose(fid);
        for jj = 1:length(C1{1,1})
            if ~strcmpi(C1{1,1}{jj}(1), '%')  % this allows to comment out animals in the list
                
                % Get the normalization of the duration 
                disp(['Reading worm: ' C1{1,1}{jj}]);
                worm = read_single_worm_molting_data(C1{1,1}{jj});
                % read molting times
                molts = get_molt_times(worm);
                
                stage_lengths = molts(2:4) - molts(1:3);
                stage_lengths = [NaN stage_lengths];
                
                
                [~,~,AuC] = get_peak_onset_offsets(dirnames{ii}, strainnames{ii}, listnames{ii}, comparison_name, jj, lp_cutoff);
                
                if ~isnan(AuC(2)) 
                    AuCs = [AuCs AuC(2)];
                end
            end
        end
        if isempty(identifiers)
            identifiers = repmat({genotypes{ii}}, [1, length(AuCs)]);
            data2plot = AuCs;
        else
            data2plot = [data2plot, AuCs];
            tmp = repmat({genotypes{ii}}, [1, length(AuCs)]);
            identifiers = {identifiers{:}, tmp{:}};

        end
       
    end
    [ss] = FscatJit2(identifiers', data2plot', [0 300]);
    modify_figure_properties(fig, FaceColor, EdgeColor, lw, ms, fs,tl);
    
    ax = get(fig, 'Children'); % Gets all axes, ax(1) is the lower axis for the differences        
    set(gcf,'color', 'w');
    set(ax(:),'color', 'w');
    
            
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function modify_figure_properties(fig, FaceColor, EdgeColor, lw, ms, fs,tl)

    
    ax = get(fig, 'Children'); % Gets all axes, ax(1) is the lower axis for the differences

    set(ax([1,2,3]), 'box','off', 'tickdir', 'out', 'linewidth', lw, 'fontsize', fs, 'ticklength', tl);    
    ax(3).YLabel.String = '';     
    ax(3).XTickLabel = {''};     
    ax(3).YLabel.FontSize = fs; 
    ax(3).XTick = [1,2];
    ax(3).XLim = [0.5 3.5];
    
    % Find objects in upper axis (objects are scatter plots whose colors
    % can be changed
    CircleSize = 500;
    tt = ax(3).findobj;
    
    % doulbe mutant dots
    tt(4).MarkerEdgeColor = EdgeColor(2,:);
    tt(4).MarkerFaceColor = FaceColor(2,:);
    
    % WT dots
    tt(5).MarkerEdgeColor = EdgeColor(1,:);
    tt(5).MarkerFaceColor = FaceColor(1,:);
    
    tt(4).LineWidth = lw;tt(5).LineWidth = lw;
    tt(4).SizeData = CircleSize; tt(5).SizeData = CircleSize;

    
    % Change the appearance of the right most errorbar
    a = get(ax(2),'children');a.LineWidth = lw; a.MarkerSize = ms;
    % Change appearance of errorbars for the WT and double mutant dots
    a = get(ax(3),'children');    a(2).LineWidth = lw; a(2).MarkerSize = ms;
    
    % Change appearance of the dotted lines
    a  = findall(fig, 'linestyle', ':');
    a(1).LineStyle = '--'; a(2).LineStyle = '--';
    a(1).LineWidth = lw;a(2).LineWidth = lw;
    a(1).Color = [0.7 0.7 0.7];a(2).Color = [0.7 0.7 0.7];
    
end

%%%%%%%%%%%%%%%%%%% USED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function replace_marker_color(fig_handle, orig_color,  Edgecolor, Facecolor,lw,ms)
    set(findobj(fig_handle,'color',orig_color),'Color',Edgecolor ,'MarkerFaceColor',Facecolor,'Linewidth', lw,'Markersize', ms);
    set(findobj(fig_handle,'MarkerEdgeColor',orig_color),'MarkerEdgeColor',Edgecolor ,'MarkerFaceColor',Facecolor,'Linewidth', lw);
end
