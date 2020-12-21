    function plot_worm_overview(dirname,strain, list_file,channel, imaging_interval,bg_mode)
%
%   dirname ... directory where the strain_list is in
%   strain_list ... text-file with a list of worm-lineage files
%   imaging_interval ...  in minutes (15min in our data set)
%   channel ... typically 'GFP'
% 
%   bg_mode ... this is a string, determining, what background mode we are
%   using. Specify either 'spline_bg', 'constant_bg', or 'linear_bg'
%
%
%
% 
% DESCRIPTION: this function just plots all worms in a given list, no scaling etc. is done.
%              this is to give an overwiew over the data
% 
%
% see also get_avg_fluorescence_timecourse.m , fit_worm_peaks
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    addpath('peakfit');

    if nargin < 6
        bg_mode = 'constant_bg'; % default is a constant background, manual determination
    end

    
    %%%%%% set these  these flags to '1' forces re-doing of the individual steps
    compute_background = 1; 
    redo_fitting = 1;
        
    % Figure properties
    fs = 36;
    tl = [0.02 0.02];
    lw = 2;
    % seed the random number generator
    rand('seed', 1234);

        
    % FOR FIGURE
    overlay_fit = 1;  % This flag overlays the fitted peak onto the data
        
    fid = fopen([dirname list_file]);
    C1 = textscan(fid, '%s');
    
    %%% Sort out the worms that are commented out
    ind = [];
    for ii  = 1:length(C1{1,1})
        if ~strcmpi(C1{1,1}{ii}(1), '%')  % this allows to comment out animals in the txt-list by putting a '%' at the beginning of the line
            ind = [ind, ii];
        end
    end
    C1{1,1} = {C1{1,1}{ind}};
    
    % for plot
    no_cols = 4;
    % no rows
    no_rows = max(3,ceil(length(C1{1,1})/no_cols));    
    max_t = 0;
    max_fluo = 0;
    min_fluo = Inf;
    
    
    % Gets the average WT lethargi lengths (this is what we are going to
    % scale everything to for one of the plots)
    WT_lengths = get_average_larval_stage_length('data/', 'all_WT_list.txt');
  
    
    fh = figure(10);
    set(gcf, 'units', 'normalized', 'outerposition',[0.1 0.1 0.9 0.9]);
    clf;
    
    % this determines the y-axis limits for each strain, depending on
    % marker
    ylim = get_ylim(strain);
    
    % Go over the list of worms
    for ii = 1:length(C1{1,1})
        
        if ~strcmpi(C1{1,1}{ii}(1), '%')  % this allows to comment out animals in the list

            disp(['Reading worm: ' C1{1,1}{ii}]);
            worm = read_single_worm_molting_data(C1{1,1}{ii});

            tmp = strfind(C1{1,1}{ii},'_');
            tmp1 = strfind(C1{1,1}{ii},'.');

            worm_index = round(str2double(C1{1,1}{ii}(tmp(end)+1:tmp1(end)-1)));

            % read molting times
            molts = get_molt_times(worm);

            % generate a filename for the summary data, filename depends on
            % method with which fluorescence was extracted
            tmp1 = strfind(C1{1,1}{ii},'/');
            tmp2 = strfind(C1{1,1}{ii},'.');


            summary_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_prob_thresh_0.5.mat'];
            peaks_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_peaks_prob_thresh_0.5.mat'];
            bg_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_bg_prob_thresh_0.5.mat'];

            % This reads the fluorescence time traces after worm
            % segmentation with Ilastik
            if exist(summary_filename, 'file')
                load(summary_filename); 
            else
                break;
            end

            ax(ii) = subplot(no_rows,no_cols, ii);
            hold on;
            if strcmpi(worm.worm_development_arrested, 'no')
                mf_color = [0.3 1 0.3];
                ml_color = [0.3 0.7 0.3];
                fitline_color = [0.7 1 0.7];
            elseif strcmpi(worm.worm_development_arrested, 'yes')
                mf_color = [1 0.3 0.3];
                ml_color = [0.7 0.3 0.3];
                fitline_color = [1 0.7 0.7];
            else 
                mf_color = [0.6 0.6 1];
                ml_color = [0.3 0.3 0.7];
                fitline_color = [0.8 0.8 1];
            end
            plot(t_fluo, avg_fluo(1,:), 'o', 'Color',ml_color, 'MarkerFaceColor', mf_color, 'markersize',5);

            % Plot the molts
            plot([molts(1) ,molts(1)], ylim, '--', 'color', [0.5 0.5 0.5]);
            plot([molts(2) ,molts(2)], ylim, '--', 'color', [0.5 0.5 0.5]);
            plot([molts(3) ,molts(3)], ylim, '--', 'color', [0.5 0.5 0.5]);
            plot([molts(4) ,molts(4)], ylim, '--', 'color', [0.5 0.5 0.5]);
            hold off;

            if compute_background || ~exist(bg_filename, 'file')
                % Manually input a few points as background, a function is
                % fit to the points, depending on background mode
                if strcmpi(bg_mode, 'constant_bg')
                    [fitted_background, manual_bg_points, avg_fluo_filt] = ...
                                    determine_background(t_fluo,avg_fluo,'polynomial', 0,molts);
                elseif strcmpi(bg_mode, 'linear_bg')
                    [fitted_background, manual_bg_points, avg_fluo_filt] = ...
                                    determine_background(t_fluo,avg_fluo,'polynomial', 1,molts);
                elseif strcmpi(bg_mode, 'spline_bg')
                    [fitted_background, manual_bg_points, avg_fluo_filt] = ...
                                    determine_background(t_fluo,avg_fluo,'spline', [],molts);
                else
                    error('background mode needs to be ''constant_bg'',''linear_bg'' or ''spline_bg'' ');
                end
                save(bg_filename, 'fitted_background', 'manual_bg_points', 'avg_fluo_filt');            
            else
                disp('Reading previously selected background points...');
                load(bg_filename, 'fitted_background', 'manual_bg_points', 'avg_fluo_filt'); 
            end        

            if  redo_fitting || ~exist(peaks_filename,'file')
                window_adds = get_window_adds(strain);% window adds is good if the peaks are close to the molts!
                %
                    [FitResults,FitError, baseline,peaks_properties, FitParameters,FitRanges] = ...
                        fit_worm_peaks(worm, molts, WT_lengths,...
                            t_fluo, avg_fluo_filt - fitted_background, window_adds);
                % Save peaks_properties
                save(peaks_filename, 'FitResults', 'baseline','peaks_properties',...
                        'bg_mode', 'FitParameters', 'FitRanges');            
            else
                disp('Reading previously obtained fit to peaks...');
                load(peaks_filename,'FitResults', 'baseline','peaks_properties',... 
                        'bg_mode', 'FitParameters');            
            end
    
            if overlay_fit
                % overlay the plot
                figure(10);
                hold on;
                plot_single_worm_peaks(gca, t_fluo, peaks_properties, fitted_background, fitline_color);
                hold off;
            end

            if max_t < max(t_fluo)
                max_t = max(t_fluo);
            end
            if max_fluo < max(avg_fluo(:))
                max_fluo = max(avg_fluo(:));
            end
            if min_fluo >  min(avg_fluo(:))
                min_fluo = min(avg_fluo(:));
            end

            xlabel('time [h]'), ylabel('Fluo Intensity [u.a.]');
            if ~ismac
                title(strrep(C1{1,1}{ii}, '_', '\_'), 'fontsize', 6);
            else
                title(strrep(C1{1,1}{ii}, '_', '\_'));
            end
        end
    end
    
    
    
    % Sets all axis limits the same
    set(0,'CurrentFigure',fh);
    linkaxes(ax,'xy');
    ax(1).XLim = [0 max_t];
    ax(1).YLim = ylim;
    %ax(1).YLim = [min_fluo 400];
    set(gcf,'color', 'w');
    
    
end

% This specifies the y-axis limits for each strain (makes it easier to
% compare strains with the same marker 
function ylim = get_ylim(strain)
    switch strain
        case {'HML620', 'HML696', 'HML697', 'HML698'}% zk180.5 fluo marker
            ylim = [0, 3000];
        case {'GR1395','HML274','HML474', 'HML475'} % mlt-10
            ylim = [0, 7000];
        case {'HML692','HML693','HML694', 'HML695'} % dpy-6
            ylim = [100, 300];
        case {'HML699','HML701','HML702'} % let-7
            ylim = [100, 400];
        case {'HML661'} % lin-4
            ylim = [0, 1200];
        case {'HML619', 'HML625', 'HML635'} % mir-47
            ylim = [0, 1000];
        otherwise
            ylim = [0, 1000];
    end
end
