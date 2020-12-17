function [onsets, offsets, AuC] = get_peak_onset_offsets(dirname, strain, list_file, marker_gene, ii,lp_cutoff)
%
% FUNCTION [onsets, offsets] = ...
%               get_peak_onset_offsets(dirname, strain, list_file, marker_gene, ii, lp_cutoff)
%
% PARAMETERS
% marker_gene can be either zk180.5, mlt-10
%
% lp_cutoff ... choose 0.3 (20min smoothing)
% (zk180.5, mlt-10)
%
%
%
% EXAMPLE: get_peak_onset_offsets('data/','GR1395','GR1395_list.txt','mlt-10', 1, 0.3);
%
%%%%%%%%%%%%%%%%%%%%%%% Coded by Wolfgang Keil, Institute Curie 2020

    verbose = 0;

    if exist([dirname list_file], 'file')
        fid = fopen([dirname list_file]);
        C1 = textscan(fid, '%s');
        fclose(fid);
    else
        warning('Cannot find list of animals.');
        return;
    end

    noise_level_cutoff = 0.98;% this is to ignore outliers in the noise levels stemming from rapid onsets of transcription
        
    % Marker and line properties      
    if contains(marker_gene, 'mlt-10')
        division_factor = 120; % this is the factor to make all axis the same in the final figure, units are arbitrary anyways
    elseif contains(marker_gene, 'zk180.5')
        division_factor = 80; % this is the factor to make all axis the same in the final figure, units are arbitrary anyways
    end    

    % generate a filename for the summary data, filename depends on
    % method with which fluorescence was extracted
    tmp1 = strfind(C1{1,1}{ii},'/');
    tmp2 = strfind(C1{1,1}{ii},'.');
    
    summary_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_prob_thresh_0.5.mat'];
    bg_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_bg_prob_thresh_0.5.mat'];
    
    if ~strcmpi(C1{1,1}{ii}(1), '%')  % this allows to comment out animals in the list
        worm = read_single_worm_molting_data(C1{1,1}{ii});
        tmp = strfind(C1{1,1}{ii},'_');
        tmp1 = strfind(C1{1,1}{ii},'.');
    end    
    
    molts = get_molt_times(worm);
 
    load(summary_filename, 'avg_fluo', 't_fluo'); % loads t_fluo, avg_fluo
    load(bg_filename, 'fitted_background', 'manual_bg_points', 'avg_fluo_filt'); % Loads background

    %%% THIS SUBTRACTS THE BACKGROUND!! --> ONly plots the oscillatory
    %%% component
    avg_fluo = avg_fluo - fitted_background;
    %%% normalize the signals so all panels of fig 7 have the same
    %%% y-axis
    avg_fluo = avg_fluo/division_factor;
    
    % Filter the signal 
    filter_shape = 'Gaussian';
    filter_type = 'highpass';
    frame_rate = 1/(t_fluo(2) - t_fluo(1)); %% in 1/h
    hp_cutoff = 0.3; %

    % interpolate with minimum value if any of the entries in avg_fluo are NaN
    % this is just for filtering, later, the NaNs are going to be
    % assigned back
    if sum(isnan(avg_fluo)) > 0
        t_OK = t_fluo(~isnan(avg_fluo));
        avg_fluo_OK = avg_fluo(~isnan(avg_fluo));
        avg_fluo = interp1(t_OK,avg_fluo_OK, t_fluo, 'spline', min(avg_fluo(:)));
    end

    %%% Gaussian Filter the signal
    [avg_fluo_filt, ~] = filter_signal(avg_fluo, frame_rate,filter_shape,filter_type, hp_cutoff);    
    avg_fluo_hp_filt = avg_fluo_filt';
    % Put the NaNs back
    avg_fluo_hp_filt(isnan(avg_fluo_hp_filt)) = NaN;

    tmp = sort(avg_fluo_hp_filt);
    noise_amp = tmp(round(noise_level_cutoff*length(tmp)));
    
    
    % Now, filter the signal again, now with lowpass cutoff specified as input parameter (this is the signal we
    % will then look for threshold crossings
    filter_shape = 'Gaussian';
    filter_type = 'lowpass';
    frame_rate = 1/(t_fluo(2) - t_fluo(1)); %% in 1/h

    % interpolate with minimum value if any of the entries in avg_fluo are NaN
    % this is just for filtering, later, the NaNs are going to be
    % assigned back
    if sum(isnan(avg_fluo)) > 0
        t_OK = t_fluo(~isnan(avg_fluo));
        avg_fluo_OK = avg_fluo(~isnan(avg_fluo));
        avg_fluo = interp1(t_OK,avg_fluo_OK, t_fluo, 'spline', min(avg_fluo(:)));
    end
    
    %%% Gaussian Filter the signal
    [avg_fluo_filt, ~] = filter_signal(avg_fluo, frame_rate,filter_shape,filter_type, lp_cutoff);    
    avg_fluo_lp_filt = avg_fluo_filt';
    % Put the NaNs back
    avg_fluo_lp_filt(isnan(avg_fluo_lp_filt)) = NaN;
    
    % Find upward threshold crossing times
    idxl = avg_fluo_lp_filt>=noise_amp;
    idxl(1) = 0;
    idx = find(idxl);
    yest = avg_fluo_lp_filt(idx-1)<noise_amp; 
    upward_crossings_t = t_fluo(idx(yest));  % gives times of upward threshold crossings   

    % Find downward threshold crossing times
    idxl = avg_fluo_lp_filt<noise_amp;
    idxl(1) = 0;
    idx2 = find(idxl);
    yest2 = avg_fluo_lp_filt(idx2-1)>=noise_amp; 
    downward_crossings_t = t_fluo(idx2(yest2)); % gives times of downward threshold crossings   

    if verbose
        figure(100);
        box off;
       
        plot(t_fluo,avg_fluo_lp_filt, 'linewidth', 2,'Color', [0.1 0.1 0.7 ]);
        hold on;
        %plot(t_fluo,avg_fluo_hp_filt);
        plot([min(t_fluo), max(t_fluo)], [noise_amp,noise_amp], '--', 'Color', [0.7 0.7 0.7]);
        plot([min(t_fluo), max(t_fluo)], [0 0], '-', 'Color', [0.7 0.7 0.7]);
        %plot([min(t_fluo), max(t_fluo)], -[noise_amp,noise_amp], '--', 'Color', [0.7 0.7 0.7]);

        plot(upward_crossings_t, avg_fluo_lp_filt(idx(yest)), 'o', 'Color', [0 0.7 0], 'linewidth',2);
        plot(downward_crossings_t, avg_fluo_lp_filt(idx2(yest2)), 'o', 'Color', [0.7 0 0], 'linewidth',2);

        hold off;
        set(gca, 'linewidth' ,2, 'fontsize', 16, 'tickdir', 'out');
        xlabel('time [h]','fontsize', 20);
        ylabel('Fluorescence [u.a.]','fontsize', 20);
        box off;
        drawnow;
    end
    
    ylim = get(gca, 'ylim');
    onsets = NaN*zeros(1,4);
    offsets = NaN*zeros(1,4);
    AuC = NaN*zeros(1,4);
    timestep = t_fluo(2)-t_fluo(1);
    
    if str2double(worm.fluorescence_peaks(1)) ~= 0
    
        % Get the onset and offset of the peaks for each larval stage
        if ~isnan(molts(1))

            if verbose; hold on; plot([molts(1),molts(1)], [0 ylim(2)], '-.','Color', [0.4 0.4 0.4]);hold off;end;
            ind = find(upward_crossings_t < molts(1),1, 'last');
            if ~isempty(ind)
                onsets(1) = upward_crossings_t(ind);
                if ~isnan(molts(2))
                    ind = find(downward_crossings_t > onsets(1) & downward_crossings_t < molts(2) ,1, 'first');
                else
                    ind = find(downward_crossings_t > onsets(1) ,1, 'first');
                end
                if ~isempty(ind)
                    offsets(1) = downward_crossings_t(ind);
                end
            end
        end  
    end
    
    if str2double(worm.fluorescence_peaks(2)) ~= 0
        % Get the onset and offset of the peaks for each larval stage
        if ~isnan(molts(2))
            if verbose; hold on; plot([molts(2),molts(2)], [0 ylim(2)], '-.','Color', [0.4 0.4 0.4]);hold off; end;

            if ~isnan(molts(1))
                ind = find(upward_crossings_t >  molts(1) & upward_crossings_t < molts(2),1, 'last');
            else
                ind = find(upward_crossings_t < molts(2),1, 'last');
            end

            if ~isempty(ind)
                onsets(2) = upward_crossings_t(ind);
                if ~isnan(molts(3))
                    ind = find(downward_crossings_t > onsets(2) & downward_crossings_t < molts(3) ,1, 'first');
                else
                    ind = find(downward_crossings_t > onsets(2) ,1, 'first');
                end
                if ~isempty(ind)
                    offsets(2) = downward_crossings_t(ind);
                end
            end
        end
    end
    
    if str2double(worm.fluorescence_peaks(3)) ~= 0
        % Get the onset and offset of the peaks for each larval stage
        if ~isnan(molts(3))
            if verbose; hold on; plot([molts(3),molts(3)], [0 ylim(2)], '-.','Color', [0.4 0.4 0.4]); hold off; end;
            if ~isnan(molts(2))
                ind = find(upward_crossings_t >  molts(2) & upward_crossings_t < molts(3),1, 'last');
            else
                ind = find(upward_crossings_t < molts(3),1, 'last');
            end
            if ~isempty(ind)
                onsets(3) = upward_crossings_t(ind);
                if ~isnan(molts(4))
                    ind = find(downward_crossings_t > onsets(3) & downward_crossings_t < molts(4) ,1, 'first');
                else
                    ind = find(downward_crossings_t > onsets(3) ,1, 'first');
                end
                if ~isempty(ind)
                    offsets(3) = downward_crossings_t(ind);
                end
            end
        end
    end
    if str2double(worm.fluorescence_peaks(4)) ~= 0
        % Get the onset and offset of the peaks for each larval stage
        if ~isnan(molts(4))
            if verbose; hold on; plot([molts(4),molts(4)], [0 ylim(2)], '-.','Color', [0.4 0.4 0.4]); hold off; end;
            if ~isnan(molts(3))
                ind = find(upward_crossings_t >  molts(3) & upward_crossings_t < molts(4),1, 'last');
            else
                ind = find(upward_crossings_t < molts(4),1, 'last');
            end
            if ~isempty(ind)
                onsets(4) = upward_crossings_t(ind);
                ind = find(downward_crossings_t > onsets(4),1, 'first');
                if ~isempty(ind)
                    offsets(4) = downward_crossings_t(ind);
                end
            end
        end
    end
    
    
    % Calculate the areas under the curve
    for jj = 1:4
        if (~isnan(onsets(jj)) && ~isnan(offsets(jj)) && str2double(worm.fluorescence_peaks(jj)) ~= 0)
            AuC(jj) = sum(avg_fluo_lp_filt(t_fluo>onsets(jj) & t_fluo<offsets(jj)))*timestep;
        end
    end
    
end