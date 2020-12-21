function [FitResults,FitError, baseline,peaks_properties,fit_params, fit_ranges] = ...
    fit_worm_peaks(worm, molts, WT_lengths, t_time, fluorescence,window_adds)
%
%
%
% DESCRIPTION: This function fits Gaussian peaks to a background-subtracted signal
%
%
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    fit_stages_individually = 1;
    manual_subranges = 1; % this will let the user select the time-range used for fitting,
                          % particularly helpful with zk180.5 marker
    
    fit_ranges = []; % will be filled in, if subranges are chosen manually
    plot_fit = 0;
    extra_peak_shape_param = 0; % we are only fitting Gaussians so no additional parameters
    
    
    % The following is just to determine the fit parameters
    % fit_parameters (all in one structure), are returned in function for
    % re-using!
    fit_params.baseline_autozero = 0; % means no background is assumed (the functions is fed with a signal 
    % that already has background substracted)
    fit_params.no_iterations = 300;
    fit_params.positions = [];
    fit_params.peak_variance = 10; %1.0 is default in routine
    fit_params.bipolar = 0;  % only positive peaks
    fit_params.fit_shape = 1; % 16 means fixed-position Gaussians
                              % 1 means Gaussians
                              % 10 means fixed-width Gaussians

    if ~isfield(worm, 'fluorescence_peaks')
        % then number of peaks is the number of molts that are scorable
        peak_vector = [0 0 0 0];
        if ~isnan(molts(1))
            peak_vector(1) = 1;
        end
        if ~isnan(molts(2))
            peak_vector(2) = 1;
        end
        if ~isnan(molts(3))
            peak_vector(3) = 1;
        end
        if ~isnan(molts(4))
            peak_vector(4) = 1;
        end
    else
        tmp = worm.fluorescence_peaks;
        for jj = 1:4
            peak_vector(jj) = round(str2double(tmp(jj)));
        end

    end

    % Determine peak start positions
    positions = [];
    peak_width = 5;%;h for a start should be good
    
    indiv_peaks_properties = [];
    
    
    if peak_vector(1) == 1
        if ~isnan(molts(1))

            % Place initial peak position at the maximum within the
            % given larval stage (if maximum is not at the edge)
            inds  = find(t_time > window_adds(1) & t_time <= molts(1) + window_adds(2));           
            stage_fluorescence = fluorescence(inds);
            stage_t_time = t_time(inds);
            [pks,locs] = findpeaks(stage_fluorescence);

            % this selects the biggest peak in the range
            if length(locs) > 1
                locs = locs(find(pks == max(pks),1));
            end
            
            positions = [positions t_time(inds(locs)) peak_width];


            if fit_stages_individually
                disp('Fitting L1 stage ...');

                if fit_params.fit_shape ==1
                    fixed_parameters =[];
                elseif fit_params.fit_shape ==16 % fixed position Gaussians
                    fixed_parameters = t_time(inds(locs));
                elseif fit_params.fit_shape ==17 % fixed width Gaussians
                    fixed_parameters = peak_width;
                else
                    error('fit_params.fit_shape must be either 1, 16, or 17.');
                end
                if manual_subranges
                    [min_range, max_range] = select_range_to_fit('L1', t_time,fluorescence, molts);     
                    
                    inds  = find(t_time > min_range & t_time <= max_range);           
                    stage_fluorescence = fluorescence(inds);                    
                    stage_t_time = t_time(inds);

                    fit_ranges = [fit_ranges; min_range, max_range];
                    
                    [pks,locs] = findpeaks(stage_fluorescence);
                    if length(locs) > 1
                        locs = locs(find(pks == max(pks),1));
                    end
                    initial_peak_pos = stage_t_time(locs);                        

                end
                % ACTUAL FITTING PROCEDURE FOR EACH LARVAL STAGE
                [FitResults,FitError, baseline] = peakfit([stage_t_time; stage_fluorescence],...
                    (min(stage_t_time) + max(stage_t_time))/2,(max(stage_t_time) - min(stage_t_time)),... 
                                1,fit_params.fit_shape,extra_peak_shape_param,...
                                fit_params.no_iterations,...
                                [initial_peak_pos, peak_width],fit_params.baseline_autozero, ...
                                fixed_parameters, plot_fit, fit_params.bipolar,[],...
                                fit_params.peak_variance);

                if isempty(indiv_peaks_properties)
                    % Initialize the individial peaks property matrix, only used if we try
                    % to fit each stage separately
                    indiv_peaks_properties = NaN*zeros(4,size(FitResults,2));
                    indiv_peaks_properties(:,1) = 1:4;
                end
                
                % Put the fit results in an array for later processing
                indiv_peaks_properties(1,2:end) = FitResults(1,2:end);
                        
                disp('...done.');
            end
                            
        
        else
            % Otherwise just put it in the middle of the stage, don't do
            % stage-specific fit
            positions = [positions WT_lengths(1)/2 peak_width];
        end
       
    end

    if peak_vector(2) == 1
        if ~isnan(molts(2))
            if ~isnan(molts(1))

                % Place initial peak position at the maximum within the
                % given larval stage (if maximum is not at the edge)
                inds  = find((t_time > (molts(1)+ window_adds(1))) ...
                                & (t_time <= (molts(2)+ window_adds(2))));

            else
                % assume worm is already in L2 when experiment starts
                inds  = find((t_time > window_adds(1)/2) ...
                                & (t_time <= (molts(2)+ window_adds(2))));
            end
            stage_fluorescence = fluorescence(inds);
            [pks,locs] = findpeaks(stage_fluorescence);

            if length(locs) > 1
                locs = locs(find(pks == max(pks),1));
            end
            positions = [positions t_time(inds(locs)) peak_width];

            if fit_stages_individually

                disp('Fitting L2 stage ...');
                if fit_params.fit_shape ==1
                    fixed_parameters =[];
                elseif fit_params.fit_shape ==16 % fixed position Gaussians
                    fixed_parameters = t_time(inds(locs));
                elseif fit_params.fit_shape ==17 % fixed width Gaussians
                    fixed_parameters = peak_width;
                else
                    error('fit_params.fit_shape must be either 1, 16, or 17.');
                end
                if manual_subranges
                    [min_range, max_range] = select_range_to_fit('L2', t_time,fluorescence, molts);     
                    inds  = find(t_time > min_range & t_time <= max_range);           
                    stage_fluorescence = fluorescence(inds);                    
                    stage_t_time = t_time(inds);
                    fit_ranges = [fit_ranges; min_range, max_range];

                    [pks,locs] = findpeaks(stage_fluorescence);
                    if length(locs) > 1
                        locs = locs(find(pks == max(pks),1));
                    end
                    initial_peak_pos = stage_t_time(locs);                        

                end

                % ACTUAL FITTING PROCEDURE FOR EACH LARVAL STAGE
                [FitResults,FitError, baseline] = peakfit([stage_t_time; stage_fluorescence],...
                    (min(stage_t_time) + max(stage_t_time))/2,(max(stage_t_time) - min(stage_t_time)),... 
                                1,fit_params.fit_shape,extra_peak_shape_param,...
                                fit_params.no_iterations,...
                                [initial_peak_pos, peak_width],fit_params.baseline_autozero, ...
                                fixed_parameters, plot_fit, fit_params.bipolar,[],...
                                fit_params.peak_variance);

                if isempty(indiv_peaks_properties)
                    % Initialize the individial peaks property matrix, only used if we try
                    % to fit each stage separately
                    indiv_peaks_properties = NaN*zeros(4,size(FitResults,2));
                    indiv_peaks_properties(:,1) = 1:4;
                end

                % Put the fit results in an array for later processing
                indiv_peaks_properties(2,2:end) = FitResults(1,2:end);
                
            else
                positions = [positions molts(2) - WT_lengths(2)/2 peak_width];
            end
        else
            positions = [positions WT_lengths(1) + WT_lengths(2)/2 peak_width];
        end
    end

    if peak_vector(3) == 1
        if ~isnan(molts(3))
            if ~isnan(molts(2))
                % Place initial peak position at the maximum within the
                % given larval stage (if maximum is not at the edge)
                inds  = find((t_time > (molts(2) + window_adds(1))) & ...
                                    (t_time <= (molts(3) + window_adds(2))));

            else
                % assume worm is already in L3 when experiment starts
                inds  = find((t_time > window_adds(1)/2) & ...
                                    (t_time <= (molts(3) + window_adds(2))));
            end
            stage_fluorescence = fluorescence(inds);
            [pks,locs] = findpeaks(stage_fluorescence);

            if length(locs) > 1
                locs = locs(find(pks == max(pks),1));
            end
            positions = [positions t_time(inds(locs)) peak_width];

            if fit_stages_individually
                disp('Fitting L3 stage ...');
                if fit_params.fit_shape ==1
                    fixed_parameters =[];
                elseif fit_params.fit_shape ==16 % fixed position Gaussians
                    fixed_parameters = t_time(inds(locs));
                elseif fit_params.fit_shape ==17 % fixed width Gaussians
                    fixed_parameters = peak_width;
                else
                    error('fit_params.fit_shape must be either 1, 16, or 17.');
                end
                if manual_subranges
                    [min_range, max_range] = select_range_to_fit('L3', t_time,fluorescence, molts);     
                    inds  = find(t_time > min_range & t_time <= max_range);           
                    stage_fluorescence = fluorescence(inds);                    
                    stage_t_time = t_time(inds);
                    fit_ranges = [fit_ranges; min_range, max_range];

                    [pks,locs] = findpeaks(stage_fluorescence);
                    if length(locs) > 1
                        locs = locs(find(pks == max(pks),1));
                    end
                    initial_peak_pos = stage_t_time(locs);                        

                end

                % ACTUAL FITTING PROCEDURE FOR EACH LARVAL STAGE
                [FitResults,FitError, baseline] = peakfit([stage_t_time; stage_fluorescence],...
                    (min(stage_t_time) + max(stage_t_time))/2,(max(stage_t_time) - min(stage_t_time)),... 
                                1,fit_params.fit_shape,extra_peak_shape_param,...
                                fit_params.no_iterations,...
                                [initial_peak_pos, peak_width],fit_params.baseline_autozero, ...
                                fixed_parameters,plot_fit, fit_params.bipolar,[],...
                                fit_params.peak_variance);

                if isempty(indiv_peaks_properties)
                    % Initialize the individial peaks property matrix, only used if we try
                    % to fit each stage separately
                    indiv_peaks_properties = NaN*zeros(4,size(FitResults,2));
                    indiv_peaks_properties(:,1) = 1:4;
                end

                % Put the fit results in an array for later processing
                indiv_peaks_properties(3,2:end) = FitResults(1,2:end);
            else
                positions = [positions molts(3) - WT_lengths(3)/2 peak_width];
            end
        else
            positions = [positions WT_lengths(1)+ WT_lengths(2) + WT_lengths(3)/2 peak_width];
        end
    end

    if peak_vector(4) == 1
        if ~isnan(molts(4))
            if ~isnan(molts(3))
                % Place initial peak position at the maximum within the
                % given larval stage (if maximum is not at the edge)
                inds  = find((t_time > (molts(3) + window_adds(1))) & ...
                                    (t_time <= (molts(4)+ window_adds(2))));

            else
                inds  = find((t_time > window_adds(1)/2) & ...
                                    (t_time <= (molts(4)+ window_adds(2))));
            end
            stage_fluorescence = fluorescence(inds);
            stage_t_time = t_time(inds);
            [pks,locs] = findpeaks(stage_fluorescence);

            if length(locs) > 1
                locs = locs(find(pks == max(pks),1));
            end
            positions = [positions t_time(inds(locs)) peak_width];

            if fit_stages_individually
                disp('Fitting L4 stage ...');
                if fit_params.fit_shape ==1
                    fixed_parameters =[];
                elseif fit_params.fit_shape ==16 % fixed position Gaussians
                    fixed_parameters = t_time(inds(locs));
                elseif fit_params.fit_shape ==17 % fixed width Gaussians
                    fixed_parameters = peak_width;
                else
                    error('fit_params.fit_shape must be either 1, 16, or 17.');
                end
                if manual_subranges
                    [min_range, max_range] = select_range_to_fit('L4', t_time,fluorescence, molts);     
                    inds  = find(t_time > min_range & t_time <= max_range);           
                    stage_fluorescence = fluorescence(inds);                    
                    stage_t_time = t_time(inds);

                    [pks,locs] = findpeaks(stage_fluorescence);
                    if length(locs) > 1
                        locs = locs(find(pks == max(pks),1));
                    end
                    initial_peak_pos = stage_t_time(locs);                        

                end

                % ACTUAL FITTING PROCEDURE FOR EACH LARVAL STAGE
                [FitResults,FitError, baseline] = peakfit([stage_t_time; stage_fluorescence],...
                    (min(stage_t_time) + max(stage_t_time))/2,(max(stage_t_time) - min(stage_t_time)),... 
                                1,fit_params.fit_shape,extra_peak_shape_param,...
                                fit_params.no_iterations,...
                                [initial_peak_pos, peak_width],fit_params.baseline_autozero, ...
                                fixed_parameters,plot_fit, fit_params.bipolar,[],...
                                fit_params.peak_variance);

                if isempty(indiv_peaks_properties)
                    % Initialize the individial peaks property matrix, only used if we try
                    % to fit each stage separately
                    indiv_peaks_properties = NaN*zeros(4,size(FitResults,2));
                    indiv_peaks_properties(:,1) = 1:4;
                end

                % Put the fit results in an array for later processing
                indiv_peaks_properties(4,2:end) = FitResults(1,2:end);
            
            else
                positions = [positions (molts(4) - WT_lengths(4)/2) peak_width];
            end
        else
            positions = [positions (WT_lengths(1)+ WT_lengths(2) + WT_lengths(3) + WT_lengths(4)/2) peak_width];
        end
    end

    fit_params.positions = positions;


    if exist('background_points', 'var')
        window_center = (min(background_points(:,1)) + max(background_points(:,1)))/2;
        window_width  = max(background_points(:,1)) - min(background_points(:,1));
    else
        window_center = (min(t_time) + max(t_time))/2;
        window_width  = max(t_time) - min(t_time);
    end

    fit_params.window_center = window_center;
    fit_params.window_width = window_width;
    fit_params.num_peaks = sum(peak_vector);


    if ~fit_stages_individually
        disp('Fitting peaks globally...');

        if fit_params.fit_shape ==1
            fixed_parameters =[];
        elseif fit_params.fit_shape ==16 % fixed position Gaussians
            fixed_parameters = positions(1:2:end);
        elseif fit_params.fit_shape ==17 % fixed width Gaussians
            fixed_parameters = positions(2:2:end);
        else
            error('fit_params.fit_shape must be either 1, 16, or 17.');
        end


        % ACTUAL FITTING PROCEDURE
        [FitResults,FitError, baseline] = peakfit([t_time; fluorescence],fit_params.window_center,fit_params.window_width,...
                        fit_params.num_peaks,fit_params.fit_shape,extra_peak_shape_param,fit_params.no_iterations,...
                        fit_params.positions,fit_params.baseline_autozero, ...
                        fixed_parameters,plot_fit, fit_params.bipolar,[], fit_params.peak_variance);



        % Put the fit results in an array for later processing
        peaks_properties = NaN*zeros(4,size(FitResults,2));
        peaks_properties(:,1) = 1:4;
        peaks_properties(peak_vector > 0,2:end) = FitResults(:,2:end);

        % Check if we want to fit individual peaks and replace with individual peak
        % parameters, if we actually fitted them
    else
       peaks_properties = NaN*zeros(4,size(FitResults,2));
       for ii = 1:4
            if ~isnan(indiv_peaks_properties(ii,2))
                peaks_properties(ii,:) = indiv_peaks_properties(ii,:);
            end
        end
    end
end


function [min_range, max_range] = select_range_to_fit(stage, t_fluo,avg_fluo_filt, molts)
%
%
% stage is a string,  'L1', 'L2' ...
%
%

    min_fluo = min(avg_fluo_filt);
    max_fluo = max(avg_fluo_filt);

    hFH = [];
    
    while isempty(hFH)
        
        figure(8);
        clf;
        title(['Select range (two points) to fit ' stage ' peaks to...']);
        hold on;
        plot(t_fluo,avg_fluo_filt, '-', 'linewidth', 2,'color', [0 1 0]);
        if ~isempty(molts) %plot molts
            for ii=1:4
                plot([molts(ii),molts(ii)], [min_fluo max_fluo],'--', 'color', [0.7 0.7 0.7]);
            end
        end
        hold off;

        hFH = impoly();
        if ~isempty(hFH) %check whether the user has put in anything
            range_points = hFH.getPosition;

            if size(range_points,1) == 2     % check whether used hasn't put in more than two points  
                min_range = range_points(1,1);
                max_range = range_points(2,1);
            else
                hFH = [];
            end
        end
    end
    close(gcf);

end
