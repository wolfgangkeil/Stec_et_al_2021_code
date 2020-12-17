function [t_time, fluorescence, success] = get_avg_fluorescence_timecourse(experiment_folder, worm_index,channel, imaging_interval) 
%
%
% PARAMETERS: 
% channel ... typically 'GFP'
% imaging_interval ... in minutes, in our data 15 minutes
%
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    fluorescence = [];
    t_time = [];
    success = 0;
        
    % This reads the actual fluorescence data, different methods used
    % depending on computation_method

    [fluorescence, t_time] = ...
        read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

    %%%%%%%%%%%%% SPECIAL TREATMENT FOR INTERRUPTED EXPERIMENTS
    
    if contains(experiment_folder, '20-Jul-2017')
        disp('Reading second part of the experiment')
        t_offset = 23.05;% Offset is 23h03min for this experiment
        
        % Read the second half of the experiment
        experiment_folder = strrep(experiment_folder,'20-Jul-2017', '21-Jul-2017');
        
        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];
    

    elseif contains(experiment_folder, '13-Mar-2018')
        disp('Reading second part of the experiment')
        t_offset = 8.75;% Offset is 08h45   min for this experiment
        
        % Read the second half of the experiment 
        experiment_folder = strrep(experiment_folder,'13-Mar-2018_2', '13-Mar-2018_3');
        
        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel, imaging_interval,computation_method, thresh);

        % Append the second part of the experiment
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];
    

    elseif contains(experiment_folder, '16-Mar-2018')
        disp('Reading second part of the experiment')
        t_offset = 23.25;% Offset is 08h45min for this experiment
        
        % Read the second half of the experiment (restarted 21-Jul-2017)
        experiment_folder = strrep(experiment_folder,'16-Mar-2018', '17-Mar-2018');

        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];


    elseif contains(experiment_folder, '23-Apr-2018')
        disp('Reading second part of the experiment')
        t_offset = 45.7833333;% Offset is 45h47min for this experiment
        
        % Read the second half of the experiment
        experiment_folder = strrep(experiment_folder,'23-Apr-2018', '25-Apr-2018');

        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment (returns empty array if
        % nothing is found 
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];

    
    %%%%%%%%%%%%% SPECIAL TREATMENT FOR INTERRUPTED EXPERIMENTS
    elseif contains(experiment_folder, '27-Apr-2018')
        disp('Reading second part of the experiment')
        t_offset = 88.0833333;% Offset is 88h05min for this experiment
        
        % Read the second half of the experiment
        experiment_folder = strrep(experiment_folder,'27-Apr-2018_2', '01-May-2018_1');

        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment (returns empty array if
        % nothing is found 
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];

    %%%%%%%%%%%%% SPECIAL TREATMENT FOR INTERRUPTED EXPERIMENTS
    elseif contains(experiment_folder, '04-May-2018')
        disp('Reading second part of the experiment')
        t_offset = 75.9;% Offset is 75h54min for this experiment
        
        % Read the second half of the experiment
        experiment_folder = strrep(experiment_folder,'04-May-2018_2', '07-May-2018_2');

        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment (returns empty array if
        % nothing is found 
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];
    %%%%%%%%%%%%% SPECIAL TREATMENT FOR INTERRUPTED EXPERIMENTS
    elseif contains(experiment_folder, '29-May-2018')
        disp('Reading second part of the experiment')
        t_offset = 96.6;% Offset is 96h36min for this experiment
        
        % Read the second half of the experiment
        experiment_folder = strrep(experiment_folder,'29-May-2018_1', '02-Jun-2018_2');

        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment (returns empty array if
        % nothing is found 
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];
        
        
    %%%%%%%%%%%%% SPECIAL TREATMENT FOR INTERRUPTED EXPERIMENTS
    elseif contains(experiment_folder, '06-Mar-2019')
        
        disp('Reading second part of the experiment')
        t_offset = 23.65;% Offset is 23h39min for this experiment
        
        % Read the second half of the experiment
        experiment_folder = strrep(experiment_folder,'06-Mar-2019_3', '07-Mar-2019_1');

        [fluorescence_2, t_time_2] = ...
            read_single_worm_avg_fluorescence_data(experiment_folder,worm_index,channel,imaging_interval, computation_method, thresh);

        % Append the second part of the experiment (returns empty array if
        % nothing is found 
        t_time = [t_time, (t_time_2 + t_offset)];
        fluorescence = [fluorescence, fluorescence_2];
        
    end
    
    
    
    if ~isempty(fluorescence)
        success = 1;
    end
end