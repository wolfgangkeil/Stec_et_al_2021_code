function molts = get_molt_times(worm)
%
%  FUNCTION molts = get_molt_times(worm)
%
% returns a vector for molt times rather than individual output arguments
% 
% molts(5) is first egg laying,  if scored
% molts(6) is exp_start_time
%
%
%
% see also read_single_worm_molting_data.m
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    molts = NaN*ones(1,6);

    if strcmpi(worm.synch_time, 'NA') || strcmpi(worm.synch_date, 'NA')...
            || strcmpi(worm.synch_time, 'NC') || strcmpi(worm.synch_date, 'NC') ...
                    || strcmpi(worm.mount_date, 'NA') || strcmpi(worm.mount_time, 'NA')...
                        || strcmpi(worm.mount_date, 'NC') || strcmpi(worm.mount_time, 'NC') % Means the worms weren't taken from a hatch-off 

        molts(6) = 0;
    else
        if strcmpi(worm.experiment_type, 'L3_arrested') %% whenever the experiment is with arrested worms, the stage lengths have to be adjusted later!!
            molts(6) = 0;
        else
            % Convert time from synch to experiment start to hours
            molts(6) = etime(datevec([worm.mount_date ' ' worm.mount_time]),datevec([worm.synch_date ' ' worm.synch_time]))/3600;
        end
    end

    % Add the frame time to this duration to determine the overall
    % hours to the molts
    if isfield(worm, 'l1_molt_time')
        if ~isempty(worm.l1_molt_time)
            if ~strcmpi(worm.l1_molt_time, 'NC')
                molts(1)  = molts(6) + to_seconds(worm.l1_molt_time)/3600; 
            else
                molts(1) = NaN;
            end
        else
             molts(1) = NaN;
        end
    else
         molts(1) = NaN;
    end

    % Add the frame time to this duration to determine the overall
    % hours to the molts
    if isfield(worm, 'l2_molt_time')
        if ~isempty(worm.l2_molt_time)
            if ~strcmpi(worm.l2_molt_time, 'NC')
                molts(2) = molts(6) + to_seconds(worm.l2_molt_time)/3600; 
            else
                molts(2) = NaN;
            end
        else
            molts(2) = NaN;
        end
    else
        molts(2) = NaN;
    end
    
    
    if isfield(worm, 'l3_molt_time')
        if ~isempty(worm.l3_molt_time)
            if ~strcmpi(worm.l3_molt_time, 'NC')
                molts(3) = molts(6) + to_seconds(worm.l3_molt_time)/3600;                               
            else
                molts(3) = NaN;
            end
        else
            molts(3) = NaN;
        end
    else
        molts(3) = NaN;
    end

    if isfield(worm, 'l4_molt_time')
        if ~isempty(worm.l4_molt_time)
            if ~strcmpi(worm.l4_molt_time, 'NC')
                molts(4) = molts(6) + to_seconds(worm.l4_molt_time)/3600;                               
            else
                molts(4) = NaN;
            end
        else
            molts(4) = NaN;
        end
    else
        molts(4) = NaN;
    end

    
    if isfield(worm, 'first_egg_time')
        if ~isempty(worm.first_egg_time)
            if ~strcmpi(worm.molts(5), 'NC')
                molts(5) = molts(6) + to_seconds(worm.molts(5))/3600;                               
            else
                molts(5) = NaN;
            end
        else
            molts(5) = NaN;
        end
    end
end
    
 