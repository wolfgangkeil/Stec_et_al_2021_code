function [stage_lengths, avg_used] = ...
         calculate_larval_stage_lengths(worm,avg_strain_stage_lengths)
%                               
%
% INPUT ARGUMENTS:
% worm ... structure returned from read_single_worm_lineage_data.m
%
% avg_strain_stage_lengths is returned from 
% 
% stage_lengths = get_average_lethargi_length_new(dirname, list_file)
% for ii = 1:size(stage_lengths,1); 
%   avg_avg_strain_stage_lengths(ii) = mean(stage_lengths(ii,~isnan(stage_lengths(ii,:))));
% end
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % read molting info for this particular worm    
    molts = get_lethargus_time_new(worm);
    
    avg_used = ones(1,4);
    % If the animal was been arrested in development, we have to 
    % renormalize the times, otherwise we cannot rescale
    if strcmpi(worm.experiment_type, 'L3_arrested')
        molts(2)  = molts(2) + avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2);
        molts(3)  = molts(3) + avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2);
        molts(4)  = molts(4) + avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2);
        molts(5)  = molts(5) + avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2);
    end
    
    
    % Get the length of the larval stages for this animal
    % if we can't determine it, we make assumptions based on overall developmental pace
    % until the first recorded times average length

    if ~isempty(molts(1))  % L1 molts is recorded
        stage_lengths(1) = molts(1);
        avg_used(1) = 0;
        if ~isempty(molts(2)) % means L1 and L2 molts are recorded
            stage_lengths(2) = molts(2) - molts(1);
            avg_used(2) = 0;
            if ~isempty(molts(3)) % means L1, L2 molts and L3 were recorded
                stage_lengths(3) = molts(3) - molts(2);
                avg_used(3) = 0;
                if ~isempty(molts(4)) % means all molts were recorded
                    stage_lengths(4) = molts(4) - molts(3);
                    avg_used(4) = 0;
                else % means L1, L2 molts and L3 were recorded but not L4 molt
                    % normalization factor based on overall developmental speed
                    normalization_factor = (molts(3))/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3)); 
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end
            else
                % means L1, L2 molts but not L3 molt was recorded
                if ~isempty(molts(4)) % means L1, L2 and L4 molts but not L3 molt was recorded
                    normalization_factor = (molts(4) - molts(2))/(avg_strain_stage_lengths(4) + avg_strain_stage_lengths(3)); 
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;

                else % means L1, L2 molts but any subsequent molts
                    normalization_factor = (molts(2))/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2)); 
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end                    
            end

        else % means L1 molt was recorded but L2 molt wasn't
            if ~isempty(molts(3)) % L1 and L3 molt was recorded so normalize to calculate the likely L2 molt
                normalization_factor = (molts(3)-molts(1))/(avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3));
                stage_lengths(2) = avg_strain_stage_lengths(2) * normalization_factor;
                stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;

                if ~isempty(molts(4)) % L1 L3 and L4 molt were recorded
                    stage_lengths(4) = molts(4) - molts(3);
                else % only L1, L3 recorded                         
                    normalization_factor = (molts(3))/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3));
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end
            else
                if ~isempty(molts(4)) % L1 and L4 molt were recorded so normalize to calculate the likely L2, and L3 molts
                    normalization_factor = (molts(4)-molts(1))/(avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3) + avg_strain_stage_lengths(4));
                    stage_lengths(2) = avg_strain_stage_lengths(2) * normalization_factor;
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(3) * normalization_factor;
                else
                    normalization_factor = (molts(1))/(avg_strain_stage_lengths(1));
                    stage_lengths(2) = avg_strain_stage_lengths(2) * normalization_factor;
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end                    
            end                
        end

    else % means no l1 molt was recorded

        if ~isempty(molts(2)) % This means the first recorded molt is L2
            % Calculate the normalization factor 
            normalization_factor = molts(2)/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2));
            % Define L1 length, and L2 length
            stage_lengths(1) = avg_strain_stage_lengths(1) * normalization_factor;
            stage_lengths(2) = avg_strain_stage_lengths(2) * normalization_factor;

            if ~isempty(molts(3)) % means L2 and L3 molt are recorded
                stage_lengths(3) = molts(3) - molts(2);
                if ~isempty(molts(4)) % means L2, L3 and L4 molt are recorded
                    stage_lengths(4) = molts(4) - molts(3);
                else % means L2, L3 we recorded but not L1 and L4
                    normalization_factor = molts(3)/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3));
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end
            else % means L2 but not L3 were recorded
                if ~isempty(molts(4)) % means only L2 and L4 molt are recorded
                    normalization_factor = (molts(4) - molts(2))/(avg_strain_stage_lengths(3) + avg_strain_stage_lengths(4));
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                else % means only L2 was recorded
                    normalization_factor = (molts(2))/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2));
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end                    
            end

        else

            if ~isempty(molts(3))
                % This means the first recorded molt is L3
                % Calculate the normalization factor 
                normalization_factor = molts(3)/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3));
                % Define L1 length, and L2 length, and L4 length
                stage_lengths(1) = avg_strain_stage_lengths(1) * normalization_factor;
                stage_lengths(2) = avg_strain_stage_lengths(2) * normalization_factor;
                stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;

                if ~isempty(molts(4)) % means L3 and L4 were recorded
                    stage_lengths(4) = molts(4) - molts(3);
                else
                    normalization_factor = molts(3)/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3));
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;
                end

            else

                if ~isempty(molts(4))
                    % This means the first recorded molt is L4
                    % Calculate the normalization factor 
                    normalization_factor = molts(3)/(avg_strain_stage_lengths(1) + avg_strain_stage_lengths(2) + avg_strain_stage_lengths(3) + avg_strain_stage_lengths(4));
                    % Define L1 length, and L2 length, and L4 length
                    stage_lengths(1) = avg_strain_stage_lengths(1) * normalization_factor;
                    stage_lengths(2) = avg_strain_stage_lengths(2) * normalization_factor;
                    stage_lengths(3) = avg_strain_stage_lengths(3) * normalization_factor;
                    stage_lengths(4) = avg_strain_stage_lengths(4) * normalization_factor;

                else
                    % This means no molt time was recorded,
                    % just assume, all is WT
                    stage_lengths(1) = avg_strain_stage_lengths(1);
                    stage_lengths(2) = avg_strain_stage_lengths(2);
                    stage_lengths(3) = avg_strain_stage_lengths(3);
                    stage_lengths(4) = avg_strain_stage_lengths(4);
                end
            end
        end
    end
end