function stage_lengths = get_average_larval_stage_length(dirname, list_file)
%
%
% this function outputs the entire array of average stage lengths as a
% vector, vector entry is NaN, if it cannot be determined because molts
% could not be scored 
%
%
% EXAMPLE: stage_lengths = get_average_larval_stage_length('data/', 'GR1395_list.txt')
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~strcmpi(dirname(end),'/')
        dirname = [dirname '/'];
    end


    fid = fopen([dirname list_file]);
    C1 = textscan(fid, '%s');
    
    
    
    stage_lengths = NaN*ones(4, length(C1{1,1}));
    
    % Go over all the worms! 
    for ii  = 1:length(C1{1,1})
        % this enables skipping certain worms with a '%' sign in the last
        filename = C1{1,1}{ii};
        if ~strcmpi(filename(1), '%')
            worm = read_single_worm_molting_data(C1{1,1}{ii});
            
            if ~isempty(worm)
                molts = get_molt_times(worm);
            else
                molts = NaN*ones(1,6);
            end
            stage_lengths(2:4,ii) = molts(2:4) - molts(1:3);
            stage_lengths(1,ii) = molts(1) - molts(6);
        end
    end
    
end