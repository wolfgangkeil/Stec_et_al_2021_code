function window_adds = get_window_adds(strain)
%
%
%
%
% DESCRIPTION: This specifies the windows in which we are looking for a peak within a
% larval stage, for some markers (e.g. zk180.5 and mlt-10) peaks are actually slightly after the
% molt, so this means 
%
%
%%%%%%%%%%%%%%%%%% code by Wolfgang Keil, Rockefeller University 2018

    switch strain
        case {'HML620', 'HML696', 'HML697', 'HML698'}% zk180.5 fluo marker
            window_adds = [8, 10];
        case {'HML699','HML701','HML702'} % let-7
            window_adds = [0, 0];
        case {'GR1395','HML274','HML474','HML475'} % mlt-10
            window_adds = [0, 4];
        otherwise
            window_adds = [0, 0];
    end

end
