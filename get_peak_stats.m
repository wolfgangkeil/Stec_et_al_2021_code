function peak_stats = get_peak_stats(dirname,strain, list_file)
%
% FUNCTION peak_stats = get_peak_stats(dirname,strain, list_file)
%
% 
%
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    if nargin < 4
        thresh = [];
    end

    fid = fopen([dirname list_file]);
    C1 = textscan(fid, '%s');

    % Pre-define all the fields of the return structure
    peak_stats.L1_molttimes = [];
    peak_stats.L2_molttimes = [];
    peak_stats.L3_molttimes = [];
    peak_stats.L4_molttimes = [];
    
    peak_stats.L1_peaktimes = [];
    peak_stats.L2_peaktimes = [];
    peak_stats.L3_peaktimes = [];
    peak_stats.L4_peaktimes = [];

    peak_stats.L1_peakheights = [];
    peak_stats.L2_peakheights = [];
    peak_stats.L3_peakheights = [];
    peak_stats.L4_peakheights = [];
    
    peak_stats.L1_peakwidths = [];
    peak_stats.L2_peakwidths = [];
    peak_stats.L3_peakwidths = [];
    peak_stats.L4_peakwidths = [];

    peak_stats.L1_peakareas = [];
    peak_stats.L2_peakareas = [];
    peak_stats.L3_peakareas = [];
    peak_stats.L4_peakareas = [];

    % Unfortunately, we don't know how long L1 takes
    peak_stats.L2_lengths = [];
    peak_stats.L3_lengths = [];
    peak_stats.L4_lengths = [];

    disp(['Getting peak and molting statistis from strain ' strain '...']);
    
    % Go over the list of worms
    for ii  = 1:length(C1{1,1})
        if ~strcmpi(C1{1,1}{ii}(1), '%')  % this allows to comment out animals in the list

            worm = read_single_worm_molting_data(C1{1,1}{ii});

            tmp = strfind(C1{1,1}{ii},'_');
            tmp1 = strfind(C1{1,1}{ii},'.');

            worm_index = str2num(C1{1,1}{ii}(tmp(end)+1:tmp1(end)-1));

            % read molting times
            molts = get_molt_times(worm);

            % generate a filename for the summary data
            tmp1 = strfind(C1{1,1}{ii},'/');
            tmp2 = strfind(C1{1,1}{ii},'.');

            % Generate the appropriate file name
            summary_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_prob_thresh_0.5.mat'];
            peaks_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_peaks_prob_thresh_0.5.mat'];
            bg_filename = ['data/' strain '/' C1{1,1}{ii}(tmp1(end)+1:tmp2(end)-1) '_ilastik_bg_prob_thresh_0.5.mat'];

            % 
            if exist(peaks_filename, 'file')
                load(peaks_filename); % loads peaks_properties

                % Pre-define all the fields of the return structure
                peak_stats.L1_peaktimes = [peak_stats.L1_peaktimes, peaks_properties(1,2)];
                peak_stats.L2_peaktimes = [peak_stats.L2_peaktimes, peaks_properties(2,2)];
                peak_stats.L3_peaktimes = [peak_stats.L3_peaktimes, peaks_properties(3,2)];
                peak_stats.L4_peaktimes = [peak_stats.L4_peaktimes, peaks_properties(4,2)];

                peak_stats.L1_peakheights = [peak_stats.L1_peakheights, peaks_properties(1,3)];
                peak_stats.L2_peakheights = [peak_stats.L2_peakheights, peaks_properties(2,3)];
                peak_stats.L3_peakheights = [peak_stats.L3_peakheights, peaks_properties(3,3)];
                peak_stats.L4_peakheights = [peak_stats.L4_peakheights, peaks_properties(4,3)];

                peak_stats.L1_peakwidths = [peak_stats.L1_peakwidths, peaks_properties(1,4)];
                peak_stats.L2_peakwidths = [peak_stats.L2_peakwidths, peaks_properties(2,4)];
                peak_stats.L3_peakwidths = [peak_stats.L3_peakwidths, peaks_properties(3,4)];
                peak_stats.L4_peakwidths = [peak_stats.L4_peakwidths, peaks_properties(4,4)];

                peak_stats.L1_peakareas = [peak_stats.L1_peakareas, peaks_properties(1,5)];
                peak_stats.L2_peakareas = [peak_stats.L2_peakareas, peaks_properties(2,5)];
                peak_stats.L3_peakareas = [peak_stats.L3_peakareas, peaks_properties(3,5)];
                peak_stats.L4_peakareas = [peak_stats.L4_peakareas, peaks_properties(4,5)];

                
                % molts array values are NaN if molt hasn't been scored
                peak_stats.L1_molttimes = [peak_stats.L1_molttimes, molts(1)];
                peak_stats.L2_molttimes = [peak_stats.L2_molttimes, molts(2)];
                peak_stats.L3_molttimes = [peak_stats.L3_molttimes, molts(3)];
                peak_stats.L4_molttimes = [peak_stats.L4_molttimes, molts(4)];


                % Unfortunately, our data doesn't tell us  how long L1 takes
                peak_stats.L2_lengths = [peak_stats.L2_lengths, molts(2)-molts(1)];
                peak_stats.L3_lengths = [peak_stats.L3_lengths, molts(3)-molts(2)];
                peak_stats.L4_lengths = [peak_stats.L4_lengths, molts(4)-molts(3)];
            else
                disp(['Couldn''t find ' peaks_filename '! Analyze fluorescence time traces with analyze_avg_fluorescence.m!']);
            end
        end
    
    end    
    disp('...done');

end