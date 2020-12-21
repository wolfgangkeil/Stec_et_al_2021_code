function [bg, background_points,avg_fluo_filt] = determine_background(t_fluo, avg_fluo,...
                    background_function, bg_polynomial_degree, molts)
%
% t_fluo ... time of the experiment, loaded from the analyis-file
% avg_fluo ... fluorescence values, extracted from experiment
% background_function ... a string specifying the interpolation of the
%                   backgruond, choose either 'polynomial' or 'spline' 
% bg_polynomial_degree ... degree of the polynomial, if specified as the
%                           background function
% molts ... they are passed on for plotting reasons (optional)
%
% constant_bg_value (optional) ... this is a constant background value that
% will be subtracted from the signal in case polynomial of 0 zero is
% selected, this essentially means no fitting of the bg, just a
% subtraction, this was used for mlt-10 and zk180.5 markers in Stec et al. 2020
% 
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if nargin < 4
        molts =[];
    end    

    %----------------- Filter settings for Gaussian low pass --------------
    filter_shape = 'Gaussian';
    filter_type = 'lowpass';
    frame_rate = 1/(t_fluo(2) - t_fluo(1)); %% in 1/h
    cut_off = 0.5; %% Means 0.5h smoothing

    % interpolate if any of the entries in avg_fluo are NaN
    if sum(isnan(avg_fluo)) > 0
        t_OK = t_fluo(~isnan(avg_fluo));
        avg_fluo_OK = avg_fluo(~isnan(avg_fluo));

        avg_fluo = interp1(t_OK,avg_fluo_OK, t_fluo, 'spline', 0);
    end            

    %%% Gaussian Filter the signal
    [avg_fluo_filt, ~] = filter_signal(avg_fluo, frame_rate,filter_shape,filter_type, cut_off);    
    avg_fluo_filt = avg_fluo_filt';
    
    min_fluo = 0;
    max_fluo = max(avg_fluo(:));
    
    
    hFH = [];
    while isempty(hFH)
        figure(8);
        clf;
        hold on;
        plot(t_fluo,avg_fluo, '-', 'linewidth', 0.5,'color', [0.5 1 0.5]);
        plot(t_fluo,avg_fluo_filt, '-', 'linewidth', 2,'color', [0 1 0]);
        if ~isempty(molts)
            plot([molts(1),molts(1)], [min_fluo max_fluo],'--', 'color', [0.7 0.7 0.7]);
            plot([molts(2),molts(2)], [min_fluo max_fluo],'--', 'color', [0.7 0.7 0.7]);
            plot([molts(3),molts(3)], [min_fluo max_fluo],'--', 'color', [0.7 0.7 0.7]);
            plot([molts(4),molts(4)], [min_fluo max_fluo],'--', 'color', [0.7 0.7 0.7]);
        end
        hold off;

        if strcmpi(background_function,'polynomial')
            if nargin < 4
                bg_polynomial_degree = 6;
            end

            if bg_polynomial_degree == 1
                title('Select background points for LINEAR background function');
            elseif bg_polynomial_degree == 0
                title('Select background points for CONSTANT background function');
            else
                title(['Select background points for polynomial of degree ' num2str(bg_polynomial_degree) ]);
            end      
        
            hFH = impoly();

            if ~isempty(hFH)
                background_points = hFH.getPosition;

                if size(background_points,1) > bg_polynomial_degree
                    p = polyfit(background_points(:,1), background_points(:,2),bg_polynomial_degree);

                    bg = zeros(size(t_fluo));

                    for ii = 1:bg_polynomial_degree + 1
                        bg = bg + p(ii)*t_fluo.^(bg_polynomial_degree - (ii-1)); 
                    end
                else
                    hFH = [];
                end
                
            end
        
        elseif strcmpi(background_function,'spline') || strcmpi(background_function,'splines')
            % 
            title('Select background points for SPLINE background function');
            hFH = impoly();
            background_points = hFH.getPosition;

            bg = spline(background_points(:,1), background_points(:,2),t_fluo);
        else
            error('Unknown background function! Choose either spline or polynomial!');
        end
    end
    close(gcf);
end