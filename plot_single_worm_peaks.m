function [x,y] = plot_single_worm_peaks(axis_handle, x, peak_properties,bg, fitline_color)
%
% x ... this is actually t_fluo
% peak_properties ... array, output from fit_worm_peaks.m
% bg ... background, if background determined automatically from the fit this has 
%         either size 1x2 (linear bg) or 1x1 (constant bg)
%         if this is a longer vector, it is a background manually
%         determined from inputs to a figure
% fitline_color ... 3x1 vector of color
%
%
% by Wolfgang Keil, Institut Curie 2020, wolfgang.keil@curie.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~nargin ==5
        fitline_color = [0.2 0.7 0.2];
    end
    
    y = zeros(size(x));
    
    for ii = 1:size(peak_properties,1)
      % means this peak is present
      if ~isnan(peak_properties(ii,2))
          y = y + peak_properties(ii,3) *exp(-((x-peak_properties(ii,2))./(0.60056120439323.*peak_properties(ii,4))).^2);         
%          exp(-(x-peak_properties(ii,2)).^2*1/peak_properties(ii,4).^2);
      end
    end
        
    hold on;
    
 
    if nargin == 3 % no background specified at all
        plot(axis_handle,x,y, '-', 'linewidth', 1.5, 'Color',fitline_color);
    else
        
        if length(bg) == 2
            y = y + bg(1)*x + bg(2); % from fit, assuming linear background
        else
            y = y + bg; % bg is either a scalar (when automatically fitted as constant bg) 
                        % or a vector, the same size as y (when manually
                        % determined and subtracted
        end
        plot(axis_handle,x,y, '-', 'linewidth', 1.5, 'Color', fitline_color);
    end
    hold off;
    
    

end