function [FaceColor, EdgeColor] = generate_cmh_colors()
%
%
% WT, blmp-1;elt-3  colors for all fluorescent markers are defined here
%
% define_marker_colors.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%
    FaceColor(1,:) = [255,160,110]/255; % dpy-6 WT color
    EdgeColor(1,:) = [252,141,89]/255; % dpy-6 WT color

    FaceColor(2,:) = [240,160,150]/255; % mlt-10 WT color
    EdgeColor(2,:) = [220,48,50]/255; % mlt-10 WT color

    FaceColor(3,:) = [255,120,80]/255; % let-7 WT color
    EdgeColor(3,:) = [239,101,72]/255; % let-7 WT color    

    FaceColor(4,:) = [240,160,150]/255; % zk180.5 WT color
    EdgeColor(4,:) = [227,50,50]/255; % zk180.5 WT color
    
    %%%%%
    FaceColor(5,:) = [110,150,240]/255; % dpy-6-double color
    EdgeColor(5,:) = [60,80,225]/255; % dpy-6-double color

    FaceColor(6,:) = [180,240,250]/255; % mlt-10 double color
    EdgeColor(6,:) = [5,112,176]/255; % mlt-10 double color

    FaceColor(7,:) = [130,190,247]/255; % let-7 double color
    EdgeColor(7,:) = [60,70,166]/255; % let-7 double color

    FaceColor(8,:) = [180,240,250]/255; % zk180.5 double color
    EdgeColor(8,:) = [5,112,176]/255;  % zk180.5 double color


end