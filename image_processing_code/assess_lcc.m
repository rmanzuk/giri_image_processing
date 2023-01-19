function [adjustments] = assess_lcc(lcc_im_path,plot_flag)
%Calculates needed adjustments for 'gray/white' card images to have uniform
%image intensity despite directionality of light source.
%
% IN:
% lcc_im_path: string with the full path to the image of the lcc card or
% gray card
% plot_flag: logical indicating if you would like a visual of the
% adjustments or not
%
% OUT:
% adjustments: 2d double matrix of equal size to the input image with the needed
% adjustments to achieve even lighting
%
% Written by Cedric Hagen, edited by R.A. Manzuk
% 10/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in the image
lcc_im = im2double(imread(lcc_im_path));

    if plot_flag
        %Plotting switch
        figure
        subplot(3, 1, 1);
        imshow(lcc_im);
        colorbar;
        axis('on', 'image');
        title('Original Image');
        subplot(3, 1, 2);
        % Display image with that colormap.
        imagesc(lcc_im) 
        colorbar;
        axis('on', 'image');
        title('Colormapped Image');
        %Plots the original 'gray' image and a colormapped version that highlights
        %differences better
    end

    % just in case there is some smallscale features in the lcc image,
    % let's just blur it a bit so we don't have to worry about those
    blurred_lcc = imgaussfilt(lcc_im,size(lcc_im,1)/20);
    
    % and see what the mean intensity is
    avg_int = mean(blurred_lcc(:));

    % the adjustments are just the local deviations from the
    % mean 
    adjustments = blurred_lcc/avg_int;
    
    % if there's a plotting flag, we can show what those adjustments look
    % like
    if plot_flag
        %Plotting switch
        new_I=blurred_lcc./adjustments;
        %Applies adjustments to the 'gray' card image
        
        subplot(3,1,3)
        imagesc(new_I);
        colorbar;
        axis('on', 'image');
        title('Adjusted Image');
        %Plots the adjusted image with uniform intensity
    end

end
