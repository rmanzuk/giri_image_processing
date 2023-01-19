function [] = response_overlay(grayscale_image, responses, opacity)
% This function overlays filter responses or other scaled values that you
% would like to visualized on top of an image.
%
% INPUT
% grayscale_image: 2D matrix with a grayscale image.
% responses: 2D matrix of equal size to grayscale_image with the responses
% or other variable you would like to overlay on the image.
% opacity: opacity of the overlay, should be between 0 and 1.
%
% OUTPUT
% no outputs, but the figure will pop up.
%
% Written by R.A. Manzuk
% 11/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % expand the grayscale image to be 3 channels
    grayscale_rep = repmat(grayscale_image,[1 1 3]);
    % show it
    imshow(grayscale_rep); hold on
    % and overlay the scaled responses
    overlay = imagesc(responses);
    % scale the opacity
    overlay.AlphaData = ones(size(grayscale_image,1),size(grayscale_image,2))*opacity;
    hold off
end

