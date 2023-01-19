function [adjustments] = assess_lcc2(im_path,im_ext)
%Calculates needed adjustments for uniform image intensity despite
%directionality of light source. This version (different from assess_lcc.m)
%does not use a white or grey card and in stead uses the mean image of the
%whole data set for that wavelength.
%
% IN:
% im_path: string with the full path to the images of the wavelength for
% which you'd like to assess the lcc.
% im_ext: String with the extension for the images
%
% OUT:
% adjustments: 2d double matrix of equal size to the input image with the needed
% adjustments to achieve even lighting
%
% Written by R.A. Manzuk
% 10/17/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end

    % make a directory for reading files
    in_dir = dir(fullfile(im_path,im_ext));
    in_dir(strncmp({in_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'

    % read in the first image so that we know the size for the mean image
    test_im = im2double(imread(fullfile(im_path,in_dir(1).name)));

    % set up a mean image that is just 0s for now
    mean_im = zeros(size(test_im,1), size(test_im,2));
    
    % taking the mean of all images could take a bit, so set up a waitbar
    f = waitbar(0, 'assessing lighting field');
    % now we can loop through the whole channel image set and update the
    % mean. We use a weighted mean so we don't have to upload all images
    % at once, causing memory issues.
    for i = 1:numel(in_dir)
        this_im = im2double(imread(fullfile(im_path,in_dir(i).name)));
        % stretch the image so it contributes evenly to the mean.
        stretched_im = imadjust(this_im);
        % calculate the weight base upon which number image this is
        weight = 1/i;
        % take the weighted mean
        mean_im = mean_im.*(1-weight) + stretched_im.*weight;

        % update the waitbar
        waitbar(i/numel(in_dir), f, 'assessing lighting field')
    end
    
    % just in case there is some smallscale features in the mean image,
    % let's just blur it a bit so we don't have to worry about those
    mean_im2 = imgaussfilt(mean_im,size(mean_im,1)/10);
    
    % and finally, the adjustments are just that blurred mean image, scaled
    % up such that the max is 1.
    adjustments = mean_im2 + (1-max(mean_im2(:)));

    % close the waitbar
    close(f);
end
