function [mean_im] = calc_mean_im(im_path,im_ext)
% One method for calcluating the needed adjustments for lcc requires the
% mean image. This function give the mean of a grayscale image dataset (all
% of equal size).
%
% IN:
% im_path: string with the full path to the images of the wavelength for
% which you'd like to take the mean.
% im_ext: String with the extension for the images
%
% OUT:
% mean_im: 2d double matrix of equal size to the input images with the mean
% of the input image set
%
% Written by R.A. Manzuk
% 10/18/2022
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
    f = waitbar(0, 'getting mean image');
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
        waitbar(i/numel(in_dir), f, 'getting mean image')
    end

    % close the waitbar
    close(f);
end
