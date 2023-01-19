function [normalized_im] = normalize_image(input_im, method, window_size)
% This function performs image normalization on a double grayscale image
% given the method specified by the user.
%
% INPUT 
% input_im: double grayscale image to be normalized.
% method: string matching the method to be used. 'zeroCenter' for a mean of
% zero and std of one. '0to1' for values to be spread such the the minimum
% is 0 and the maximum is 1. 'windows' each pixel is normalized with
% respect to the maximum and minimum in its local neighborhood.
% window_size: only used for the windows method. This input specifies the
% size of the square used as the local neighborhood. So 3 would use a 3x3
% neighborhood around each image
%
% OUTPUT
% normalized_im: double grayscale matrix of the normalized image.
%
% Written by R.A. Manzuk
% 11/07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % we might normalize a couple of different ways depending on what the
    % user wants, so start with if statements to get to the method of
    % choice
    if strcmp(method, 'zeroCenter')
        % in this method, we'll centering the image at zero and giving it a
        % standard deviation of 1
        % start by reshaping the image into a column vector
        reshaped_im = reshape(input_im,[],1);

        % and subtract the mean
        centered_im = reshaped_im - mean(reshaped_im);

        % divide by the standard deviation to give a standard deviation of
        % 1
        normalized_column = centered_im/std(centered_im);

        % reshape back into a column
        normalized_im = reshape(normalized_column,size(input_im,1),size(input_im,2));
     
    elseif strcmp(method, '0to1')
        % in this method, we'll just be spreading the image between 0 and 1
        % start by reshaping the image into a column vector
        reshaped_im = reshape(input_im,[],1);
        
        % subtract the min
        min_sub = reshaped_im - min(reshaped_im);

        % and divide by the max
        normalized_column = min_sub/max(min_sub);

        % reshape back into a column
        normalized_im = reshape(normalized_column,size(input_im,1),size(input_im,2));
    elseif strcmp(method, 'windows')
        % to not impose square features, we'll use a circular neighborhood
        % filter of radius specified by the user
        % first make that filter
        circ_filter = fspecial("disk", window_size);
        circ_filter(circ_filter > 0) = 1;
        min_filtered = ordfilt2(input_im,1,circ_filter);

        % subtract the min filtered image
        min_subbed = input_im - min_filtered;

        % max filter the same windows
        max_filtered = ordfilt2(min_subbed,sum(circ_filter,'all'),circ_filter);

        % if the max filter is 0 for a pixel, just replace it with a very
        % small amount
        max_filtered(max_filtered == 0) = 1e-14;

        % and divide by the max
        normalized_im = min_subbed./max_filtered;
    else 
        error('not a valid normalization method string')
    end
    
end






