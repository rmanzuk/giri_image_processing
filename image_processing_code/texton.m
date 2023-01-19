function [texton_map, centroids, filter_responses] = texton(image, scaling, filter_bank, n_features)
% This function uses a prescribed filter bank to perform texton analysis of
% an image with no pre or post processing
%
% IN:
% image: double grayscale image matrix for texton analysis
% scaling: fraction by which you would like the image to be resized prior
% to being filtered. Likely will want to downsample for efficiency.
% filter_bank: mxnxn_filter matrix with the set of filters you would like
% to use for feature extraction. Alternatively you can use a gabor filter
% bank object.
% n_features. The number of clusters you would like to have in the kmeans
% clustering of the filter responses for texton mapping
%
% OUT:
% texton_map: indexed image of equal size to the input image where each
% pixel's index gives the assigned kmeans category.
% centroids: matrix with row vectors for the centroid locations of the k
% means clusters.
% filter responses: matrix with row vectors for every pixel's responses to
% all filters 
%
% Written by R.A. Manzuk
% 10/27/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first, scale the image if we need to
    scaled_image = imresize(image,scaling, 'bicubic');
    
    % and we can go ahead and do the filtering depending on if this is a
    % gabor bank or not
    if isa(filter_bank,'gabor')
        % if it's a gabor bank, we can just apply directy with the matlab
        % function
        filter_responses = imgaborfilt(scaled_image,filter_bank);
        % we'll want to store the number of filters for later
        n_filters = numel(filter_bank);
    else
        n_filters = size(filter_bank,3);
        % set up an empty recepticle for the filter responses
        filter_responses = zeros(size(scaled_image,1),size(scaled_image,2),n_filters);
        % go through and get those responses
        for i = 1:size(filter_bank,3)
            filter_responses(:,:,i) = conv2(scaled_image,filter_bank(:,:,i),'same');
        end
    end

    % and we'll want to k means these responses, first put them in a 2d matrix    
    responses_matrix = reshape(filter_responses,[size(scaled_image,1)*size(scaled_image,2), n_filters]);

    % do the k means with the input number of clusters
    [k_inds,centroids] = kmeans(responses_matrix,n_features,'MaxIter',1000);
    
    % and the texton map or segmentation is just those clustering results
    % reshaped into the image
    texton_map = reshape(k_inds,[size(scaled_image,1),size(scaled_image,2)]);
end
