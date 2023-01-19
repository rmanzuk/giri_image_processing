function [texton_image] = texton_apply_clusters(image, scaling, filter_bank, centroids)
%
% IN:
% image: double grayscale image matrix for texton analysis
% scaling: fraction by which you would like the image to be resized prior
% to being filtered. Likely will want to downsample for efficiency.
% filter_bank: mxnxn_filter matrix with the set of filters you would like
% to use for feature extraction. Alternatively you can use a gabor filter
% bank object.
% n_features: The number of clusters you would like to have in the kmeans
% clustering of the filter responses for texton mapping
% centroids: n_centroids x n_filters matrix with cluster centroid locations in
% filter response space, likely the output of kmeans.
%
% OUT:
% texton_image: indexed image of equal size to the input image where each
% pixel's index gives the assigned cluster. 
%
% Written by R.A. Manzuk
% 11/14/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first, scale the image if we need to
    scaled_image = imresize(image,scaling, 'bicubic');

    % then we want to normalize within a 3x3 window for each image
    % following Jemwa, et. al.
    % use the normalize image function properly to do that
    %window_normalized = normalize_image(scaled_image,'windows',6);
    
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

    % and to normalize the responses, jemwa suggests normalizing the filter
    % responses according to F(x) = log(1+ ||F(x||)/0.03/||F(x)|| where
    % ||F(x)|| is the magnitude of each pixel's filter response
    response_magnitudes = sqrt(sum(filter_responses.^2,3));
    norm_factor = (log(1 + response_magnitudes)./0.03)./response_magnitudes;
    normalized_responses = filter_responses./norm_factor;
    
    % and now we reshape into a 2d matrix of columns for assigning to
    % clusters

    responses_matrix = reshape(normalized_responses,[], n_filters);

    % and we just have to find the nearest neighbor in the centroids for
    % each response
    centroid_idx = knnsearch(centroids,responses_matrix);

    % reshape the indices into an image for to output as a texture
    % segmentation
    texton_image = reshape(centroid_idx,[size(scaled_image,1),size(scaled_image,2)]);
end

