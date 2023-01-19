function [normalized_responses] = texton_jemwa(image, scaling, filter_bank)
% This function uses a prescribed filter bank to perform texton analysis of
% an image using pre and post-processing steps described by Jemwa, et. al.,
% 2021 (1).
%
% Jemwa, Gorden T., and Chris Aldrich. "Estimating size fraction categories
% of coal particles on conveyor belts using image texture modeling
% methods." Expert Systems with Applications 39.9 (2012): 7947-7960.
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
% dim_red: dimensionality of final output you would like, after reduction.
% So if you would like a 50x50xn_filter final output, input 50.
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
% 11/07/2022
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

    % for export to build the library, we can downsample with a box filter
    % to just a 100x100 image
    %responses_downsampled = imresize(normalized_responses,[dim_red,dim_red], "box");
    
    % and we'll spit out a matrix of these responses
    %responses_matrix = reshape(responses_downsampled,[], n_filters);


    % next step is to do a neighborhood normalization of the responses,
    % just like we did with the input image
    %responses_window_normalized = zeros(size(filter_responses,1),size(filter_responses,2),size(filter_responses,3));
    %for i = 1:size(responses_window_normalized,3)
        %responses_window_normalized(:,:,i) = normalize_image(filter_responses(:,:,i),'windows',6);
    %end

    % Jemwa says to apply an activation and thresholding nonlinearity to
    % get values of the responses between 0 and 1, so we'll use a sigmoid
    % to accomplish that. Extra lines are to create and extract data from a
    % deep learning array, which is the data format needed for the sigmoid
    % function.
    %responses_window_normalized = dlarray(responses_window_normalized,'SSCB');
    %activated_responses = sigmoid(responses_window_normalized);
    %activated_responses = extractdata(activated_responses);

    % now we want to use the l1 norm of given neighborhoods across the filter
    % space to normalize.
    % the l1 norm of 3x3xn_filter neighborhoods is just the convolution
    % of the filter responses with a 1s filter of that size
    %neighborhood_size = 10;
    %local_norms = zeros(size(activated_responses,1),size(activated_responses,2),size(activated_responses,3));
    %for i = 1:size(local_norms,3)
        % first subtract the mean
        %demeaned = activated_responses(:,:,i) - mean(activated_responses(:,:,i),"all");

        % and convolve absolute values with a ones spatial filter to get
        % the local l1 norm
        %local_norms(:,:,i) = conv2(abs(demeaned),ones(neighborhood_size,neighborhood_size),'same');
    %end
    
    % to get a single set of norms for each pixel, take the 3rd dimension
    % sum
    %l1_norms = sum(local_norms,3);

    % and if we want those norms to equal 1, simply divide by the norms
    %normalized_responses = activated_responses./l1_norms;

    % and we'll want to k means these responses, first put them in a 2d matrix    
    %responses_matrix = reshape(normalized_responses,[], n_filters);

    % do the k means with the input number of clusters
    %[~,centroids] = kmeans(responses_matrix,n_features,'MaxIter',1000);
    
    % and the texton map or segmentation is just those clustering results
    % reshaped into the image
    %texton_map = reshape(k_inds,[size(scaled_image,1),size(scaled_image,2)]);

    % for export, and reduced dimensionality use a boxcar filter to get it
    % down to the final proper dmension. 
    %boxcar_size = size(texton_map)-dim_red+1;
    %box_kernel = ones(boxcar_size(1),boxcar_size(2),n_filters);

    % and convolve to get the downsampled set
    %reduced_responses = convn(normalized_responses,box_kernel);
    
end

