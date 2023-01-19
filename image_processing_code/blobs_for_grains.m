function [positions, rads, response_values] = blobs_for_grains(img,radii,local_thresh,global_thresh, plot_flag)
% This function uses blob detection at multiple scales input by the user to
% detect grains, provided they can be represented as blobs.
%
% IN:
% img: 2d grayscale image in double format.
% radii: vector of desired radii (in number of pixels) to test for blob
% detection
% local_thresh: threshold for local nonmaximum supression. 3 tends to be a
% good value.
% global_thresh: global response threshold, below which a response is not
% considered to be a blob. 0.05 may be a good starting value.
% plot_flag: logical if you would automatically like a figure with the
% image and detected blobs at the end.
%
% OUT: 
% positions: 2 column matrix with the row and column positions of the blob
% centers on the full resolution images.
% rads: column vector of same number of rows as positions with the radius
% of each blob given in positions. 
% resoponse_values: column vector of same number of rows as positions with
% the filter response value of each blob given in positoins.
% Ryan A. Manzuk 11/03/2022
    %%
    % to be computationally efficient, we will scale the image to match the
    % relative sizes of the radii requested, not scale the filter.
    % so start with a basic filter radius.
    base_rad = 2*sqrt(2);

    % and the relative image scalings to test for all of these radii
    im_scalings = base_rad./radii;
    
    % set up empty stuff to hold results
    positions = [];
    rads = [];
    response_values = [];

    % now we're free to iterate through the set of radii, resize the image,
    % and detect blobs
    % set up a variable for all responses
    test_size = imresize(img,max(im_scalings));
    response_space = zeros(size(test_size,1),size(test_size,2),numel(radii));
    for i = 1:numel(im_scalings)

        % first resize
        im_resized = imresize(img, im_scalings(i),'bicubic');
        % need a sigma for our LoG from the radius
        sigma = base_rad/sqrt(2);
        %define and implement Laplacian of Gaussian
        filter_size = 2*ceil(3*sigma)+1;
        LoG = sigma^2 * fspecial('log',filter_size,sigma);
        % then run the LoG over the image
        filtered_img = imfilter(im_resized,LoG,'same','replicate');
        % square the responses and scale to the size of the largest image
        % we'll use
        response_space(:,:,i) = imresize(filtered_img.^2,[size(test_size,1),size(test_size,2)],'bicubic');
    end

    % now we need to do nonmaximum supression 
    max_space = zeros(size(test_size,1),size(test_size,2),numel(radii));
    for i = 1:numel(radii)
        max_space(:,:,i) = ordfilt2(response_space(:,:,i), local_thresh^2, ones(local_thresh));
    end
    
    % designate the maxima
    for i = 1:numel(radii)
        max_space(:,:,i) = max(max_space(:,:,max(i-1,1):min(i+1,numel(radii))),[],3);
    end
    max_space = max_space .* (max_space == response_space);

    % now go into the max space and grab the blobs
    rows = [];   
    cols = []; 
    rads = [];
    response_values = [];
    for i = 1:numel(radii)
        [r,c,value] = find(max_space(:,:,i).*(max_space(:,:,i) >= global_thresh));
        numBlobs = length(r);
        if(numBlobs >0)
            rads = [rads; repmat(radii(i),numBlobs,1)];
            rows = [rows; r];
            cols = [cols; c];
            response_values = [response_values ;value];
        end
    end
    
    positions = round([rows,cols]./max(im_scalings));

    % once we're done we can go through and plot if the user wants
    if ~isempty(positions) && plot_flag
        figure
        % show the pic
        imshow(img)
        hold on
        % define a set of angles
        theta = linspace(0, 2*pi, 24);
        % and go through and draw each blob
        for i = 1:numel(rads)
            r = rads(i);
            plot(positions(i,2) + r*cos(theta), positions(i,1) + r*sin(theta), 'r-', 'linewidth',2);
        end
    end
end