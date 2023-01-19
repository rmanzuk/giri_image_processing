function [blobs] = detect_blobs(img,radius,local_thresh,global_thresh)
% This function takes a grayscale double image and detects the blobs of a
% given radius using a Laplacian of a Gaussian and non-maximum suppression.
%
% IN:
%
% img: m x n image matrix (grayscale and double).
%
% radius: radius, in pixels, of the blobs to be detected. I find that 16
% works well for the distortion target dots.
%
% local_thresh: threshold for non-maximum supression. I find that a value
% of 3 works well for most images.
%
% global_thresh: threshold above which blobs that make it past non-maximum
% supression will be considered blobs in the image. I find that a value of
% 0.05 works well to the distortion target dot.
%
% OUT:
%
% blobs: n_blobs x 3 matrix. The first two columns are the row and column
% indices (in image coords) for the locations of the blob centers. The
% third column is the response value of the LoG for that blob.
%

% Ryan A. Manzuk 01/05/2021
    %%
    sigma = radius/sqrt(2);
    %define and implement Laplacian of Gaussian
    filter_size = 2*ceil(3*sigma)+1;
    LoG = sigma^2 * fspecial('log',filter_size,sigma);
    
    % then run the LoG over the image
    filtered_img = imfilter(img,LoG,'same','replicate');
    % square the responses
    filtered_img = filtered_img.^2;
    % resize the responses back to the image size.
    response_space = imresize(filtered_img,size(img),'bicubic');
    
    %nonmaximum suppression
    max_space = ordfilt2(response_space, local_thresh^2, ones(local_thresh));
    
    max_space = max_space .* (max_space == response_space);
    %Define blobs
    [rows,cols,values] = find(max_space.*(max_space >= global_thresh));
    
    blobs = [rows,cols,values];

end