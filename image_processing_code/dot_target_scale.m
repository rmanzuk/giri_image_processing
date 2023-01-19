function [per_pixel_scale1,per_pixel_scale2] = dot_target_scale(target_im_path)
% This function takes in an image of the dot target used as a scale and
% distortion standard for GIRI, and gives back the size of the pixels in
% um, using two different methods. This function uses the known distance
% between dot centers: 250 um.
%
% IN:
%
% target_im_path: string with the full path to the dot target image to be
% used. Should be just a single channel, grayscale image (of any reasonable
% file type). 
%
% OUT:
%
% per_pixel_scale1: size of pizels in the image in um. This one is
% calculated by taking distance between each dot and its 4 neighbors above,
% below, left, and right and is given as the mean of all of those
% measurements over the entire image.
%
% per_pixel_scale2: size of pizels in the image in um. For this one, we
% figure out the headings that the rows of dots follow in the image and
% takes an improfile along each row of dots. From those improfiles we can
% then find the numer of dots in each row and the length of the row to get
% the average number of pixels between dots.
%
% Created: 09/15/2022, R. A. Manzuk
    %% begin the function
    % read in the image 
    target_im = im2double(imread(target_im_path));

    % we'll assume we don't know the scale of the image at all, so we'll
    % take a small portion of the image from the middle to test different
    % blob sizes for the dots. 

    % extract the center of the image. 400x400 pixels should be large
    % enough.
    image_center = round(size(target_im)/2);
    test_im = target_im((image_center(1)-200):(image_center(1)+199), (image_center(2)-200):(image_center(2)+199));
    
    % and iterate through a set of blob radii to see the responses at many
    % sizes
    radii = [1:20];
    mean_responses = zeros(numel(radii),1);
    for i = 1:numel(radii)
        these_blobs = detect_blobs(test_im,radii(i),3,0.05);
        mean_responses(i) = mean(these_blobs(:,3));
    end

    % fingers crossed the radius with the maximum mean response is the best
    % blob size to use for actual blob detection
    [~,ind] = max(mean_responses,[],'omitnan');
    blob_radius = radii(ind);

    % now detect blobs for real
    dot_blobs = detect_blobs(target_im, blob_radius, 3, 0.05);

    % one way to measure the pixel scale is to trace along rows/columns
    % of dots and average the length of the whole row/column given the
    % number of dots
    % first we need to know the relative orientation of rows/columns in the
    % image.
    % get the 5 nearest neighbors for each dot (nn1 will be the dot being
    % matched with itself. The other 4 should be up, down, left, right)
    [nn_inds,nn_dists] = knnsearch(dot_blobs(:,1:2),dot_blobs(:,1:2),'K',5);

    % throw out the first columns
    nn_inds = nn_inds(:,2:end);
    nn_dists = nn_dists(:,2:end);

    % okay, so the modal nn_dist is probably a good generalization for the
    % number of pixels between dots, but ones near the mode are also
    % probably valid and important, so lets get all the distances near the
    % mode
    modal_dist = mode(nn_dists(:));
    in_range_inds = nn_dists <= (modal_dist + 1.5) & nn_dists >= (modal_dist - 1.5);
    all_nn_dists = nn_dists(in_range_inds);
    
    % and their mean is probably one okay measure of the pixel-wise distance
    % between all dot centers
    mean_dot_dist = mean(all_nn_dists);

    % the dots should be 250um apart
    per_pixel_scale1 = 250/mean_dot_dist;
    % going back to tracing along rows and columns, we can use
    % in_range_inds to filter the nn_inds and figure out the general
    % vector headings for rows/columns in the image.
    % iterate throught the columns of nn_dists
    headings = [];
    for i = 1:size(nn_dists,2)
        % use written function to calculate vector headings between points
        headings = [headings ; pwise_headings(dot_blobs(in_range_inds(:,i),1:2), dot_blobs(nn_inds(in_range_inds(:,i),i), 1:2))];
    end


    % and we'll want to group anything near 2pi (360 degrees) with those
    % near 0, so adjust
    near_2pi_ind = headings > 11*pi/6;
    headings(near_2pi_ind) = headings(near_2pi_ind) - (2 *pi);
    
    % in order to properly take the average headings, we'll need to group
    % those that represent the same general heading but aren't equal
    % start with the unique headings
    unique_headings = unique(headings);

    % sort first to make later steps easier
    sorted_unique = sort(unique_headings);
   
    % and take the pairwise distances between all of those
    heading_dists = pdist2(sorted_unique,sorted_unique);

    % we can assume that the headings that should be grouped together
    % should all be within 0.5 radians (probably even closer, but doesn't
    % matter).
    close_enough_inds = heading_dists < 0.5;

    % and we should only have 4 true headings, so we can find those
    heading_group_flags = zeros(numel(sorted_unique),1);
    for i = 1:4
        if i == 1
            heading_group_flags(close_enough_inds(:,1) == 1) = 1;
        else
            heading_group_flags(close_enough_inds(:,next_start_col) == 1) = i;
        end
        next_start_col = find(heading_group_flags == 0,1,'first');
    end

    % now we can go through and take the mean of the 4 headings from teh
    % headings vector.
    final_headings = zeros(4,1);
    for i = 1:4
        % take the headings that correspond to the current one
        current_headings = ismember(headings,sorted_unique(heading_group_flags==i));
        % and take the mean
        final_headings(i) = mean(headings(current_headings));
    end

    % and now apply those headings to put down some line segements and
    % measure the distances between a bunch of dots

    % first, to go along rows of dots, we need a set of the leftmost dot
    % centers. We dont need to have every single one, but we probably can
    % get most of them if we divide the image up vertically by the distance
    % between dots that we calculated above and look for the leftmost in
    % each window
    % first, adjust the dot distance by the known skew from the headings
    cosines = cos(final_headings);
    % the cosine we want will be close to 1
    cos_dists = 1-cosines;
    [~,cos_ind] = min(cos_dists);
    % and apply it to the dot distance
    window_size = mean_dot_dist * cosines(cos_ind);
    % and get the windows
    windows = [0:window_size:size(target_im,1)];
    % bin the points
    bin_ids = discretize(dot_blobs(:,1),windows);
    % go through and pick out the leftmost in each bin
    leftmost_inds = zeros(max(bin_ids),1);
    for i = 1:max(bin_ids)
        in_this_bin = bin_ids == i;
        % just grab the col coords of the dots
        col_coords = dot_blobs(:,2);
        % set those outside of the window to nan
        col_coords(~in_this_bin) = nan;
        % and get the ind
        [~,leftmost_inds(i)] = min(col_coords);
    end

    % and go through and improfile the images along all these segments and
    % get the distances between many dot centers

    % we'll use the tangent of the heading to figure out where we're going
    tangents = tan(final_headings);
    % we'll be using the one nearest 0, and make sure it's positive
    this_tangent = min(abs(tangents));
    
    dists_betw_dots = zeros(numel(leftmost_inds),1);
    for i = 1:numel(leftmost_inds)
        start_point = dot_blobs(leftmost_inds(i),1:2);
        % figure out how much room we have to work with
        cols_available = size(target_im,2)-start_point(2);
        
        % use the tangent to figure out how far down we need to go for the
        % end point's row coordinate
        rows_below = this_tangent * cols_available;
        % grab the im profile
        dot_profile = improfile(target_im,[start_point(2),start_point(2)+cols_available], [start_point(1),start_point(1)+rows_below]);

        % and the local minima should be the dot centers (we can use half
        % the already known distance between dots to avoid any funny
        % business
        minima = islocalmin(movmean(dot_profile,mean_dot_dist/4),'MinSeparation',mean_dot_dist/2,'MinProminence',0.5,'FlatSelection','center');

        % take the first and last min
        first_min = find(minima, true, 'first');
        last_min = find(minima, true, 'last');

        % the average distance is the distance between the first and last
        % divided by the number of minima-1
        dists_betw_dots(i) = (last_min-first_min)/(sum(minima) - 1);
    end
    % get the pixel scale with the mean, trimming a few outliers
    per_pixel_scale2 = 250/trimmean(dists_betw_dots,10);
end
