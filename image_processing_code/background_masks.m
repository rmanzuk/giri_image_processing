function background_masks(in_dirs,out_dir,im_ext)
% This function takes a set of multispectral images and makes masks to
% remove the background based upon ratios between the channels. Unrealistic
% ratios between channels are taken to be signs of uneven lighting/shadows,
% so those areas are to be masked.
%
% INPUT
% in_dirs: cell array containing the string full paths to the directories for each
% channel to be considered when making the masks.
% out_dir: string full path where you would like the masks to be saved.
% im_ext: extention for the input images. 
%
% OUTPUT
% no outputs, but the masks will go in the out_dir specified.
% Written by R.A. Manzuk
%
% 11/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end

    % figure out the number of images we're working with by looking at the
    % first directory
    test_dir = dir(fullfile(in_dirs{1}, im_ext));
    test_dir(strncmp({test_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
    num_ims = numel(test_dir);

    % and test for image size
    test_im = im2double(imread(fullfile(test_dir(1).folder, test_dir(1).name)));
    im_height = size(test_im,1);
    im_width = size(test_im,2);
    
    keyboard
    % Make a waitbar to keep track of how long it's taking
    f = waitbar(0, 'Making masks');

    % now we can just loop through and handle each image
    for i = 1:num_ims
        
        % sub loop to read in each channel
        this_im = zeros(im_height,im_width,numel(in_dirs));
        for j = 1:numel(in_dirs)
            this_dir = dir(fullfile(in_dirs{j},im_ext));
            this_dir(strncmp({this_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
            this_im(:,:,j) = im2double(imread(fullfile(this_dir(i).folder, this_dir(i).name)));
        end

        % superpixels will help us make a nice, continuous mask
        label_mat = superpixels(this_im(:,:,1:3), 10000);

        % for stats of the superpixels, easiest to work in coulmns with the
        % image and lablematrix
        im_reshaped = reshape(this_im,[],numel(in_dirs));
        label_col = reshape(label_mat,[],1);
        
        % take the mean color values of all superpixels
        mean_colors = splitapply(@mean,im_reshaped,label_col);

        % expand out to the size of the whole image again
        mean_colors_expanded = mean_colors(label_col,:);
        mean_color_image = reshape(mean_colors_expanded,size(this_im));

        % we want to get the ratios of all columns, so we'll use n choose k
        % to get all pairings
        pairings = nchoosek([1:numel(in_dirs)], 2);

        % grab all ratios indicated by the pairings
        ratios = mean_colors_expanded(:,pairings(:,1))./mean_colors_expanded(:,pairings(:,2));

        % define an initial mask where all of the ratios look pretty good
        initial_mask = reshape(all(ratios > 0.5 & ratios < 1.8,2), size(this_im,1),size(this_im,2));

        % and decide which ratios are too big or small: 
        too_small_big = ratios < 0.5 | ratios > 2;
        bad_ratio = any(too_small_big,2);

        % and reshape into a mask
        initial_mask = ~reshape(bad_ratio,size(test_im,1),size(test_im,2));

        % open and close to clean it up we'll use the number of pixels to
        % approximate island and hole sizes
        mask2 = bwareaopen(initial_mask,round(numel(test_im)/100),4);
        mask3 = ~bwareaopen(~mask2,round(numel(test_im)/200),4);
        % write the mask with the name of the image in the test_dir
        imwrite(mask3,fullfile(out_dir,test_dir(i).name));


        % update the waitbar
        waitbar(i/num_ims, f, 'Making_masks')
    end

    % close the waitbar
    close(f);
end







