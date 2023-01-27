function save_superpixels(in_dirs,out_dir,im_ext, n_superpixels)
% This function takes a set of three channel images (input as separate
% paths per channel) and saves .tif images in the output directory. The
% stored tif images are label matrices for the superpixel oversegmentation
% of the image.
%
%
% INPUT
% in_dirs: cell array containing the string full paths to the directories for each
% channel to be considered when making superpixels.
% out_dir: string full path where the superpixel label matrices should be
% saved.
% im_ext: extention for the input images. 
% n_superpixels: number of superpixels desired in the oversegmentation
%
% OUTPUT
% no outputs, but the superpixels will go in the out_dir specified.
% Written by R.A. Manzuk
%
% Thursday, January 26, 2023 at 5:03:35 PM
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
    
    % Make a waitbar to keep track of how long it's taking
    f = waitbar(0, 'Making superpixels');

    % now we can just loop through and handle each image
    for i = 1:num_ims
        
        % sub loop to read in each channel
        this_im = zeros(im_height,im_width,numel(in_dirs));
        for j = 1:numel(in_dirs)
            this_dir = dir(fullfile(in_dirs{j},im_ext));
            this_dir(strncmp({this_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
            this_im(:,:,j) = im2double(imread(fullfile(this_dir(i).folder, this_dir(i).name)));
        end

        % make those superpixels
        label_mat = superpixels(this_im(:,:,1:3), n_superpixels);
        
        % save the image with the name of the input 
        imwrite(uint16(label_mat),fullfile(out_dir,test_dir(i).name));

        % update the waitbar
        waitbar(i/num_ims, f, 'Making superpixels')
    end

    % close the waitbar
    close(f);
end