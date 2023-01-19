function apply_lcc(in_dir,out_dir,im_ext,adjustments)
%Applies gray or lcc/white card image adjustments to a set of sample images
%taken at the same wavelength
%
% INPUT 
% in_dir: String with full path to the directory the input images, unadjusted
% out_dir: String with full path to where you'd like the adjusted images
% sent
% im_ext: String with the extension for the images
% adjustments: 2d double matrix of same size as the input images with the
% lcc adjustments, output from assess_lcc.m
%
% OUTPUT
% there are no outputs, but the adjusted images will be in the indicated
% output folder
% Written by Cedric Hagen, edited by R.A. Manzuk
% 10/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end

    % get all of the images into a directory structure
    input_files = dir(fullfile(in_dir, im_ext));
    input_files(strncmp({input_files.name}, '.', 1)) = []; %remove files in dir starting with '.'

    % This will check if the directory exists, if not, it will make them 
    mkdir(out_dir);

    % for writing the images later, we want the input bit depth
    img_info = imfinfo(fullfile(in_dir,input_files(1).name));
    bit_depth = img_info(1).BitDepth;

    % now we just have to loop through all of the images, read them, and
    % export them
    % Make a waitbar to keep track of how long it's taking
    f = waitbar(0, 'Correcting the images');
    
    for i = 1:numel(input_files)
        
        % read the im
        this_im = im2double(imread(fullfile(in_dir,input_files(i).name)));
        
        % apply the adjustments
        new_im = this_im./adjustments;

        % and we need to recenter the histogram by subtracting the mean
        % increase
        mean_increase = mean((new_im-this_im),"all");

        new_im = new_im - mean_increase;

        % write the new image
        if bit_depth == 16
            imwrite(uint16(new_im.*2^16),fullfile(out_dir,input_files(i).name))
        else
            imwrite(new_im,fullfile(out_dir,input_files(i).name))
        end


        % update the waitbar
        waitbar(i/numel(input_files), f, 'Still Correcting')
    end

    % close the waitbar
    close(f);
end







