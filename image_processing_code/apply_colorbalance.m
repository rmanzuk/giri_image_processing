function  apply_colorbalance(gray_vals,wavelengths, bit_depth, in_dir,out_dir,im_ext)
% Assesses colorbalance of all images in a set based upon gray boxes in the
% input color card images. For this function, be sure to put the
% wavelenghts of each image in the image names.
%
% INPUT 
% gray_vals: number wavelengths x number gray boxes matrix with the average
% value within that box for each channel output from assess_colorbalance.m
% wavelenghts: cell array with the character vectors with the wavelengts of
% each channel as gleaned from the file names.
% bit_depth = bit depth of the images you are working with.
% in_dir: String with full path to the directory the input images, unadjusted
% out_dir: String with full path to where you'd like the adjusted images
% sent
% im_ext: String with the extension for the images
%
% OUTPUT
% nothing but the color corrected will be placed in the proper folders
% 
% Written by R.A. Manzuk
% 10/10/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end

    % okay, first, we need to figure out the gamma-corrected mapping
    % and for some future housekeeping, we need a set of values spanning
    % the whole bit space and linearly spaced indices for the gray values
    % covering that space
    bit_steps = linspace(0,1,2^bit_depth);
    gray_scaled = linspace(1/2^bit_depth,1,size(gray_vals,2));
    

    % start also with a directory of all the image folders
    whole_in_dir = dir(in_dir);
    whole_in_dir(strncmp({whole_in_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
    
    % we need to go through every channel and fit a 2nd order polynomial to
    % all of the gray values this will give us an idea of how to map all of
    % the images to the color space indicated by the color card.
    params = zeros(size(gray_vals,1),3);

    for i = 1:size(gray_vals,1)
        params(i,:) = polyfit(gray_scaled,gray_vals(i,:),2);
    end

    % and we'll make the mean gray value we're going for from the luminance
    % of rgb
    red_ind = contains(wavelengths,'625');
    green_ind = contains(wavelengths,'530');
    blue_ind = contains(wavelengths,'470');
    just_rgb = cat(3,gray_vals(red_ind,:),gray_vals(green_ind,:),gray_vals(blue_ind,:));
    gray_params = polyfit(gray_scaled,rgb2gray(just_rgb),2);

    % let's plot up all the fits and let the user select the curve to serve
    % as the reference. 
%     figure
%     hold on
%     for i = 1:size(gray_vals,1)
%         y = polyval(params(i,:),bit_steps);
%         plot(bit_steps,y,'DisplayName',[num2str(i) '. ' wavelengths{i}])
%     end
%     drawnow
%     legend
%     prompt = 'Please enter the index of the wavelength you would like to be the reference.\n';
%     wave_ind = input(prompt);
    
    % use that wavelength to set up the reference curve
    ref_poly = polyval(gray_params,bit_steps);

    if bit_depth == 16
        ref_int = uint16(ref_poly.*2^16);
    elseif bit_depth == 8
        ref_int = uint8(ref_poly.*2^8);
    else 
        disp('sorry, I do not understand that bit depth. Try again.')
        return
    end
        
    f = waitbar(0, 'Correcting the images');

    % now that we've fit the color card gray values, we can go in and see
    % how much these curves deviate from the best channel selected by the user.
    for i = 1:size(params,1)
        % first evaluate our polynomial so we can get the difference
        % between our channel's current span of color and perfect gamma
        % corrected values
        this_poly = polyval(params(i,:),bit_steps);

        % and we want to make sure it's scrunched between 0 and 1
        if min(this_poly) < 0
            this_poly = this_poly - min(this_poly);
        end

        if max(this_poly) > 1
            this_poly = this_poly./max(this_poly);
        end

        % and we can use this polygon evaluation to map into gamma
        % corrected space wiht integers.
        if bit_depth == 16
            poly_int = uint16(this_poly.*2^16);
        elseif bit_depth == 8
            poly_int = uint8(this_poly.*2^8);
        end

        % and make a multiplication key between this curve and the
        % reference curve
        mult_key = double(ref_int)./double(poly_int);

        % we're now ready to correct the images
        % start by selecting the proper folder in the input directory and
        % making an image directory
        this_folder = contains({whole_in_dir(:).name},wavelengths(i));
        image_dir = dir(fullfile(in_dir,whole_in_dir(this_folder).name,im_ext));
        image_dir(strncmp({image_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
        keyboard
        % make sure we have an output_directory
        mkdir([out_dir '/' wavelengths{i} 'nm'])
        % now we can loop through and correct the images
        for j = 1:numel(image_dir)
            this_im = imread(fullfile(image_dir(j).folder,image_dir(j).name));
            % we'll do this by setting up a multiplier image
            mult_im = mult_key(this_im(:)+1);
            mult_im = reshape(mult_im,size(this_im,1),size(this_im,2));
            
            % and all we need to do is subtract the difference image to be
            % ready for export
            if bit_depth == 16
                for_export = im2double(this_im).*mult_im;
                imwrite(uint16(for_export.*2^16),fullfile([out_dir '/' wavelengths{i} 'nm'],image_dir(j).name))
            elseif bit_depth == 8
                for_export = im2double(this_im).*mult_im;
                imwrite(uint8(for_export.*2^8),fullfile([out_dir '/' wavelengths{i} 'nm'],image_dir(j).name))
            end

            % update the waitbar
            total_done = (i-1)*numel(image_dir) + j;
            to_do = numel(wavelengths) * numel(image_dir);
            waitbar(total_done/to_do, f, 'Correcting color')
        end
    end
    % close the waitbar
    close(f);
end