function superpixel_trainer(obj, n_supers)
% This function takes in a petro_image object and a specified number of
% superipixels, and runs through some user input protocols to allow the
% user to specify the classes present in the image (from
% master_class_list.mat) and then click example superpixels for each class.
% The resulting clicked superpixels are stored with the identified class
% name for each index of identified superpixels in the class_labels
% property of the object. The trainer also updates the class_regions with
% the outlines of regions specified as each class for making training masks
% later. If the superpixel label matrix has not already been created and
% saved, this function creates and saves in the specified path.
%
% IN:
% obj: instance of the petro_image object class
% n_supers: single value of the number of superpixels in the
% oversegmentation to be used for clicking
%
% OUT:
%
% none, but the object from the input is updated
%
% R. A. Manzuk 
% written: Friday, February 17, 2023 at 3:08:23 PM
% R. A. Manzuk 11/30/2021
% Major edit: Tuesday, February 7, 2023 at 11:10:29 AM
    %% begin the function
    % start by reading in red green and blue images, and we'll actually
    % just make them 8 bit versions because be don't need all that color
    % depth for now
    disp('reading in red, green, and blue images');
    red_im = uint8(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 625}, [obj.sample_name, obj.default_ext]))/256);
    green_im = uint8(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 530}, [obj.sample_name, obj.default_ext]))/256);
    blue_im = uint8(imread(fullfile(obj.main_path, obj.im_subpaths{obj.wavelengths == 470}, [obj.sample_name, obj.default_ext]))/256);
    
    % concatenate into rgb and clear individual images
    rgb_im = im2double(cat(3,red_im,green_im,blue_im));
    clear red_im green_im blue_im

    % now we have to check if there are superpixels already made
    superpix_fname = fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers), [obj.sample_name, obj.default_ext]);
    if isfile(superpix_fname)
        disp('loading in superpixel indices');
        % if they are there we can just load them
        label_mat = imread(superpix_fname);
    else
        label_mat = obj.get_superpix(n_supers, rgb_im);
    end
    
    % hard coding in the scaling here because all images are the same.
    % Could bring outside function for future*******
    scaling = 0.25;
    
    % we'll need a superpixel boundary mask for later
    bound_mask = boundarymask(imresize(label_mat,scaling,'nearest'));

    % because we scale, the boundary mask is a bit chunky, so we can erode
    % it with  a small disk
    erode_element = strel('square',2);
    bound_mask = imerode(bound_mask,erode_element);

    % scale the rgb image as well
    rgb_im = imresize(rgb_im,scaling,"bicubic");
    
    % now we want to assess classes.
    % start by displaying the image
    just_im = figure();
    imshow(rgb_im)

    % and load the master class list
    load("master_class_list.mat");

    % if no classes currently present, we neet to ask the user to give
    % classes
    if isempty(obj.classes_present)
        fprintf('no classes currently entered for this image. \n');
        fprintf('please enter a list of indices for the classes present in this image. \n');
        fprintf('possible classes are: \n');
        for i = 1:numel(master_class_list)
            fprintf([num2str(i) '. ' master_class_list{i} '\n']);
        end
        % user can now input which classes they see
        class_inds = input('enter indices as a vector \n');
        % and make that this instance's class list 
        obj.classes_present = master_class_list(class_inds);
    else
        % otherwise we can display the current class list and ask if the
        % user would like to add any
        fprintf('Current class list is : \n');
        for i = 1:numel(obj.classes_present)
            fprintf([obj.classes_present{i} '\n']);
        end
        user_desire = input('Would you like to add any more [y/n] \n','s');

        if strcmpi(user_desire,'y')
            fprintf('please enter a list of indices for the classes present in this image. \n');
            fprintf('possible classes are: \n');
            for i = 1:numel(master_class_list)
                fprintf([num2str(i) '. ' master_class_list{i} '\n']);
            end
             % user can now input which classes they see
            class_inds = input('enter indices as a vector \n');
            % and add to this instance's class list 
            obj.classes_present = [obj.classes_present, master_class_list(class_inds)];
        end
    end

    close(just_im)

    % okay, ready to move on! might as well set the number of classes 
    n_classes = numel(obj.classes_present);
    
    % make empty string array to hold class assignments for each superpixel
    % only if there isn't one already
    if ismember(n_supers, obj.n_superpixels)
        n_superpix_ind = find(obj.n_superpixels == n_supers);
    else
        obj.n_superpixels = [obj.n_superpixels, n_supers];
        n_superpix_ind = find(obj.n_superpixels == n_supers);
    end

    if numel(obj.class_labels) <= n_superpix_ind || isempty(obj.class_labels{n_superpix_ind})
        obj.class_labels{n_superpix_ind} = cell(max(label_mat(:)),1);
    end

    % now we should ginput the superpixels that represent our classes
    % that we want
    % gotta loop through all classes
    for i = 1:n_classes
        
        % make an empty array to hold the pixel numbers for this class
        this_class = [];

        % show the image 
        fig = figure();
        % have the figure fill the screen outright
        fig.WindowState = 'maximized';
        imshow(imoverlay(rgb_im,bound_mask,'cyan'));
        title("Click superpixels for " + obj.classes_present{i} + ". Press enter to end.");
        hold on
        
        % to indicate clicked superpixels, we'll need an alpha mask of same
        % size as the image
        alpha_mask = ones(size(rgb_im,1), size(rgb_im,2));
        mask_influencer = zeros(size(rgb_im,1), size(rgb_im,2));
        mask_show = imshow(alpha_mask);

        % start with the mask alphadata as the zeros influence array
        mask_show.AlphaData = mask_influencer;

        % we might have some existing ids to block off already
        if sum(strcmp(obj.class_labels{n_superpix_ind},obj.classes_present{i})) > 0 
            % get the list of indices that have this label already
            to_fill = find(strcmp(obj.class_labels{n_superpix_ind},obj.classes_present{i}));
            % gotta look at the superpixels and update the mask influencer
            mask_influencer(find(ismember(imresize(label_mat,scaling,'nearest'), to_fill))) = 1;
            mask_show.AlphaData = mask_influencer;
        end
        
        % make a breakable loop so we can stop clicking
        while true

            % ginput a super pixel
            [col,row,button] = ginput(1);

            % if we didn't click, break
            if isempty(button)
                break; 
            % if the user presses i, they want to zoom in     
            elseif button == 105
                % get the axes and dimensions
                ax = axis; 
                width=ax(2)-ax(1);
                height=ax(4)-ax(3);
                % reset the dimensions
                axis([col-width/2 col+width/2 row-height/2 row+height/2]);
                % and zoom in    
                zoom(2); 
            
            % if the user presses o, they want to zoom out
            elseif button == 111
                % get the axes and dimensions
                ax = axis; 
                width=ax(2)-ax(1);
                height=ax(4)-ax(3);
                % reset the dimensions
                axis([col-width/2 col+width/2 row-height/2 row+height/2]);
                % and zoom out    
                zoom(1/2); 
            
            % otherwise the user wants to click a superpixel, so handle it    
            else
                % have to add an extra if statement here in case the user
                % clicked outside the image 
                if round(row)/scaling < size(label_mat,1) && round(col)/scaling < size(label_mat,2) && round(row)/scaling > 0 && round(col)/scaling > 0
                    % figure out which number superpixel was clicked
                    clicked_pixel = label_mat(round(row)/scaling,round(col)/scaling);
        
                    % update pixel list and mask only if this is a new pixel
                    if ~ismember(clicked_pixel,this_class)
                        this_class = [this_class;clicked_pixel];
        
                        % gotta look at the superpixels and update the mask influencer
                        mask_influencer(find(ismember(imresize(label_mat,scaling,'nearest'), clicked_pixel))) = 1;
                        mask_show.AlphaData = mask_influencer;
                    else
                        % if it's already in the list, the user wants to unclick,
                        % so remove it
                        this_class(this_class == clicked_pixel) = [];
        
                        % and update the alpha mask stuff
                        mask_influencer(find(ismember(imresize(label_mat,scaling,'nearest'), clicked_pixel))) = 0;
                        mask_show.AlphaData = mask_influencer;
                    end
                end
            end
        end

        % let the user know we're exporting everything 
        fprintf('Handling the data.\n')

        %so we can clear the figure for the next class
        hold off
        close(fig)
        % put the numbers for this class into the string array for output
        obj.class_labels{n_superpix_ind}(this_class) = obj.classes_present(i);

        %% last step is to take superpixel groups and make them labeled regions stored in json format
    
        % make the mask from all currently labeled superpixels
        region_mask = ismember(label_mat, this_class);

        % get the boundaries
        region_bounds = bwboundaries(region_mask);
        
        % gotta worry about holes, but only if there is more than one
        % boundary
        if numel(region_bounds) > 1
            % need to decide if any of these are holes by looking for overlap
            % first construct a polyshape vector
            polyvec = repmat(polyshape,1,numel(region_bounds));
            warning('off','all')
            for j = 1:numel(region_bounds)
                polyvec(j) = polyshape([region_bounds{j}(:,1), region_bounds{j}(:,2)]);
            end
            warning('on','all');

            % and test for overlap
            is_overlapping = overlaps(polyvec);
            
            % the diagonal of the overlap matrix is automatically true, but
            % we want it to be false. so get the linear indices of the
            % diagonal and set them to false
            diag_inds = linspace(1,numel(is_overlapping),size(is_overlapping,1));
            is_overlapping(diag_inds) = false;
    
            % if we do have overlap, need to do something about it
            if sum(is_overlapping,'all') > 0
                holes = false(1,numel(region_bounds));
                % because the overlapping matrix is symmetric, we only need to
                % iterate over one dimension
                hole_pairs = nchoosek(find(any(is_overlapping)),2);
                for j = 1:size(hole_pairs,1)
                    % we'll use areas to decide which one is the hole
                    area1 = area(polyvec(hole_pairs(j,1)));
                    area2 = area(polyvec(hole_pairs(j,2)));
                    [~,hole_ind] =  min([area1, area2]);
                    % and set the proper one as the hole in the logical
                    holes(hole_pairs(j,hole_ind)) = true;
                end
            else
                % in this case no holes
                holes = false(1,numel(region_bounds));
            end
        else
            holes = false(1,1);
        end

        % we have boundaries and holes delineated...ready to put the
        % attributes in the class. Do so as a struct for each region
        for j = 1:numel(region_bounds)
            region_prop = struct;
            region_prop.class = obj.classes_present{i};
            region_prop.col_coords = region_bounds{j}(:,2);
            region_prop.row_coords = region_bounds{j}(:,1);
            region_prop.is_hole = holes(j);
    
            % and put the struct into the overall struct array for the instance
            obj.class_regions = [obj.class_regions, region_prop];
        end
    end
end