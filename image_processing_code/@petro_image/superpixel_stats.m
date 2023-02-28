function superpixel_stats(obj, n_supers)
% This function takes in a petro_image object and a specified number(s) of
% superipixels, and fills in many statistics for in the superpix_stats
% property of the object. The stats are a matrix in a cell array, which has
% a row for each superpixel in the image, and a column for each statistic.
% Each entry in the cell array is for a different defined number of
% superpixels on the same image, matching the indices stored in the
% n_superpixels property. If the indices for the superpixel haven't been
% made and stored yet, the function makes based upon an RGB image and saves
% them in the specified path for the object.
%
% IN:
% obj: instance of the petro_image object class
% n_supers: vector containing the desired number of superpixels in the
% oversegmentation to get statistics for
%
% OUT:
%
% none, but the object from the input is updated
%
% R. A. Manzuk 
% written: Friday, February 17, 2023 at 3:08:23 PM
    %% begin the function
    % process is a bit lengthy, so don't want to recreate stats that
    % already exist. So check for that and return
    [~, super_inds] = intersect(obj.n_superpixels, n_supers);
    if numel(obj.superpix_stats) >= max(super_inds)
        has_stats = ~cellfun(@isempty, obj.superpix_stats(super_inds));
        if sum(has_stats) == numel(n_supers)
            disp('All of those stats already exist for that object.')
            return
        end
    end

    % start by reading in all available wavelengths. 
    % read in the first for dimensions. then we can concatenate all others
    disp('reading in all available image channels');
    multispec_im = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{1}, [obj.sample_name, obj.default_ext])));
    for i = 2:numel(obj.im_subpaths)
        this_channel = im2double(imread(fullfile(obj.main_path, obj.im_subpaths{i}, [obj.sample_name, obj.default_ext])));
        multispec_im = cat(3,multispec_im,this_channel);
    end

    % the multispec_im by itself will be good enough for color stuff. but
    % for GLCM, we need a binned, gray level image
    disp('Calculating gray level co-occurrences');
    n_gray_levels = 8; % can make this more flexible in the future
    bin_edges = linspace(0,1,n_gray_levels+1);

    binned_image = discretize(multispec_im,bin_edges);

    % and we'll be comparing each pixel to the pixel to the right, so we
    % can make a comparison matrix that's the same thing but shifted over,
    % padded with ones
    comp_image = [binned_image(:,2:end,:), ones(size(binned_image,1),1,obj.num_channels)];

    % to easily get the linear indices of each pixel's GLCM assignment, we
    % can just make a 2 column matrix of the binned and comp images side by
    % side, and get its subindices
    glcm_inds = sub2ind([n_gray_levels n_gray_levels],binned_image(:),comp_image(:));
 
    % and just reshape that to be the original size of the image
    for_glcm = reshape(glcm_inds,size(multispec_im,1),size(multispec_im,2),[]);

    % last thing to get stats from is activations of a filter bank on a
    % normalized, downsampled version of the image
    disp('Getting filter bank responses');
    for_filt = imresize(multispec_im,0.1);
    for i = 1:size(for_filt,3)
        for_filt(:,:,i) = normalize_image(for_filt(:,:,i), 'windows', 6);
    end

    
    % and set up an empty response matrix for all filters and all 8
    % channels, then iterate through and get the filter responses
    f_responses = zeros(size(for_filt,1),size(for_filt,2),size(obj.filter_bank,3)*obj.num_channels);
    for i = 1:size(f_bank,3)
        f_responses(:,:,((i-1)*size(for_filt,3)+1):i*size(for_filt,3)) = imfilter(for_filt,obj.filter_bank(:,:,i));
    end
    
    % time to iterate through all the superpixel numbers the user asked for
    for i = 1:numel(n_supers)
        % takes some time to assemble superpix stats, so check they don't
        % exist for this n_superpix
        if ismember(n_supers(i), obj.n_superpixels)
            this_ind = find(obj.n_superpixels == n_supers(i));
            if this_ind <= numel(obj.superpix_stats) && ~isempty(obj.superpix_stats(this_ind))
                disp(['Stats already exist for ' num2str(n_supers(i)) '. Moving on.'])
                continue
            end
        end

        % now we have to check if there are superpixels already made
        superpix_fname = fullfile(obj.main_path, obj.superpixel_subpath, num2str(n_supers(i)), [obj.sample_name, obj.default_ext]);
        if isfile(superpix_fname)
            disp('loading in superpixel indices');
            % if they are there we can just load them
            label_mat = imread(superpix_fname);
        else
            % otherwise gotta make them
            % first get the red green and blue inds for the multispec im to
            % input to superpixel making
            red_ind = find(obj.wavelengths == 625);
            green_ind = find(obj.wavelengths == 530);
            blue_ind = find(obj.wavelengths == 470);
            label_mat = obj.get_superpix(n_supers(i), multispec_im(:,:,[red_ind,green_ind,blue_ind]));
        end
        
        % now we just loop through each superpixel and gather its stats
        disp(['Assembling all of the stats for ' num2str(n_supers(i)) ' supers.']);
        % preallocate the size of the stat mat for speed. 2 elements for
        % centroid, obj.num_channels elements for color, nchan*4 elements for glcm, and
        % 3rd dimensions size of f_responses for mean color response.
        n_glc = 4;
    
        % to remove as much from the main loop as possible, we'll work with the
        % color and filter response means in columns   
        label_vec = label_mat(:);
        col_image = reshape(multispec_im,[],obj.num_channels);
        color_means = splitapply(@(x)mean(x, 1), col_image,label_vec);
        
        % work with the label vector to get centroids
        [label_rows,label_cols] = ind2sub(size(label_mat), 1:numel(label_vec));
        centroids = [splitapply(@(x)mean(x, 1),label_rows', label_vec), splitapply(@mean,label_cols', label_vec)];
        
        % downsample the label vector to get the mean filter responses on a
        % column version of the filter responses
        response_col = reshape(f_responses,[],size(f_responses,3));
        label_vec_small = imresize(label_vec,[size(response_col,1),1],'nearest');
        mean_responses = splitapply(@(x)mean(x, 1),response_col,label_vec_small);
    
        % all we need is glcm stuff
        glc_col = reshape(for_glcm,[],obj.num_channels);
        glc_stats = splitapply(@get_grayprops, glc_col, label_vec);

        % to output things, gotta make sure the number of superpixels used
        % is in the object's list
        if ismember(n_supers(i), obj.n_superpixels)
            n_superpix_ind = find(obj.n_superpixels == n_supers(i));
        else
            obj.n_superpixels = [obj.n_superpixels, n_supers(i)];
            n_superpix_ind = find(obj.n_superpixels == n_supers(i));
        end
    
        % and wrap it all up in one matrix output
        obj.superpix_stats{n_superpix_ind} = [centroids, color_means, glc_stats, mean_responses];

    end

    function glc_stats = get_grayprops(gray_counts)
        glc_stats = zeros(1,size(gray_counts,2)*4);
        for k = 1:size(gray_counts,2)
            glc_counts = histcounts(gray_counts(:,k),0.5:1:(8^2)+1);
            this_glcm = reshape(glc_counts,8,8);
            glc_stats_stuct = graycoprops(this_glcm);
            glc_stats(((k-1)*4)+1:k*4) = [glc_stats_stuct.Contrast, glc_stats_stuct.Correlation, glc_stats_stuct.Energy, glc_stats_stuct.Homogeneity];
        end
    end
    

end


