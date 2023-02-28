classdef petro_superpix < handle
    % straight up properties the user can declare
    properties
        % cell array of sample names in order
        sample_names {mustBeText} = '';
        % wavelenghts for the images, in order that the paths occur
        wavelengths(1, :) {mustBeNumeric} = [365, 470, 505, 530, 590, 625, 760, 940];
        % cell array with the sub paths to the image folders corresponding
        % the strike of the image planes in the real world
        im_strikes {mustBeNumeric}
        % dip of the image planes in the real world
        im_dips {mustBeNumeric}
        % the field lithology ids of the sample
        field_lithos
        % pointing of the arrow of strike in the image (NSEW style)
        strike_im_headings {mustBeNumeric}
        % the strike of local bedding planes in the field
        local_strikes {mustBeNumeric}
        % dip of the local bedding planes in the field
        local_dips {mustBeNumeric}
        % latitude and longitude location of the image
        utm_xyz (:,3)
        % subpath to superpixel label matrices
        superpixel_ind_path {mustBeText} = '';
        % path to low resolution images for display
        low_res_path {mustBeText} = '/Users/ryan/Dropbox (Princeton)/reef_survey_project/nevada/nevada_qgis/rgb_jpgs';
        % how many superpixels we use
        n_superpixels {mustBeInteger}
        % stats for all of the superpixels 
        superpix_stats cell
        % the titles for all of the stats
        stat_titles {mustBeText} = ''
        % the filter bank used for responses when getting stats
        filter_bank = makeLMfilters(6);
    end

    % all of the functions we can perform with this class
    methods
        % a function that makes a really nice plot for exploring statistics
        function explore_plot(obj, im_num, stat_num)
            % the big plot on the left will be a basemap with all of the
            % samples scattered, color coded by the value of the stat
            % assigned
            % need the mean of all stats to color code the dots 
            mean_stat = cell2mat(cellfun(@(x)mean(x(:,stat_num)), obj.superpix_stats, 'UniformOutput', false));
            % to plot just this sample as a different symbol but within the
            % same color scheme, where it falls in the range of this
            % dataset
            stat_frac = (mean_stat(im_num)-min(mean_stat))/range(mean_stat);
            figure()
            subplot(1,2,1)
          
            % make the scatter
            hold on
            h = scatter(obj.utm_xyz(:,1), obj.utm_xyz(:,2),40,mean_stat,'filled', ...
                'MarkerEdgeColor',[.1,.1,.1]);
            xlabel('easting [m]')
            ylabel('northing [m]')
            % don't want scientific notation
            ax = ancestor(h, 'axes');
            ax.YAxis.Exponent = 0;
            ax.XAxis.Exponent = 0;
            xtickformat('%.0f')
            ytickformat('%.0f')
            axis image
            colorbar

            % and add in the single point for the current image
            cmap = colormap;
            cmap_ind = round(size(cmap,1)*stat_frac);
            scatter(obj.utm_xyz(im_num,1), obj.utm_xyz(im_num,2),500,cmap(cmap_ind,:),...
                'filled','pentagram','MarkerEdgeColor',[.1,.1,.1])
            hold off
            
            % upper right will be the rgb image with some high and low
            % superpixels for the stat at hand marked
            rgb_im = imread(fullfile(obj.low_res_path,[obj.sample_names{im_num} '.jpg']));
            superpixel_inds = imresize(imread(fullfile(obj.superpixel_ind_path, [obj.sample_names{im_num}, '.tif'])),0.05,'nearest');
            ind_col = superpixel_inds(:);
            stat_col = obj.superpix_stats{im_num}(ind_col,stat_num);
            stat_im = reshape(stat_col,size(superpixel_inds,1),size(superpixel_inds,2));
            
            [~, top_5] = maxk(obj.superpix_stats{im_num}(:,stat_num), 5);
            [~, bot_5] = mink(obj.superpix_stats{im_num}(:,stat_num), 5);
            scale_ratio = size(superpixel_inds,1)/size(rgb_im,1);
            subplot(2,2,2)
            imshow(rgb_im)
            hold on
            top_mask = ismember(superpixel_inds,top_5);
            top_bounds = bwboundaries(top_mask);
            for i = 1:numel(top_bounds)
                patch(top_bounds{i}(:,2)./scale_ratio,top_bounds{i}(:,1)./scale_ratio,'white','FaceColor', 'none', 'EdgeColor', cmap(end-10,:),'LineWidth', 2);
                drawnow
            end
            bot_mask = ismember(superpixel_inds,bot_5);
            bot_bounds = bwboundaries(bot_mask);
            for i = 1:numel(bot_bounds)
                patch(bot_bounds{i}(:,2)./scale_ratio,bot_bounds{i}(:,1)./scale_ratio,'white','FaceColor', 'none', 'EdgeColor', cmap(10,:),'LineWidth', 2);
                drawnow
            end
            hold off

            % lower right is the superpixel values for the stats as
            % calculated earlier
            subplot(2,2,4)
            imagesc(stat_im)
            axis image
            grid off
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            colorbar
        end
       
        % a function to make all of the stat titles. A bit fixed based upon
        % the superpixel stat function. Didn't know where to put it, but
        % it's here for now
        function obj = make_stat_titles(obj)
            % start with centroids
            stat_titles = {'centroid row'; 'centroid column'};

            % add in the wavelenghth means
            for i = 1:numel(obj.wavelengths)
                stat_titles = [stat_titles; [num2str(obj.wavelengths(i)) 'nm mean']];
            end

            % and the GLCM stats
            glcm_names = {'GLCM contrast', 'GLCM correlation', 'GLCM energy', 'GLCM homogeneity'};
            for i = 1:numel(obj.wavelengths)
                for j = 1:numel(glcm_names)
                    stat_titles = [stat_titles; [num2str(obj.wavelengths(i)) 'nm ' glcm_names{j}]];
                end
            end

            % and the filter responses
            n_filters = size(obj.filter_bank,3);
            for i = 1:n_filters
                for j = 1:numel(obj.wavelengths)
                    stat_titles = [stat_titles; [num2str(obj.wavelengths(j)) 'nm - filter ' num2str(i) ' response']];
                end
            end

            obj.stat_titles = stat_titles;
        end
    end
end