function explore_plot(obj, stat_num, color_scheme)

    % the big plot on the left will be a basemap with all of the
    % samples scattered, color coded by the value of the stat
    % assigned

    % and we don't want to include background superpixels in this, so id
    % the background superpixels for removal from mean
    is_background = cellfun(@(x)strcmpi(x,'background'), obj.class_labels, 'UniformOutput', false);

    % need the mean of all stats to color code the 
    not_background = cellfun(@(a,b) a(~b,:), obj.superpix_stats, is_background, 'UniformOutput', false);
    mean_stat = cell2mat(cellfun(@(x)mean(x(:,stat_num)), not_background, 'UniformOutput', false));
    
    % and call in a figure of the right size 
    figure('units','normalized','outerposition',[0 0 0.9 1])
    % we want the colormap information for later, depending on if the user
    % wants a diverging or discrete colormap
    if strcmpi(color_scheme,'diverging')
        burd_raw = [33, 102, 172; 67, 147, 195; 146, 197, 222; 209, 229, 240; ...
            247, 247, 247; 253, 219, 199; 244, 165, 130; 214, 96, 77; 178, 24, 43]./255;
        cmap = imresize(burd_raw,[256, 3], 'bilinear');
    elseif strcmpi(color_scheme,'discrete')
        disc_raw = [209, 187, 215; 174, 118, 163; 136, 46, 114; 25, 101, 176; ...
            82, 137, 199; 123, 175, 222; 78, 178, 101; 144, 201, 135; ...
            202, 224, 171; 247, 240, 86; 246, 193, 65; 241, 147, 45; 232, 96, 28;...
            220, 5, 12]./255;
        cmap = imresize(disc_raw,[256, 3], 'bilinear');
    else
        error('not a valid color scheme request')
    end

    % we'll want the field facies loaded in for later
    load('field_facies.mat');

    % before we start plotting, need to get the dem loaded and adjusted
    % a bit.
    [dem_im, dem_ref] = readgeoraster(obj.dem_path);

    % everything likely is in lat long, so we need to transform some parts
    % of our reference frame to utm for easy plotting
    % get the geoid of the zon and construct its projection structure
    utm_zone = utmzone(dem_ref.LatitudeLimits(1), dem_ref.LongitudeLimits(1));
    [ellipsoid,estr] = utmgeoid(utm_zone);
    utmstruct = defaultm('utm');
    utmstruct.zone = utm_zone;
    utmstruct.geoid = ellipsoid;
    utmstruct = defaultm(utmstruct);
    % get converted image limits in utm
    [utm_xlims, utm_ylims] = mfwdtran(utmstruct,dem_ref.LatitudeLimits,dem_ref.LongitudeLimits);

    % can also use those limits and the size of the dem to get scale of the
    % dem
    dem_xscale = size(dem_im,2)/range(utm_xlims);
    dem_yscale = size(dem_im,1)/range(utm_ylims);

    % we don't want to show all of the dem likely, so we need the limits of
    % the data to be plotted, with a 10m buffer on all sides
    buffer_size = 10;
    data_xmin = min(obj.utm_xyz(:,1)) - buffer_size;
    data_xmax = max(obj.utm_xyz(:,1)) + buffer_size;
    data_ymin = min(obj.utm_xyz(:,2)) - buffer_size;
    data_ymax = max(obj.utm_xyz(:,2)) + buffer_size;

    % now just need to get where those numbers exist in pixel space on the
    % dem
    imshow_start_row = round((max(utm_ylims) - data_ymax)*dem_yscale);
    imshow_end_row = round((max(utm_ylims) - data_ymin)*dem_yscale);
    imshow_start_col = round((data_xmin- min(utm_xlims))*dem_xscale);
    imshow_end_col = round((data_xmax - min(utm_xlims))*dem_xscale);

    % in the dem, background pixels are set to -32767, so reset those to
    % nan
    dem_im(dem_im == -32767) = NaN;
    
    % and grab the part to use as the basemap and rescale it
    to_show = imresize(dem_im(imshow_start_row:imshow_end_row, imshow_start_col:imshow_end_col), 1/mean([dem_xscale,dem_yscale]),'bicubic');
    % last thing is to set up an alpha channel to block NaNs
    dem_alpha = ones(size(to_show));
    dem_alpha(isnan(to_show)) = 0;


    % start by showing the image on custom axes to mirror aspect ratio of the image
    max_height = 0.6;
    max_width = max_height * (size(to_show,2)/size(to_show,1));
    axes_positions = [0.1 0.15 max_width max_height];
    ax1 = axes('InnerPosition', axes_positions);
    imagesc(ax1, to_show, 'AlphaData',dem_alpha);
    set(ax1,'color',0*[1 1 1])
    axis image
    hold on
    % make the scatter 
    
    ax2 = axes('InnerPosition',axes_positions);
    hold(ax2,"on");
    % and before we put the markers for samples down, we can place rough
    % polygons for geologic boundaries that exist
    geo_vert_dir = dir(fullfile(obj.geo_vertex_path, '*.csv'));
    % and we'll need the color order to keep polygons of the same type of
    % the same color
    c_ord = get(gca,'colororder');
    % only try if there are files there
    if ~isempty(geo_vert_dir)
        for i = 1:numel(geo_vert_dir)
            these_verts = readmatrix(fullfile(obj.geo_vertex_path, geo_vert_dir(i).name));
            % might have multiple polygons identified for these boundaries,
            % so loop through the ids stored in the 3rd column
            for j = 1:max(these_verts(:,3))
                % and we should already be in utm, so just plot, adding in
                % the first point at the end to close the loop
                this_poly = find(these_verts(:,3) == j);
                geo_poly = plot([these_verts(this_poly, 1); these_verts(this_poly(1), 1)], [these_verts(this_poly, 2); these_verts(this_poly(1), 2)]-203,'Color',c_ord(i,:));
                set(get(get(geo_poly,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end
    
    % we'll loop through the field facies and scatter each one with a
    % different marker
    filled_markers = {"o", "square", "diamond", "^", "v", ">", "<"};
    for i = 1:numel(field_facies)
        stat_frac = (mean_stat(obj.field_lithos == i)-min(mean_stat))/range(mean_stat);
        cmap_ind = round(size(cmap,1)*stat_frac);
        cmap_ind(cmap_ind == 0) = 1;
        scatter(ax2,obj.utm_xyz(obj.field_lithos == i,1), obj.utm_xyz(obj.field_lithos == i,2),40,cmap(cmap_ind, :),'filled', ...
            filled_markers{i},'MarkerEdgeColor',[.1,.1,.1], 'DisplayName',field_facies{i});
    end
    hold(ax2,"off");

    axis image
    % refine some stuff on the second axes for visibility and geography
    ax2.Visible = 'off';
    ylim(ax2,[data_ymin data_ymax])
    xlim(ax2, [data_xmin data_xmax])
    % refine a bunch of things about the first axes to not show up
    set(ax1,'YTick',[])
    set(ax1,'XTick',[])
    set(ax1,'YTicklabel',[])
    set(ax1,'XTicklabel',[])
    grid(ax1,'off')

    % and some appearance things and individual colormaps
    xlabel(ax2,'easting [m]')
    ylabel(ax2,'northing [m]')
    % don't want scientific notation
    ax2.YAxis.Exponent = 0;
    ax2.XAxis.Exponent = 0;
    xtickformat(ax2,'%.0f')
    ytickformat(ax2,'%.0f')
    % customize the gray colormap to lighten the dem a bit
    gray_cmap = colormap(gray);
    colormap(ax1, gray_cmap(75:end,:));
    colormap(ax2, cmap);
    
    % make sure we have a color bar and legend
    c = colorbar(ax2, 'Position', [(axes_positions(2) + axes_positions(3) -0.05), ax1.InnerPosition(1)+0.05, 0.02, ax1.InnerPosition(4)]);
    c.Label.String = obj.stat_titles{stat_num};
    legend('Location','northwest')


    % Set up crosshairs on each axis at the edges
    gobj(1,1) = xline(ax2,min(xlim(ax2)), 'w-');
    gobj(1,2) = yline(ax2,min(ylim(ax2)), 'w-');

    % but don't have them show up in the legend
    set(get(get(gobj(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(gobj(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    % Assign windowbuttonmotion fcn on axis #1
    set(ax2.Parent,'windowbuttonmotionfcn', {@mouseMove, ax2, gobj, obj, cmap, mean_stat, stat_num, field_facies});

    % Assign mouse button functions to start/stop tracking
    WindowButtonMotionFcnInput = {@mouseMove, ax2, gobj, obj, cmap, mean_stat, stat_num, field_facies};
    set(ax2.Parent,'windowbuttondownfcn', {@startStopMouseMove, WindowButtonMotionFcnInput})

end


function mouseMove(~, ~, ax, gobj, ~, ~, ~, ~, ~)
    % Responds to mouse movement in axis #1
    % ax is a vector of subplot handles; ax1 tracks mouse movement, ax(2) follows.
    % gobj(1,1) is xline in ax 1
    % gobj(1,2) is yline in ax 1
    % Get mouse coordinate
    C = ax.CurrentPoint;
    x = C(1,1);
    y = C(1,2);
    % If mouse isn't on axis #1, do nothing.
    if x < ax.XLim(1) || x > ax.XLim(2) || y < ax.YLim(1) || y > ax.YLim(2)
        return
    end
    % Update crosshairs (cross hairs in ax 2 are yoked).
    gobj(1,1).Value = x;
    gobj(1,1).Label = x;
    gobj(1,2).Value = y;
    gobj(1,2).Label = y;
end

function startStopMouseMove(hobj,~,WindowButtonMotionFcnInput)
    % Turns mouse tracking off (right mouse button) and on (left mouse button)
    buttonID = hobj.SelectionType;
    switch buttonID
        case 'alt'  %right mouse button
            % Start interactivity
            set(hobj,'windowbuttonmotionfcn', WindowButtonMotionFcnInput);
        case 'normal'     % left mouse button
            % Stop interactivity
            set(hobj,'windowbuttonmotionfcn', []);

            % this also is the case where have 'selected' a point. So we
            % should do the other subplotting

            % start by grabbing the closest point. 
            % the final crosshair values are stored in the constant line
            % objects in the 3rd cell of the Window..Input array
            x_clicked = get(WindowButtonMotionFcnInput{3}(1), 'Value');
            y_clicked = get(WindowButtonMotionFcnInput{3}(2), 'Value');

            % and we've passed the petro_superpix object in as the 4th cell
            % of the input arry, so we can look at the x-y locations to see
            % which one we're closest to
            obj = WindowButtonMotionFcnInput{4};
            nearest_point = knnsearch(obj.utm_xyz(:,[1,2]),[x_clicked,y_clicked],'K',1);

            % also need to take the color map out of the input
            % array
            cmap = WindowButtonMotionFcnInput{5};
            % to plot just this sample as a different symbol but within the
            % same color scheme, where it falls in the range of this
            % dataset
            mean_stat = WindowButtonMotionFcnInput{6};
            stat_frac = (mean_stat(nearest_point)-min(mean_stat))/range(mean_stat);
            cmap_ind = round(size(cmap,1)*stat_frac);
            % sometimes the cmap_ind is zero for the min, so correct it to 1
            if cmap_ind <= 0
                cmap_ind = 1;
            end

            % we may be reentering this loop having already selected one
            % point. So if there is a pentagram on the current axes, we
            % have to delete it

            pentagram_obj = findobj('Marker','pentagram');
            if ~isempty(pentagram_obj)
                delete(pentagram_obj)
            end

            hold(WindowButtonMotionFcnInput{2},'on')
            % put a star on the initial plot given the input clicked
            % location
            title_string = char(join(split(obj.sample_names{nearest_point}, '_')));
            single_point = scatter(WindowButtonMotionFcnInput{2},obj.utm_xyz(nearest_point,1), obj.utm_xyz(nearest_point,2),500,cmap(cmap_ind,:),...
                'filled','pentagram','MarkerEdgeColor',[.1,.1,.1], 'DisplayName', title_string);
            legend('Location','northwest')
            
            % now we're also clear to make the two subplots on the right
            % upper right will be the rgb image with some high and low
            % superpixels for the stat at hand marked
            rgb_im = imread(fullfile(obj.low_res_path,[obj.sample_names{nearest_point} '.jpg']));
            superpixel_inds = imresize(imread(fullfile(obj.superpixel_ind_path, [obj.sample_names{nearest_point}, '.tif'])),0.05,'nearest');
            ind_col = superpixel_inds(:);
            stat_num = WindowButtonMotionFcnInput{7};
            stat_col = obj.superpix_stats{nearest_point}(ind_col,stat_num);

            % also need to account for background
            background_supers = find(strcmpi(obj.class_labels{nearest_point},'background'));
            stat_col(ismember(ind_col,background_supers)) = NaN;
            stat_im = reshape(stat_col,size(superpixel_inds,1),size(superpixel_inds,2));
            im_alpha = ones(size(stat_im));
            im_alpha(isnan(stat_im)) = 0;
    
   

            % also calculate the superpixels that have the highest value
            % for this stat to display on the rgb image
            % start by extracting the stat and setting background to nan
            this_stat = obj.superpix_stats{nearest_point}(:,stat_num);
            this_stat(background_supers) = NaN;
            
            % and get top and bottom 5
            [~, top_5] = maxk(this_stat, 5);
            [~, bot_5] = mink(this_stat, 5);

            % need the scale ratio between the superpixel indices and the
            % rgb image for display of the highest and lowest value
            % superpix
            scale_ratio = size(superpixel_inds,1)/size(rgb_im,1);

            % show the thing
            subplot(2,2,2)
            imshow(rgb_im)
            title(title_string)
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
            % calculated earlier. So show that as well
            super_sc = subplot(2,2,4);
            imagesc(stat_im, 'AlphaData',im_alpha)
            axis image
            grid off
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            set(gca,'color',0*[1 1 1])
            colormap(super_sc, WindowButtonMotionFcnInput{5})
            colorbar('southoutside')
            % make the field facies the title
            field_facies = WindowButtonMotionFcnInput{8};
            title(['field facies: ' field_facies{obj.field_lithos(nearest_point)}])

            % and add text to the figure for categorical angle of imgage
            % wrt bedding
            axes('Position',[.1 .8 .4 .1])
            grid off
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            set(gca,'YTick',[])
            set(gca,'XTick',[])
            xlim([0.5 10])
            ylim([0.5 1.5])
            if isnan(obj.im_bedding_angles(nearest_point))
                text(1,1,'Image orientation not available.', 'FontSize', 18)
            elseif obj.im_bedding_angles(nearest_point) < 15 || obj.im_bedding_angles(nearest_point) > 165
                text(1,1,'Image is approximately along bedding plane.', 'FontSize', 18)
            elseif obj.im_bedding_angles(nearest_point) > 75 && obj.im_bedding_angles(nearest_point) < 105
                text(1,1,'Image is a cross-section wrt bedding.', 'FontSize', 18)
            else
                text(1,1,'Image is oblique wrt bedding.', 'FontSize', 18)
            end
        
            % if this stat contains the word filter, show that filter
            if contains(obj.stat_titles{stat_num}, 'filter')
                axes('Position',[.5 .45 .1 .1])
                filter_split = split(obj.stat_titles{stat_num}, 'filter ');
                filter_num = str2num(filter_split{2}(1));
                imagesc(obj.filter_bank(:,:,filter_num))
                axis image
                grid off
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);
            end
    end
end