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
        % dip of bedding with respect to the image plane
        im_bedding_angles {mustBeNumeric} 
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
        % path to the petro im instances
        petro_ims_path {mustBeText} = '/Volumes/ryan_ims/sms_images/petro_im_instances';
        % how many superpixels we use
        n_superpixels {mustBeInteger}
        % stats for all of the superpixels 
        superpix_stats cell
        % labels for all of the superpixels that we have so far
        class_labels cell
        % the titles for all of the stats
        stat_titles {mustBeText} = ''
        % the filter bank used for responses when getting stats
        filter_bank = makeLMfilters(6);
        % path to the dem to be used as a base map
        dem_path {mustBeText} = '/Users/ryan/Dropbox (Princeton)/reef_survey_project/nevada/nevada_qgis/drone_derived_dem.tif';
        % path to csv files with vertices of geo polygons  
        geo_vertex_path {mustBeText} = '/Users/ryan/Dropbox (Princeton)/reef_survey_project/nevada/nevada_qgis/basemap_vertices';
        % file name of sheet with carbon and oxygen data
        carb_ox_file {mustBeText} = '/Volumes/ryan_ims/sms_images/geochem/carbon_oxygen_master.csv';
        % file name of sheet with element concentration data
        element_file {mustBeText} = '/Volumes/ryan_ims/sms_images/geochem/element_concentrations_master.csv';
        % a cell array with all of the clicked geochem localities
        all_geochem_locs cell
        % the drilled phases have letters as codes, so just have the key as
        % part of the object
        phase_codes cell = {'A', 'ooid';...
            'B', 'recrystalized mud';...
            'C', 'dolomitized mud';...
            'D', 'calcite vein';...
            'E', 'mudstone';...
            'F', 'regular archaeo';...
            'G', 'irregular archaeo';...
            'H', 'calcimicrobe';...
            'I', 'shell';...
            'J', 'microbial';...
            'K', 'coralomorph'};
        % the superpixel indices of each geochem loc are stored in the
        % petro_image instances, have a property to store all of them here
        geochem_superpix_inds
        % the results of a PCA on geochemical data
        geochem_pca
        % the results of a PCA on superpixel stats
        superpixel_pca
        % a master table that holds geochem data and matched superpixel
        % values
        master_table
        % a within-sample mean of each geochemical variable. rows are
        % oganized in same order as sample_names
        sample_geochem_means
        % a within-sample variance of each geochemical variable. rows are
        % oganized in same order as sample_names
        sample_geochem_vars
        % can just read in available spreadsheets of geochem and organize
        % them
        all_geochem_data
    end

    properties (Dependent) 
        % based upon the utm coordinates of each sample and the dip of
        % bedding at the locality, we can assign a relative strat height to
        % each sample
        strat_height
    end

    % all of the functions we can perform with this class
    methods
        %%
        % load geochemical data in from spreadsheeds and make it into one
        % organized array
        function set_all_geochem_data(obj)
            % start with the carbon data
            % Load CSV file
            raw_carb_data = readtable(obj.carb_ox_file);
            
            % Extract labels from CSV file and split at the underscore
            split_carb_names = cellfun(@(x) strsplit(x,'_'), raw_carb_data.export_id, 'UniformOutput',false);

            % and get the inds for non split names (likely samples not
            % standards)
            is_carb_sample = cellfun(@numel, split_carb_names) > 1;

            % and get started with the trace elements data before long.
            % This csv has some weird formatting, so gotta do some stuff
            % with options when reading
            opts = detectImportOptions(obj.element_file);
            new_formats = cell(1, numel(opts.VariableTypes)-1);
            new_formats(:) = {'double'};
            opts.VariableTypes(2:end) = new_formats;
            opts.VariableNamingRule = 'preserve';
            raw_element_data = readtable(obj.element_file, opts);
            split_element_names = cellfun(@(x) strsplit(x,'_'), raw_element_data.name, 'UniformOutput',false);
            % sometimes have a third element of the split names due to a
            % note in the spreadsheet, only want the first two
            for i = 1:numel(split_element_names)
                if numel(split_element_names{i}) ~= 2
                    split_element_names{i} = split_element_names{i}(1:2);
                end
            end


            % okay, so now we can take just the carbon isotope and trace
            % elements data and get ready to match them up. Part of that is
            % making everything lowercase
            just_carb_samples = lower(vertcat(split_carb_names{is_carb_sample}));
            just_element_samples = lower(vertcat(split_element_names{:}));
            
            % for matching, we'll treat the carbon data as the template and
            % try to give the trace elements an index that corresponts to
            % the right row in the carbon data
            [is_match,match_inds] = ismember(strcat(just_element_samples(:,1), ' ', just_element_samples(:,2)),strcat(just_carb_samples(:,1), ' ', just_carb_samples(:,2)));

            % now we're ready to assemble a final table of geochem data
            element_headers = fieldnames(raw_element_data);
            good_elements = element_headers(12:36);

            % make up a table of just the elements data we want, put in the
            % right order for the carb table
            element_table = array2table(zeros(sum(is_carb_sample),numel(good_elements)),'VariableNames',good_elements);
            element_table(match_inds,:) = raw_element_data(is_match,12:36);
            
            % and just one final thing, let's put the clicked localities of
            % the drill marks into here.

            % to start, to get rid of the smg_ from the sample numbers of
            % the clicked localities.
            split_clicked = cellfun(@(x) strsplit(x,'_'), obj.all_geochem_locs(:,1), 'UniformOutput',false);
            split_clicked = vertcat(split_clicked{:});

            % now we can iterate through and put the row index from the
            % csv file that matches in a new vector
            chem_locs = zeros(sum(is_carb_sample),2);
            for i = 1:size(chem_locs,1)
                this_sample_no = just_carb_samples{i,1};
                this_phase = just_carb_samples{i,2};

                % get vectors for where the sample number and phase match
                sample_match = strcmpi(this_sample_no,split_clicked(:,2));
                phase_match = strcmpi(this_phase,obj.all_geochem_locs(:,2));

                % and combine those to see where we got a full match
                full_match = and(sample_match,phase_match);
                
                % and maybe we didn't get a hit because the 1 was missing
                % from the clicked locality phase, so try that if needed
                if sum(full_match) == 0
                    add_one = repmat({'1'},numel(obj.all_geochem_locs(:,2)),1);
                    phase_match = strcmpi(this_phase,strcat(obj.all_geochem_locs(:,2), ' ', add_one));
                    full_match = and(sample_match,phase_match);
                end
                
                % and now fill in match inds if we have a hit
                if sum(full_match) ~= 0
                    chem_locs(i,:) = obj.all_geochem_locs{full_match,3};
                else
                    chem_locs(i,:) = NaN;
                end
            end
            
            % and just make into a table for export
            loc_table = array2table(chem_locs,'VariableNames',{'im_location_col','im_location_row'});
            
            % and we want to keep the split names in separate columns, so
            % make a table for that
            split_names_table = array2table(just_carb_samples,'VariableNames',{'sample_number', 'phase'});
            % the final table is just a joining of all the tables
            obj.all_geochem_data = [split_names_table, raw_carb_data(is_carb_sample,2:end), element_table, loc_table];
        end
        %%
        % use utm coordinates and the bedding dip to assign a strat height
        % to each sample
        function value = get.strat_height(obj)
            % we'll use all the local bedding measurements to get a mean plane
            mean_dip = mean(obj.local_dips,'omitnan');
            
            % and for the mean strike, we'll have to use the circular mean
            % https://en.wikipedia.org/wiki/Circular_mean . Just in case
            % we're near 0 degrees. Start by going to cartesian coordinates
            % around the unit circle:
            [cart_x, cart_y] = pol2cart(obj.local_strikes(~isnan(obj.local_strikes)), ones(1,sum(~isnan(obj.local_strikes))));
            
            % and just take the mean of those positions and convert back to
            % an angle
            mean_x = mean(cart_x);
            mean_y = mean(cart_y);
            mean_strike = cart2pol(mean_x, mean_y);

            % for the rotation, we actually want to convert strike to dip
            % direction 
            mean_ddir = mean_strike + 90;
            if mean_ddir > 360
                mean_ddir = mean_ddir - 360;
            end 

            % we need a lever point to rotate everything around, for now
            % we'll just use the southermost point
            [~,lever_ind] = min(obj.utm_xyz(:,2));
            
            % correct points such that lever is the origin for rotation
            corrected_locations = obj.utm_xyz - obj.utm_xyz(lever_ind,:);
            
            %rotate around the z axis such that dip points exactly north
            z_rot_mat = [cosd(mean_ddir), -sind(mean_ddir), 0; sind(mean_ddir), cosd(mean_ddir), 0; 0, 0, 1];
            points_prime = z_rot_mat * corrected_locations';
            
            % and now snap points up on a lever to correct for dip (rotate
            % around x axis)
            x_rot_mat = [1, 0, 0; 0, cosd(mean_dip), -sind(mean_dip); 0, sind(mean_dip), cosd(mean_dip)];
            points_final = x_rot_mat * points_prime;
            
            % subtract the minimum new elevation to be point 0 for the
            % strat height
            value = (points_final(3,:) - min(points_final(3,:)))';
        end

        %% 
        % go within samples to assess mean and variance of geochemical
        % variables
        function set_within_samp(obj)   
            % not all of the geochem data in the table are relevant to have
            % a sample mean or variance, so need to select a subset
            all_vars = obj.all_geochem_data.Properties.VariableNames;
            bad_vars = {'sample_number', 'phase', 'line', 'delta_13C_std', ...
                'delta_18O_std', 'date_run', 'im_location_col', 'im_location_row'};
            bad_var_inds = ismember(all_vars, bad_vars);

            % so then we can set up arrays for mean and variance within
            % each sample
            sample_geochem_means = zeros(numel(obj.sample_names), sum(~bad_var_inds));
            sample_geochem_vars = zeros(numel(obj.sample_names), sum(~bad_var_inds));

            % we'll just do this by looping through the samples
            for i = 1:numel(obj.sample_names)
                % to compare to the sample names in the geochem data table,
                % we don't want the smg_, so crop that off and figure out
                % where these data are
                these_data_inds = strcmpi(obj.sample_names{i}(5:end), obj.all_geochem_data.sample_number);

                % all good, ready to assess the mean and variance of the
                % sample and place it in the array.
                sample_geochem_means(i,:) = mean(table2array(obj.all_geochem_data(these_data_inds, ~bad_var_inds)), 'omitnan');
                sample_geochem_vars(i,:) = var(table2array(obj.all_geochem_data(these_data_inds, ~bad_var_inds)), 'omitnan');
            end
            
            % and finish up by making a table and puting it in the object
            % property
            obj.sample_geochem_means =  array2table(sample_geochem_means, "VariableNames", all_vars(~bad_var_inds));
            obj.sample_geochem_vars =  array2table(sample_geochem_vars, "VariableNames", all_vars(~bad_var_inds));
            

        end
        %%
        % make the master table
        function set_master_table(obj)
            % at the outset, the master table is mostly the geochem table,
            % then we'll paste the superpixel stuff to it
            mast_tab = obj.all_geochem_data;

            % some of the variables in there, we don't really want for this
            % stage of data analysis, so clear them.
            mast_tab = removevars(mast_tab,{'line', 'delta_13C_std','delta_18O_std',...
                'im_location_col', 'im_location_row', 'date_run'});   

            % get the list of sample names and phases from the geochem
            % table to match superpix inds to         
            samped_samps = [obj.all_geochem_data.sample_number, obj.all_geochem_data.phase];
            [~, match_ind] = ismember(strcat(samped_samps(:,1),' ',samped_samps(:,2)), strcat(obj.geochem_superpix_inds(:,1), ' ', obj.geochem_superpix_inds(:,2)));
            super_inds = cell2mat(obj.geochem_superpix_inds(match_ind,3));
            
            % gotta prep for some name comparison with the sample names of
            % the superpixel stat cells
            split_names = cellfun(@(x) strsplit(x,'_'), obj.sample_names, 'UniformOutput',false);
            split_names_cell = vertcat(split_names{:});

            % now we're ready to go through and sample the superpixels to
            % get our matrix of superpix that were sampled for geochem. And
            % while we're looping through, might as well get associated
            % strat heights
            super_stats = zeros(size(mast_tab,1), numel(obj.stat_titles));
            strat_heights = zeros(size(mast_tab,1), 1);
            for i = 1:numel(split_names)
                has_geochem = strcmpi(samped_samps(:,1), split_names_cell{i,2});
                if sum(has_geochem) > 0
                    super_stats(find(has_geochem), :) = obj.superpix_stats{i}(super_inds(has_geochem),:);
                    strat_heights(find(has_geochem), :) = repmat(obj.strat_height(i), sum(has_geochem), 1);
                end
            end   

            % add the first 4 PCs from image stat and geochem pca stuff
            % first make sure we have geochem and superpixel pcas to
            % display. if not, make them.
            if isempty(obj.geochem_pca)
                obj.set_geochem_pca;
            end
                
            if isempty(obj.superpixel_pca)
                obj.set_superpixel_pca(10000, [3, 4, 5, 6, 7, 8, 9, 10, 23, 26]);
            end
            
            % now we have to match the phases and sample numbers that were
            % inclueded in the pca to the right rows on the master table
            unique_samples = unique(obj.geochem_pca.samp_nums);
            pc_scores = zeros(size(mast_tab, 1), 8);
            for i = 1:numel(unique_samples)
                this_sample_ind_pca = strcmpi(obj.geochem_pca.samp_nums, unique_samples{i});
                this_sample_ind_table = strcmpi(mast_tab.sample_number, unique_samples{i});
                % then find where the phase labels from the pca match the
                % table
                [is_match, match_ind] = ismember(strcat(mast_tab.sample_number(this_sample_ind_table),' ',mast_tab.phase(this_sample_ind_table)),...
                    strcat(obj.geochem_pca.samp_nums(this_sample_ind_pca), ' ', obj.geochem_pca.phases(this_sample_ind_pca)));
                % then just fill
                row_ind_table = find(this_sample_ind_table);
                row_ind_pca = find(this_sample_ind_pca);
                pc_scores(row_ind_table(is_match), 1:4) = obj.geochem_pca.scores(row_ind_pca(match_ind(is_match)), 1:4);
                pc_scores(row_ind_table(is_match), 5:8) = obj.superpixel_pca.scores(row_ind_pca(match_ind(is_match)), 1:4);
            end

            % set the pc scores that are still 0 to nan
            pc_scores(pc_scores == 0) = NaN;

            % and make it a table
            pc_table = array2table(pc_scores, 'VariableNames', {'Geochem PC1', 'Geochem PC2',...
                'Geochem PC3', 'Geochem PC4', 'Image PC1', 'Image PC2', 'Image PC3', 'Image PC4'});
         
            % all good, now just make the superpixel table and strat height
            % table, and paste everything together
            im_stat_tab = array2table(super_stats, 'VariableNames', obj.stat_titles);     
            strat_height_tab =  array2table(strat_heights, 'VariableNames', {'stratigraphic height'});
            obj.master_table = [mast_tab(:,1:2), strat_height_tab, mast_tab(:,3:end), pc_table, im_stat_tab];
        end
        %%
        % GUI type plot in another file to go around and look at superpixel
        % stats. kept in separate file
        explore_plot(obj, stat_num, color_scheme)
        %%
        % GUI interface for exploring crossplots
        function crossplot_gui(obj)
            
            % extract the variable names for easy access, omitting the
            % first two
            var_names = fieldnames(obj.master_table);
            var_names = var_names(3:end);

            % this plot will have user control over both the x and y axes,
            % so we need to set up two ui sliders
            slide_step = 1/(numel(var_names)-1);
            f = figure('Visible','off');
            % and create axes where all this will go
            ax = axes('Position',[0.1 0.3 0.7 0.6]);
            px = uipanel(f,'Position',[0.1 0.1 0.2 0.1]);
            cx = uicontrol(px, 'Style', 'slider');
            py = uipanel(f,'Position',[0.4 0.1 0.2 0.1]);
            cy = uicontrol(py, 'Style', 'slider');
            
            % and axes for filters if we need them
            x_filt_ax = axes('Position',[0.65 0.05 .13 .13], 'Visible','off');
            y_filt_ax = axes('Position',[0.8 0.05 .13 .13], 'Visible','off');

            % Set up crosshairs on each axis at the edges
            gobj(1,1) = xline(ax,min(xlim(ax)), 'k-');
            gobj(1,2) = yline(ax,min(ylim(ax)), 'k-');
        
            % but don't have them show up in the legend
            set(get(get(gobj(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(gobj(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            % set some of the finer variables on those
            cx.Min = 1;
            cx.Max = numel(var_names);
            cx.Value = 1;
            cx.SliderStep = [slide_step slide_step];
            cx.Callback = {@plot_vars, obj, obj.master_table(:,3:end), ax, var_names, x_filt_ax, y_filt_ax, gobj};
            cy.Min = 1;
            cy.Max = numel(var_names);
            cy.Value = 1;
            cy.SliderStep = [slide_step slide_step];
            cy.Callback = {@plot_vars, obj, obj.master_table(:,3:end), ax, var_names, x_filt_ax, y_filt_ax, gobj};
            f.Visible='on';
            
            % and some additional ui boxes just for text labels of the
            % sliders
            ptext1 = uipanel(f,'Position',[0.1 0.02 0.2 0.1]);
            uicontrol(ptext1, 'Style', 'text', 'String', 'x variable');
            ptext2 = uipanel(f,'Position',[0.4 0.02 0.2 0.1]);
            uicontrol(ptext2, 'Style', 'text', 'String', 'y variable');

            % set up blank figure to receive the image of the clicked point
            im_fig = figure('Visible','off');
            
            % in order to click a point and see the the accompanying image,
            % need to windowbutton stuff
            % Assign windowbuttonmotion fcn on axis #1
            set(ax.Parent,'windowbuttonmotionfcn', {@mouseMove, ax, gobj, obj, cx, cy, obj.master_table(:,3:end), im_fig});
            
            % Assign mouse button functions to start/stop tracking
            WindowButtonMotionFcnInput = {@mouseMove, ax, gobj, obj, cx, cy, obj.master_table(:,3:end), im_fig};
            set(ax.Parent,'windowbuttondownfcn', {@startStopMouseMove, WindowButtonMotionFcnInput})

            % no matter which slider was moved, we need to rescatter the
            % data, so just one function to respond to movements
            function plot_vars(~, ~, obj, plotting_table, ax, var_names, x_filt_ax, y_filt_ax, gobj)
                % extract the index values for the variables to plot
                x_ind = round(cx.Value);
                y_ind = round(cy.Value);

                % and set up a logical to mark all non-0 and non-nan values
                % in these variables
                x_good_data = and(~isnan(table2array(plotting_table(:,x_ind))), table2array(plotting_table(:,x_ind)) ~= 0);
                y_good_data = and(~isnan(table2array(plotting_table(:,y_ind))), table2array(plotting_table(:,y_ind)) ~= 0);
                all_good_data = and(x_good_data, y_good_data);
                x_plot = table2array(plotting_table(all_good_data,x_ind));
                y_plot = table2array(plotting_table(all_good_data,y_ind));

                % Set up crosshairs on each axis at the edges
                gobj(1,1) = xline(ax,min(xlim(ax)), 'k-');
                gobj(1,2) = yline(ax,min(ylim(ax)), 'k-');
            
                % but don't have them show up in the legend
                set(get(get(gobj(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                set(get(get(gobj(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

                % now prep some stuff and loop through phase codes
                colors = get(gca, 'ColorOrder');
                filled_markers = {"o", "diamond", "square", "^", "v", ">", "<"};
                % and we're ready to go, clear the axes first
                all_scatters = findobj('Type', 'Scatter');
                if ~isempty(all_scatters)
                    delete(all_scatters)
                end
                dead_lines = findobj('Label','');
                if ~isempty(dead_lines)
                    delete(dead_lines)
                end
                cla(x_filt_ax, 'reset');
                cla(y_filt_ax, 'reset');
                x_filt_ax.Visible = 'off';
                y_filt_ax.Visible = 'off';
                for i = 1:size(obj.phase_codes,1)
                    % first pick a symbol based on where we are in looping
                    % through the color order
                    this_marker = filled_markers{ceil(i/size(colors,1))};
    
                    % and get the inds of the stats that plot as part of this
                    % phase
                    these_stats = find(strcmpi(obj.master_table.phase(all_good_data), obj.phase_codes{i,1}));
    
                    % ready to scatter, first getting an ind for the color row
                    % that accounts for the fact that i may at some point be
                    % larger than the color order
                    if i > size(colors,1)
                        color_ind = i+1 - floor(i/size(colors,1))*size(colors,1);
                    else
                        color_ind = i;
                    end
                    % and only scatter if we have data for this phase
                    if numel(these_stats) > 0
                        hold(ax, "on")
                        scatter(ax, x_plot(these_stats), y_plot(these_stats), 'filled', this_marker, 'MarkerFaceColor', colors(color_ind,:), 'DisplayName', obj.phase_codes{i,2});
                        hold(ax, "off")
                    end
                end

                legend(ax, 'Position', [0.8, 0.3, 0.15, 0.6])
                % sort out stuff with underscores in some variable names
                for_x_label = var_names{x_ind};
                for_x_label = strrep(for_x_label, '_', ' ');
                for_y_label = var_names{y_ind};
                for_y_label = strrep(for_y_label, '_', ' ');
                xlabel(ax, for_x_label)
                ylabel(ax, for_y_label)

                % and if we have filter in the name of a variable, show the
                % filter
                if contains(for_x_label, 'filter')
                    x_filt_ax.Visible = 'on';
                    filter_split = split(for_x_label, 'filter ');
                    filter_num = str2num(filter_split{2}(1));
                    imagesc(x_filt_ax, obj.filter_bank(:,:,filter_num))
                    axis(x_filt_ax, 'image');
                    grid(x_filt_ax, 'off');
                    set(x_filt_ax,'YTickLabel',[]);
                    set(x_filt_ax,'XTickLabel',[]);
                    xlabel(x_filt_ax, 'x filter')
                end

                 if contains(for_y_label, 'filter')
                    y_filt_ax.Visible = 'on';
                    filter_split = split(for_y_label, 'filter ');
                    filter_num = str2num(filter_split{2}(1));
                    imagesc(y_filt_ax, obj.filter_bank(:,:,filter_num))
                    axis(y_filt_ax, 'image');
                    grid(y_filt_ax, 'off');
                    set(y_filt_ax,'YTickLabel',[]);
                    set(y_filt_ax,'XTickLabel',[]);
                    xlabel(y_filt_ax, 'y filter')
                end
            end

            % a function to respond to movement on the scatterplot with
            % crosshairs
            function mouseMove(~, ~, ax, gobj, ~, ~, ~, ~, ~)
                % Responds to mouse movement in the scatterplot axis
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
            
            % a function to handle displaying an image when a point is
            % clicked on the scatter plot
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
                    % should display the image and highlight the clicked
                    % point
        
                    % start by grabbing the closest point. 
                    % the final crosshair values are stored in the constant line
                    % objects in the 3rd cell of the Window..Input array
                    x_clicked = get(WindowButtonMotionFcnInput{3}(1), 'Value');
                    y_clicked = get(WindowButtonMotionFcnInput{3}(2), 'Value');
        
                    % and we've passed the petro_superpix object in as the 4th cell
                    % of the input arry, so we can look at the x-y locations to see
                    % which one we're closest to
                    obj = WindowButtonMotionFcnInput{4};
                    cx = WindowButtonMotionFcnInput{5};
                    cy = WindowButtonMotionFcnInput{6};
                    plotting_table = WindowButtonMotionFcnInput{7};
                    % extract the index values for the variables to plot
                    x_ind = round(cx.Value);
                    y_ind = round(cy.Value);
                    these_data = [table2array(plotting_table(:,x_ind)), table2array(plotting_table(:,y_ind))];
                    nearest_point = knnsearch(these_data,[x_clicked,y_clicked],'K',1);
       
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
                    colors = get(gca, 'ColorOrder');
                    title_string = ['smg_' obj.master_table.sample_number{nearest_point}];
                    this_phase_ind = find(strcmpi(obj.phase_codes(:,1), obj.master_table.phase{nearest_point}(1)));
                    hold(WindowButtonMotionFcnInput{2}, 'on')
                    if this_phase_ind > size(colors,1)
                        color_ind = this_phase_ind+1 - floor(this_phase_ind/size(colors,1))*size(colors,1);
                    else
                        color_ind = this_phase_ind;
                    end
                    scatter(WindowButtonMotionFcnInput{2}, table2array(plotting_table(nearest_point,x_ind)), table2array(plotting_table(nearest_point,y_ind)),500,colors(color_ind,:),...
                        'filled','pentagram','MarkerEdgeColor',[.1,.1,.1], 'DisplayName', strrep(title_string,'_', ' '));
                    hold(WindowButtonMotionFcnInput{2}, 'off')
                    legend('Position', [0.8, 0.3, 0.15, 0.6])

                    % now we're clear to pop up another figure with the
                    % image of the rock that was clicked, and an outline of
                    % the superpixel that the point was drilled from
                    sample_ind = strcmpi(title_string, obj.sample_names);
                    rgb_im = imread(fullfile(obj.low_res_path,[obj.sample_names{sample_ind} '.jpg']));
                    superpixel_inds = imresize(imread(fullfile(obj.superpixel_ind_path, [obj.sample_names{sample_ind}, '.tif'])),0.05,'nearest');
                    ind_col = superpixel_inds(:);

                    % need the scale ratio between the superpixel indices and the
                    % rgb image for display of the highest and lowest value
                    % superpix
                    scale_ratio = size(superpixel_inds,1)/size(rgb_im,1);

                    % also need to account for background
                    background_supers = find(strcmpi(obj.class_labels{sample_ind},'background'));
                    background_mask = ones(size(ind_col));
                    background_mask(ismember(ind_col, background_supers)) = 0;
                    background_mask = imresize(reshape(background_mask, size(superpixel_inds)), scale_ratio);

                    % and get an outline for the superpixel from this point
                    % in particular
                    this_super = ismember(superpixel_inds, obj.geochem_superpix_inds{nearest_point,3});
                    super_bounds = bwboundaries(this_super);

                    % alright, let's show the thing
                    pos1 = [0.14 0.05 0.8 0.7];
                    pos2 = [0.1 0.85 0.8 0.15];
                    this_fig = WindowButtonMotionFcnInput{8};
                    this_fig.Visible = 'on';
                    set(0, 'currentfigure', this_fig);
                    clf(this_fig,'reset')
                    sp1 = subplot('Position', pos1);
                    hold(sp1, "on")
                    imshow(rgb_im)
                    title(strrep(title_string,'_', ' '))
                    patch(super_bounds{1}(:,2)./scale_ratio,super_bounds{1}(:,1)./scale_ratio,'white','FaceColor', 'none', 'EdgeColor', [1, 1, 1],'LineWidth', 2);
                    drawnow
                    hold(sp1, "off")
        
                    % and add text to the the right with a categorical note
                    % about image orientation
                    sp2 = subplot('Position', pos2);
                    grid(sp2,'off')
                    set(sp2,'YTickLabel',[]);
                    set(sp2,'XTickLabel',[]);
                    set(sp2,'YTick',[])
                    set(sp2,'XTick',[])
                    xlim(sp2, [0.5 10])
                    ylim(sp2, [0.5 1.5])
                    if isnan(obj.im_bedding_angles(sample_ind))
                        text(1,1,'Image orientation not available.', 'FontSize', 18)
                    elseif obj.im_bedding_angles(sample_ind) < 15 || obj.im_bedding_angles(sample_ind) > 165
                        text(1,1,'Image is approximately along bedding plane.', 'FontSize', 18)
                    elseif obj.im_bedding_angles(sample_ind) > 75 && obj.im_bedding_angles(sample_ind) < 105
                        text(1,1,'Image is a cross-section wrt bedding.', 'FontSize', 18)
                    else
                        text(1,1,'Image is oblique wrt bedding.', 'FontSize', 18)
                    end
            end
        end
            
        end
        
        %%
        % GUE interface for explorting strat sections of variables with
        % bootstrapped mean and standard deviations
        function strat_gui(obj)
            % extract the variable names for easy access, omitting the
            % first two
            var_names = fieldnames(obj.master_table);
            var_names = var_names(3:end);

            % this plot will have user control over both the x and y axes,
            % so we need to set up two ui sliders
            slide_step = 1/(numel(var_names)-1);
            f = figure('Visible','off');
            % and create axes where all this will go
            ax = axes('Position',[0.1 0.3 0.7 0.6]);
            px = uipanel(f,'Position',[0.1 0.1 0.2 0.1]);
            cx = uicontrol(px, 'Style', 'slider');
            
            % and axes for filters if we need them
            filt_ax = axes('Position',[0.65 0.05 .13 .13], 'Visible','off');
  
            % set some of the finer variables on the ui slider
            cx.Min = 1;
            cx.Max = numel(var_names);
            cx.Value = 1;
            cx.SliderStep = [slide_step slide_step];
            cx.Callback = {@plot_sec, obj, obj.master_table(:,3:end), ax, var_names, filt_ax};
            f.Visible='on';
            
            % and an additional ui box just for a slicer label
            ptext1 = uipanel(f,'Position',[0.1 0.02 0.2 0.1]);
            uicontrol(ptext1, 'Style', 'text', 'String', 'x variable');

            % no matter which slider was moved, we need to rescatter the
            % data, so just one function to respond to movements
            function plot_sec(~, ~, obj, plotting_table, ax, var_names, filt_ax)
                % extract the index values for the variable to plot
                var_ind = round(cx.Value);

                % and set up a logical to mark all non-0 and non-nan values
                % in these variables
                good_data = and(~isnan(table2array(plotting_table(:,var_ind))), table2array(plotting_table(:,var_ind)) ~= 0);
                data_plot = table2array(plotting_table(good_data,var_ind));
                strat_heights = table2array(plotting_table(good_data,1));
                
                % now, we'll need to bootstrap resample in order to plot a
                % reasonable mean and standard deviation of the data
                % start by getting a strat height distance metric, at this
                % point arbitrarility scaled by 10m.
                strat_dists = pdist2(strat_heights, strat_heights)./10;
                strat_dists(eye(size(strat_dists))==1) = nan;

                % and get a proximity value that is the reciprocal of the
                % sum of squared scaled strat dists plus 1
                prox_weights = sum(1./(1 + strat_dists.^2), 'omitnan');

                % turn the pros values into probabilities of sampling,
                % still replicating math from Akshay's paper
                prob_vals = 1 ./ ((prox_weights .* median(5 ./ prox_weights)) + 1);

                % and random sample with replacement, not worrying about
                % uncertainties for now
                sample_ids = [1:numel(data_plot)]';
                random_samples = randsample(sample_ids, 10000, true, prob_vals);
                sampled_data = data_plot(random_samples);
                sampled_strat_heights = strat_heights(random_samples);

                % all that's left to to is discretize based upon the strat
                % heights and get the stats for each bin. Discretizing
                % every 5m right now
                strat_edges = linspace(min(strat_heights), max(strat_heights), ceil(range(strat_heights)/5));
                strat_binned = discretize(sampled_strat_heights, strat_edges);

                % may have empty bins, so use findgroups to work around
                edited_bins = findgroups(strat_binned);
                
                % we want to know the empty bins at the end still, so we'll
                % make a vector of nans for the stats and just fill in the
                % splitapply stats around them
                bin_means = nan(1,numel(strat_edges)-1);
                bin_10pct = nan(1,numel(strat_edges)-1);
                bin_90pct = nan(1,numel(strat_edges)-1);
                missing_bins = setdiff([1:numel(strat_edges)-1], strat_binned);
                filled_bins = setdiff([1:numel(strat_edges)-1], missing_bins);

                bin_means(filled_bins) = splitapply(@mean, sampled_data, edited_bins);
                bin_10pct(filled_bins) = splitapply(@(x) prctile(x,10), sampled_data, edited_bins);
                bin_90pct(filled_bins) = splitapply(@(x) prctile(x,90), sampled_data, edited_bins);
                
                % before plotting, just confirm the corresponding strat
                % heights for the bins
                strat_bin_centers = (strat_edges(2:end)+strat_edges(1:end-1))./2;
               
                % when plotting, we'll want to collapse the percentile
                % window at missing bins. We'll collapse it to the center,
                % so we need to get interpolated values of the means in the
                % missing windows
                bin_means = fillmissing(bin_means,'linear');
                bin_10pct(missing_bins) = bin_means(missing_bins);
                bin_90pct(missing_bins) = bin_means(missing_bins);

                % now we're ready to actually plot
                % start by clearing the axes
                cla(ax,"reset")
                
                % make a patch in the back using the 10 and 90 percentile
                % windows
                colors = get(gca, 'ColorOrder');
                patch(ax,[bin_10pct, flip(bin_90pct)], [strat_bin_centers, flip(strat_bin_centers)], colors(2,:), 'EdgeAlpha', 0, 'FaceAlpha', 0.3)
                hold(ax, 'on')
                plot(ax, bin_means, strat_bin_centers)
                
                % just finish up with some axis labeling
                for_x_label = var_names{var_ind};
                xlabel(ax, for_x_label)
                ylabel(ax, var_names{1})

                % and if we have filter in the name of a variable, show the
                % filter
                if contains(for_x_label, 'filter')
                    filt_ax.Visible = 'on';
                    filter_split = split(for_x_label, 'filter ');
                    filter_num = str2num(filter_split{2}(1));
                    imagesc(filt_ax, obj.filter_bank(:,:,filter_num))
                    axis(filt_ax, 'image');
                    grid(filt_ax, 'off');
                    set(filt_ax,'YTickLabel',[]);
                    set(filt_ax,'XTickLabel',[]);
                    xlabel(filt_ax, 'x filter')
                end
            end
        end
        %%
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
       
        %%
        % not even using the images, can make histograms of all the
        % different phases
        function phase_histos(obj,chem_variable)
            % first decide which chemical variable we're making histograms
            % of, and choose the appropriate column in the geochem table
            chem_headers = fieldnames(obj.all_geochem_data);
            if strcmpi(chem_variable,'carb')
                % user wants carbon isotopes
                data_col = find(strcmpi('delta_13C_mean',chem_headers));
                for_title = '\delta^{13}C';
            elseif strcmpi(chem_variable,'ox')
                % user wants oxygen isotopes
                data_col = find(strcmpi('delta_18O_mean',chem_headers));
                for_title = '\delta^{18}O';
            else 
                % user wants an element concentration
                data_col = find(strcmpi(chem_variable,chem_headers));
                for_title = strrep(chem_variable,'_','/');
            end
            
            % phases are in the second column of that cell array and they
            % might have quailfiers on them, so we can get rid of those for
            % now so we just have the phase letters
            only_phases = cellfun(@(x) x(1), obj.all_geochem_data.phase);

            % and then just gather the unique phases and iterate through
            unique_phases = unique(only_phases);

            % before plotting, we'll curate the data a little bit to remove
            % zeros, nans, and outliers
            chosen_data = table2array(obj.all_geochem_data(:,data_col));
            is_zero = chosen_data == 0;
            is_nan = isnan(chosen_data);
            chosen_data = chosen_data(~or(is_nan,is_zero));
            only_phases = only_phases(~or(is_nan,is_zero));

            [~,is_outlier] = rmoutliers(chosen_data); 
            chosen_data = chosen_data(~is_outlier);
            only_phases = only_phases(~is_outlier);
            figure()
            for i = 1:numel(unique_phases)
                these_samples = strcmp(unique_phases(i), cellstr(only_phases));
                subplot(numel(unique_phases),1,i)
                histogram(chosen_data(these_samples))
                phase_ind = strcmpi(unique_phases(i), obj.phase_codes(:,1));
                xlabel([obj.phase_codes{phase_ind,2} ' (n = ' num2str(sum(these_samples)) ')'])
                xlim([min(chosen_data)-0.2*range(chosen_data), max(chosen_data)+0.2*range(chosen_data)])
            end
            sgtitle(for_title)

        end
        
        %%
        % general function to call for plotting data on a map,color coded
        % by a variable value
        function map_plot(obj, data_mat, title_string, color_scheme)
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
            buffer_size = 50;
            data_xmin = min(data_mat(:,1)) - buffer_size;
            data_xmax = max(data_mat(:,1)) + buffer_size;
            data_ymin = min(data_mat(:,2)) - buffer_size;
            data_ymax = max(data_mat(:,2)) + buffer_size;
        
            % now just need to get where those numbers exist in pixel space on the
            % dem
            imshow_start_row = round((max(utm_ylims) - data_ymax)*dem_yscale);
            imshow_end_row = round((max(utm_ylims) - data_ymin)*dem_yscale);
            imshow_start_col = round((data_xmin - min(utm_xlims))*dem_xscale);
            imshow_end_col = round((data_xmax - min(utm_xlims))*dem_xscale);
        
            % in the dem, background pixels are set to -32767, so reset those to
            % nan
            dem_im(dem_im == -32767) = NaN;
            
            % and grab the part to use as the basemap and rescale it
            to_show = imresize(dem_im(imshow_start_row:imshow_end_row, imshow_start_col:imshow_end_col), 1/mean([dem_xscale,dem_yscale]),'bicubic');
            % last thing is to set up an alpha channel to block NaNs
            dem_alpha = ones(size(to_show));
            dem_alpha(isnan(to_show)) = 0;
            
            
            figure('units','normalized','outerposition',[0 0 0.9 1])
            title(title_string)
            grid off
            axis off
            
            % start by showing the image on custom axes to mirror aspect ratio of the image
            ax1 = axes();
            title(title_string)
            imagesc(ax1, to_show, 'AlphaData',dem_alpha);
            set(ax1,'color',0*[1 1 1])
            axis image
            hold on
            % make the scatter 
            ax2 = axes('position', ax1.Position);
            axis image
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
            
            % now, when plotting the actual data, we might have been given
            % a field facies for each datapoint that we'll have to account
            % for. Otherwise we just straigth up plot on the map
            if size(data_mat,2) == 4 
                % we'll want the field facies loaded in for later
                load('field_facies.mat');
                
                % we'll loop through the field facies and scatter each one with a
                % different marker
                filled_markers = {"o", "square", "diamond", "^", "v", ">", "<"};
                for i = 1:numel(field_facies)
                    these_lithos = data_mat(:,4) == i;
                    stat_frac = (data_mat(these_lithos,3) - min(data_mat(:,3)))/range(data_mat(:,3));
                    cmap_ind = round(size(cmap,1)*stat_frac);
                    cmap_ind(cmap_ind == 0) = 1;
                    scatter(ax2,data_mat(these_lithos,1), data_mat(these_lithos,2),40,cmap(cmap_ind, :),'filled', ...
                        filled_markers{i},'MarkerEdgeColor',[.1,.1,.1], 'DisplayName',field_facies{i});
                end
            else
                    stat_fracs = (data_mat(:,3) - min(data_mat(:,3)))/range(data_mat(:,3));
                    cmap_inds = round(size(cmap,1)*stat_fracs);
                    cmap_inds(cmap_inds == 0) = 1;
                    scatter(ax2,data_mat(:,1), data_mat(:,2),40,cmap(cmap_inds, :),'filled', ...
                        'MarkerEdgeColor',[.1,.1,.1], 'DisplayName','geochem data');
            end
            hold(ax2,"off");
            % refine some stuff on the second axes for visibility and geography
            ax2.Visible = 'on';
            box(ax2,'off')
            set(ax2,'color','none')
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
            colormap(ax1, gray_cmap(1:end,:));
            
            % make sure we have a legend
            legend('Location','northwest')

            % and of course we actually need to put the colorbar on
            % separate axes so it doesn't interfere with the plot
            % overlaying the map
            ax3 = axes;
            colormap(ax3, cmap);
            c = colorbar(ax3);
            c.Ticks = linspace(0,1,9);
            c.TickLabels = num2cell(linspace(min(data_mat(:,3)), max(data_mat(:,3)), 9));
            ax3.Visible = 'off';
        end
        %% 
        % display geochemical data in map view
        function geochem_map(obj, phase, chem_variable, color_scheme)
            % first decide which chemical variable we're making a map of,
            % and load the appropriate file
            chem_headers = fieldnames(obj.all_geochem_data);
            if strcmpi(chem_variable,'carb')
                % user wants carbon isotopes
                data_col = find(strcmpi('delta_13C_mean',chem_headers));
                for_title = '\delta^{13}C';
            elseif strcmpi(chem_variable,'ox')
                % user wants oxygen isotopes
                data_col = find(strcmpi('delta_18O_mean',chem_headers));
                for_title = '\delta^{18}O';
            else 
                % user wants an element concentration
                data_col = find(strcmpi(chem_variable,chem_headers));
                for_title = strrep(chem_variable,'_','/');
            end

            % grab the just sample values for the variable of choice
            chosen_data = table2array(obj.all_geochem_data(:,data_col));

            % okay, so the user would either select a specific phase to
            % display or ask to display the bulk value
            % handle those in separate paths
            sample_nums = obj.all_geochem_data.sample_number;
            if strcmpi(phase, 'bulk')
                % user just wants a bulk value, so we'll average for every
                % sample 
                % first get the unique sample names present
                unique_samples = unique(obj.all_geochem_data.sample_number);

                % just for loop through the sample labels and take the
                % mean of all corresponding chemical values
                for_map = zeros(numel(unique_samples),1);
                for i = 1:numel(unique_samples)
                    for_map(i) = mean(chosen_data(strcmpi(unique_samples{i}, sample_nums)),'omitnan');
                end
            elseif sum(strcmpi(phase, obj.phase_codes(:,1))) == 1
                % user wants a specific phase that matches the codes we
                % have in the object
                % first narrow down to samples that have this phase, only
                % checking the first character in case there are modifying
                % numbers
                first_char = cellfun(@(x) x(1), obj.all_geochem_data.phase, 'UniformOutput',false);
                this_phase = strcmpi(phase, first_char);
                unique_samples = unique(sample_nums(this_phase,1));

                % just for loop through the sample labels and take the
                % mean of all corresponding chemical values for the right
                % phase             
                for_map = zeros(numel(unique_samples),1);
                for i = 1:numel(unique_samples)
                    for_map(i) = mean(chosen_data(and(strcmpi(unique_samples{i}, sample_nums), this_phase)),'omitnan');
                end
            else
                % we don't have a matching geochemical phase, so return an
                % error
                error('not a valid geochemical phase')
            end

            % it's possible that the only value for a sample was nan, so
            % get rid of that 
            is_nan = isnan(for_map);
            for_map = for_map(~is_nan);
            unique_samples = unique_samples(~is_nan);

            % also possible that the value was 0 and we want to get rid of
            % those
            is_zero = for_map == 0;
            for_map = for_map(~is_zero);
            unique_samples = unique_samples(~is_zero);
            
            % now we're ready to plot on the map
            % first get the lat long coordinates of the points we just
            % extracted
            all_split_samps = cellfun(@(x) strsplit(x,'_'), obj.sample_names', 'UniformOutput',false);
            all_split_samps = vertcat(all_split_samps{:});
            all_samp_nums = all_split_samps(:,2);

            % now have to get the indices of the unique_samples in the
            % overall sample names of the object
            [~, coord_inds] = ismember(lower(unique_samples), all_samp_nums);
  
            % last thing before plotting is to set up the data matrix with
            % lat, long, the geochemical data, and field facies
            data_mat = [obj.utm_xyz(coord_inds,1), obj.utm_xyz(coord_inds,2), for_map, obj.field_lithos(coord_inds)'];
            % and now call the map plotting function

            obj.map_plot(data_mat, for_title, color_scheme);
        end
        
        %%
        % take in the geochem inds from a petro im and set the
        % geochem_superpix_inds
        function set_geochem_super_inds(obj)
            % this property is going to reflect the geochem data table. But
            % it'll be a cell array for simplicity
            clicked_inds = cell(size(obj.all_geochem_data,1),3);
            clicked_inds(:,1) = table2array(obj.all_geochem_data(:,1));
            clicked_inds(:,2) = table2array(obj.all_geochem_data(:,2));

            % now we just need to make a list of unique sample numbers and
            % loop through them
            unique_samples = unique(clicked_inds(:,1));
            for i = 1:numel(unique_samples)
                % assemble a path to read in the petro_im instance for this
                % sample
                petro_path = [obj.petro_ims_path '/smg_' unique_samples{i} '.mat'];
                load(petro_path);
                % just reset the variable name of the loaded petro_im so we
                % don't have to keep using eval lines
                eval(['this_pet = smg_' unique_samples{i} ';'])
                eval(['clear  smg_' unique_samples{i}]);

                % now we're ready to fill out clicked inds
                n_superpix_ind = find(obj.n_superpixels == this_pet.n_superpixels);
                this_samp = strcmpi(unique_samples{i}, clicked_inds(:,1));
                % gotta match the phases
                [is_match,matched_inds] = ismember(clicked_inds(this_samp,2), this_pet.geochem_superpix_inds(:,1));
                % a few may not have matched because a 1 is missing in the
                % clicked inds, so try to fill that in
                if sum(is_match) < sum(this_samp)
                    add_one = repmat({'1'},numel(this_pet.geochem_superpix_inds(:,1)),1);
                    [retry_match, retry_inds] = ismember(clicked_inds(this_samp,2), strcat(this_pet.geochem_superpix_inds(:,1), ' ', add_one));
                    if sum(retry_match) > 0
                        matched_inds(retry_inds ~= 0) = retry_inds(retry_inds ~=0);
                        is_match(retry_inds ~= 0) = true;
                    end
                end
                
                % all set, just place the inds
                all_these_super_inds = cell2mat(vertcat(this_pet.geochem_superpix_inds(:,2)));
                clicked_inds(this_samp,3) = num2cell(all_these_super_inds(matched_inds(is_match), n_superpix_ind));
            end

            % export the property
            obj.geochem_superpix_inds = clicked_inds;
        end
   
        %% 
        % perform PCA on geochem data
        function set_geochem_pca(obj)
            % ugh, so many different data types floating around in this
            % class, but for this, I think it makes sense to make it a
            % structure with easily accessible fields
            obj.geochem_pca = struct;

            % and we'll actually only perform the pca on the variables that
            % ran well and mostly have data
            all_variables = fieldnames(obj.all_geochem_data);
            not_data = {'sample_number', 'phase', 'line', 'delta_13C_std', 'delta_18O_std', 'date_run',...
                'im_location_col', 'im_location_row'};
            not_data_inds = ismember(all_variables, not_data);
            just_data = table2array(obj.all_geochem_data(:,~not_data_inds(1:end-3)));
            frac_nan = sum(isnan(just_data))./size(just_data,1);

            % so now we can hold out variables based on not having enough
            % data. We'll use 85% as the threshold. Also note that we'll
            % only consider trace element variables normalized to calcium
            dont_use = frac_nan > 0.1;
            norm_mg = contains(all_variables(~not_data_inds), '_Mg');
            norm_camg = contains(all_variables(~not_data_inds), '_CaMg');
            final_not_include = any([dont_use', norm_mg(1:end-3), norm_camg(1:end-3)], 2);

            good_var_names = all_variables(~not_data_inds);
            good_var_names = good_var_names(1:end-3);
            good_var_names = good_var_names(~final_not_include);
            geochem_subset = table2array(obj.all_geochem_data(:,find(ismember(all_variables,good_var_names))));

            % and samples that haven't been run for trace elements yet have
            % zeros there, so let's leave them out. Also leave out ones
            % that ran poorly and have an Nan value
            not_run = any(geochem_subset == 0,2);
            bad_run = any(isnan(geochem_subset),2);
            leave_out = any([not_run, bad_run], 2);
            geochem_subset = geochem_subset(~leave_out,:);

            % be sure to normalize each column to be zero centered with a
            % standard deviation of 1.
            geochem_zero_center = geochem_subset - mean(geochem_subset);
            geochem_normalized = geochem_zero_center./std(geochem_zero_center);

            % ready to get pcas
            [coeff, score, latent, ~, explained] = pca(geochem_normalized);

            % fill out the structure to export
            obj.geochem_pca.scores = score;
            obj.geochem_pca.loadings = coeff;
            obj.geochem_pca.pct_explained = explained;
            obj.geochem_pca.latent = latent;
            obj.geochem_pca.input_vars = good_var_names;
            obj.geochem_pca.samp_nums = obj.all_geochem_data.sample_number(~leave_out);
            obj.geochem_pca.phases = obj.all_geochem_data.phase(~leave_out);

        end

        %% 
        % Perform a PCA on superpixel stats
        function set_superpixel_pca(obj, n_considered, stat_subset)
            % this will go down in a similar fashion to the geochem pca,
            % but it's unrealistic to include all superpixels in the pca,
            % so we'll include those that represent a drilled part for
            % geochem plus a random sampling of the others to get us up
            % to the number input by the user.

            % so start with sampling just the superpixels that were sampled
            % for geochem. Gonna need to do some name comparison
            split_names = cellfun(@(x) strsplit(x,'_'), obj.sample_names, 'UniformOutput',false);
            split_names_cell = vertcat(split_names{:});

            % we're gonna do this based upon the names in the geochem pca,
            % so if that hasn't been done, do it now
            if isempty(obj.geochem_pca)
                obj.set_geochem_pca;
            end

            % get the list of sample names and phases from the geochem pca
            % to match superpix inds to 
            samped_samps = [obj.geochem_pca.samp_nums, obj.geochem_pca.phases];
            [~, match_ind] = ismember(strcat(samped_samps(:,1),' ',samped_samps(:,2)), strcat(obj.geochem_superpix_inds(:,1), ' ', obj.geochem_superpix_inds(:,2)));
            super_inds = cell2mat(obj.geochem_superpix_inds(match_ind,3));

            % now we're ready to go through and sample the superpixels to
            % get our matrix of superpix that were sampled for geochem,
            % along the way marking that they were already sampled or
            % background and don't need to be considered for remaining
            % random sample.
            total_supers = sum(cell2mat(cellfun(@(x) size(x, 1), obj.superpix_stats, 'UniformOutput',false)));
            for_pca = zeros(n_considered, numel(obj.stat_titles));
            dont_sample = [];
            count = 0;
            for i = 1:numel(split_names)
                has_geochem = strcmpi(samped_samps(:,1), split_names_cell{i,2});
                if sum(has_geochem) > 0
                    for_pca(find(has_geochem), :) = obj.superpix_stats{i}(super_inds(has_geochem),:);
                    dont_sample = [dont_sample; super_inds(has_geochem) + count];
                end
                dont_sample = [dont_sample; find(strcmpi(obj.class_labels{i}, 'background')) + count];
                count = count + size(obj.superpix_stats{i},1);
            end

            % and now fill in the rest that we need to randomly sample
            population = 1:total_supers;
            population(dont_sample) = [];
            to_sample = randsample(population,n_considered-numel(super_inds));
            count = 0;
            for i = 1:numel(split_names)
                sample_here = to_sample(and(to_sample > count, to_sample < (count+ size(obj.superpix_stats{i},1)+1)));
                start_row = find(all(for_pca == 0,2), 1,"first");
                for_pca(start_row:start_row+numel(sample_here)-1, :) = obj.superpix_stats{i}(sample_here-count,:);
                count = count + size(obj.superpix_stats{i},1);
            end

            % be sure to normalize each column to be zero centered with a
            % standard deviation of 1.
            superpix_zero_center = for_pca - mean(for_pca, 'omitnan');
            superpix_normalized = superpix_zero_center./std(superpix_zero_center, 'omitnan');
            % ready to get pcas
            superpix_normalized = fillmissing(superpix_normalized(:,stat_subset),'nearest');
            [coeff, score, latent, ~, explained] = pca(superpix_normalized);

            % fill out the structure to export
            obj.superpixel_pca.scores = score;
            obj.superpixel_pca.loadings = coeff;
            obj.superpixel_pca.pct_explained = explained;
            obj.superpixel_pca.latent = latent;
            obj.superpixel_pca.input_vars = obj.stat_titles(stat_subset);
        end

        %%
        % show a dashboard of PCA results
        function pca_dashboard(obj)
            % first make sure we have geochem and superpixel pcas to
            % display. if not, make them.
            if isempty(obj.geochem_pca)
                obj.set_geochem_pca;
            end
                
            if isempty(obj.superpixel_pca)
                obj.set_superpixel_pca(10000, [3, 4, 5, 6, 7, 8, 9, 10, 23, 26]);
            end

            % okay, we're going to make a 2x3 set of subplots. Top row will
            % be two crossplots and PC loadings for geochem. Bottom row
            % will be two crossplots and PC loadings for superpix. We'll do
            % the crossplots phase-wise, so we'll loop through and do all
            % those right away
            figure()
            % before we loop throught the phase codes, we can just plot all
            % the superpixel data in those crossplots as gray as a
            % background of all the data
            subplot(2,3,4)
            scatter(obj.superpixel_pca.scores(:,1), obj.superpixel_pca.scores(:,2), 20, [0.4, 0.4, 0.4]);
            subplot(2,3,5)
            scatter(obj.superpixel_pca.scores(:,3), obj.superpixel_pca.scores(:,4), 20, [0.4, 0.4, 0.4]);
            
            % now prep some stuff and loop through phase codes
            colors = get(gca, 'ColorOrder');
            filled_markers = {"o", "square", "diamond", "^", "v", ">", "<"};
            for i = 1:size(obj.phase_codes,1)
                % first pick a symbol based on where we are in looping
                % through the color order
                this_marker = filled_markers{ceil(i/size(colors,1))};

                % and get the inds of the scores that plot as part of this
                % phase code
                these_scores = find(strcmpi(obj.geochem_pca.phases, obj.phase_codes{i,1}));

                % ready to scatter, first getting an ind for the color row
                % that accounts for the fact that i may at some point be
                % larger than the color order
                if i > size(colors,1)
                    color_ind = i+1 - floor(i/size(colors,1))*size(colors,1);
                else
                    color_ind = i;
                end
                % and only scatter if we have data for this phase
                if numel(these_scores) > 0
                    subplot(2,3,1)
                    hold on
                    scatter(obj.geochem_pca.scores(these_scores,1), obj.geochem_pca.scores(these_scores,2), 'filled', this_marker, 'MarkerFaceColor', colors(color_ind,:), 'DisplayName', obj.phase_codes{i,2});
                    hold off
                    subplot(2,3,2)
                    hold on
                    scatter(obj.geochem_pca.scores(these_scores,3), obj.geochem_pca.scores(these_scores,4), 'filled', this_marker, 'MarkerFaceColor', colors(color_ind,:));
                    hold off
                    subplot(2,3,4)
                    hold on
                    scatter(obj.superpixel_pca.scores(these_scores,1), obj.superpixel_pca.scores(these_scores,2), 'filled', this_marker, 'MarkerFaceColor', colors(color_ind,:));
                    hold off
                    subplot(2,3,5)
                    hold on
                    scatter(obj.superpixel_pca.scores(these_scores,3), obj.superpixel_pca.scores(these_scores,4), 'filled', this_marker, 'MarkerFaceColor', colors(color_ind,:));
                    hold off
                end
            end
            
            % set up legend and axis labels
            subplot(2,3,1)
            legend()
            xlabel(['Geochem PC1 (' num2str(round(obj.geochem_pca.pct_explained(1))) '%)'])
            ylabel(['Geochem PC2 (' num2str(round(obj.geochem_pca.pct_explained(2))) '%)'])
            subplot(2,3,2)
            xlabel(['Geochem PC3 (' num2str(round(obj.geochem_pca.pct_explained(3))) '%)'])
            ylabel(['Geochem PC4 (' num2str(round(obj.geochem_pca.pct_explained(4))) '%)'])
            subplot(2,3,4)
            xlabel(['Image PC1 (' num2str(round(obj.superpixel_pca.pct_explained(1))) '%)'])
            ylabel(['Image PC2 (' num2str(round(obj.superpixel_pca.pct_explained(2))) '%)'])
            subplot(2,3,5)
            xlabel(['Image PC3 (' num2str(round(obj.superpixel_pca.pct_explained(3))) '%)'])
            ylabel(['Image PC4 (' num2str(round(obj.superpixel_pca.pct_explained(4))) '%)'])
            
            % and just loading plots
            % first adjust some input variable names to succinctly label
            % the variables
            geochem_vars = obj.geochem_pca.input_vars;
            for i = 1:numel(geochem_vars)
                if contains(geochem_vars{i}, '_Ca')
                    geochem_vars{i} = strrep(geochem_vars{i}, '_', '/');
                elseif contains(geochem_vars{i}, 'delta_13')
                    geochem_vars{i} = '\delta^1^3C';
                elseif contains(geochem_vars{i}, 'delta_18')
                    geochem_vars{i} = '\delta^1^8O';
                end
            end

            superpix_vars = obj.superpixel_pca.input_vars;
            for i = 1:numel(superpix_vars)
                if contains(superpix_vars{i}, ' mean')
                    superpix_vars{i} = strrep(superpix_vars{i}, 'nm mean', '');
                elseif contains(superpix_vars{i}, '530nm GLCM ')
                    superpix_vars{i} = strrep(superpix_vars{i}, '530nm GLCM ', '');
                end
            end

            subplot(8,3,3)
            bar(obj.geochem_pca.loadings(:,1))
            set(gca,'xtick',[])
            chem_labels = text(1:numel(geochem_vars), ones(1,numel(geochem_vars)) .* 0.05, geochem_vars);
            set(chem_labels, 'Rotation', 90)
            ylabel('Geochem PC1', 'FontSize', 8)
            subplot(8,3,6)
            bar(obj.geochem_pca.loadings(:,2))
            set(gca,'xtick',[])
            ylabel('Geochem PC2', 'FontSize', 8)
            subplot(8,3,9)
            bar(obj.geochem_pca.loadings(:,3))
            set(gca,'xtick',[])
            ylabel('Geochem PC3', 'FontSize', 8)
            subplot(8,3,12)
            bar(obj.geochem_pca.loadings(:,4))
            set(gca,'xtick',[])
            ylabel('Geochem PC4', 'FontSize', 8)
            subplot(8,3,15)
            bar(obj.superpixel_pca.loadings(:,1))
            set(gca,'xtick',[])
            supers_labels = text(1:numel(superpix_vars), ones(1,numel(superpix_vars)) .* 0.05, superpix_vars, 'FontSize', 8);
            set(supers_labels, 'Rotation', 90)
            ylabel('Image PC1', 'FontSize', 8)
            subplot(8,3,18)
            bar(obj.superpixel_pca.loadings(:,2))
            set(gca,'xtick',[])
            ylabel('Image PC2', 'FontSize', 8)
            subplot(8,3,21)
            bar(obj.superpixel_pca.loadings(:,3))
            set(gca,'xtick',[])
            ylabel('Image PC3', 'FontSize', 8)
            subplot(8,3,24)
            bar(obj.superpixel_pca.loadings(:,4))
            set(gca,'xtick',[])
            ylabel('Image PC4', 'FontSize', 8)
        end
    end
end