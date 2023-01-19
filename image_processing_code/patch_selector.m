function[patches_table] = patch_selector(in_path, im_ext, patch_types, varargin)
% This function takes in an image data set and allows you to extract and
% label rectangular patches from each image. The labels and coordinates are
% stored in the output table for future use.
%
% INPUT 
% in_path: string with the full path to the images you would like to
% extract patches from.
% im_ext: String with the extension for the images
% patch types: cell array with strings for all the different categories of
% patches you might click in your images
% varargin: optional arguments to be put in as name pairs
%     savePath: path with file name where you would like the table to be
%     saved as an xlsx spreadhseet. Default is working directory.
%     scaleRatio: ratio of the size of your clicking images to the size of
%     the images you will be analyzing. In case you're clicking patches on
%     downsampled versions for efficiency. Default value is 1.
%     existingTable: path to existing table to be read in and added to.
%
% OUTPUT
% patches_table: table containing the image names, patch types, and
% coordinates for all clicked patches. This table also is saved to the
% input path or the current directory if no path is given.
% 
% Written by R.A. Manzuk
% 10/20/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end

    % make an input parser for varargin
    p = inputParser;
    
    % validation functions
    validPath = @(x) isfolder(fileparts(x));
    
    % default values
    defaultPath = [pwd '/patches.xlsx'];
    defaultScale = 1;
    defaultTable = cell2table(cell(0,6), 'VariableNames', {'image_name','patch_type','top_row','bot_row','left_col','right_col'});
    
    % possible parameters
    addParameter(p,'savePath',defaultPath);
    addParameter(p,'scaleRatio',defaultScale);
    addParameter(p,'existingTable',defaultTable);
    
    p.KeepUnmatched = true;
    % parse the inputs 
    parse(p,varargin{:});

    % now that we've parsed the inputs, we're ready to go through and click
    % some patches
    % extract the table from the parser to make it easier to call
    patches_table = p.Results.existingTable;

    % while loop we can get out of as a user
    clicking = true;
    while clicking
        % id the file the user would like to click
        disp('please select an image')
        [this_name] = uigetfile([in_path '/' im_ext]);
        
        % read in that image
        this_im = im2double(imread(fullfile(in_path,this_name)));
        
        % and show the image along with some info
        f1 = figure(1);
        subplot(1,3,[1,2])
        imshow(this_im)
        % if any patches for this image already exist, we should plot them
        [~,name_no_ext] = fileparts(this_name);
        existing_patch_inds = strcmp(name_no_ext,patches_table.image_name);
        if sum(existing_patch_inds) > 0 
            to_plot = find(existing_patch_inds == 1);
            hold on
            for i = 1:numel(to_plot)
                this_entry = patches_table(to_plot(i),:);
                to_plot_x = [this_entry.left_col,this_entry.right_col,this_entry.right_col,this_entry.left_col,this_entry.left_col];
                to_plot_y = [this_entry.top_row,this_entry.top_row,this_entry.bot_row,this_entry.bot_row,this_entry.top_row];
                plot(to_plot_x*p.Results.scaleRatio,to_plot_y*p.Results.scaleRatio,'r');
                % and we need to add the text for the category
                this_category = strcmp(patch_types,this_entry.patch_type);
                text(to_plot_x(1)*p.Results.scaleRatio + 20,to_plot_y(1)*p.Results.scaleRatio + 40,num2str(find(this_category ==1)),'Color','red','FontSize',14);
            end
            hold off
        end
        subplot(1,3,3)
        text(0.1,0.7,'press enter to move to the next image'); axis off
        text(0.1,0.65,'press "e" to end clicking session altogether')
        text(0.1,0.6,'watch the command window for instructions')
        text(0.1,0.55,'Patch Types:')
        for i = 1:numel(patch_types)
            text(0.1,0.55 - 0.05*i, [num2str(i), ' - ' patch_types{i}])
        end
        
        disp('Please click the upper left, followed by the lower right of the patch')
        % now we can ask the user to click the box(es)
        more_boxes = true;
        counter = 0;
        while more_boxes
            [x_i,y_i,button_i] = ginput(1);
            if isempty(x_i) % gets out of the while loop for this image because enter was pushed 
                more_boxes = false;
                close(f1); 
                break; 
            end % gets out of the while loop for this image
            if button_i == 101 % here the user pressed e, so need to get out of all loops
                clicking = false;
                more_boxes = false;
                close(f1);
                break;
            end
            % update the counter if we've made it this far
            counter = counter +1; 
            % if we didn't need to break, we're clear to add the
            % coordinates, and how we handle that will be conditional on if
            % this was the 1st or 2nd click for a particular box.
            if rem(counter,2) ~= 0
                % this is the 1st corner, so add it as such
                im_box = zeros(1,4);
                im_box(1) = round(x_i);
                im_box(2) = round(y_i);
            else
                % this is the 2nd corner, so we handle it
                im_box(3) = round(x_i);
                im_box(4) = round(y_i);
                % we'll show the user what that box looks like
                to_plot_x = [im_box(1),im_box(3),im_box(3),im_box(1),im_box(1)];
                to_plot_y = [im_box(2),im_box(2),im_box(4),im_box(4),im_box(2)];
                hold on
                h1 = plot(to_plot_x,to_plot_y,'r');
                drawnow
              
                % now we have to ask the user to give us the type or tell
                % us they would like to redraw
                prompt = "input category index for this box? Or in type 'r' to discard the box and reclick. ";
                response = input(prompt,"s");
                
                if strcmpi(response,'r')
                    % clear the coordinates
                    im_box(end,:) = [];
                    delete(h1);
                elseif str2num(response) > 0 && str2num(response) <= numel(patch_types)
                    % we'll be cute here and put the type index in the
                    % uppper corner of the box to label it
                    text(to_plot_x(1) + 20,to_plot_y(1) + 40,response,'Color','red','FontSize',14);
                  
                    % put the box values in a table for adding to output
                    [~,image_name,~] = fileparts(this_name);
                    image_name = {image_name}; % need to have a cell to put into table
                    patch_type = patch_types(str2num(response));
                    top_row = round(im_box(2) / p.Results.scaleRatio);
                    bot_row = round(im_box(4) / p.Results.scaleRatio);
                    left_col = round(im_box(1) / p.Results.scaleRatio);
                    right_col = round(im_box(3) / p.Results.scaleRatio);
                    table_to_add = table(image_name,patch_type,top_row,bot_row,left_col,right_col);

                    % place the new values in the overall table
                    patches_table = [patches_table; table_to_add];

                    % and write the table for safekeeping 
                    writetable(patches_table,p.Results.savePath);

                    % tell the user to keep going
                    disp('Awesome, good for you! Keep going!')
                end
            end
        end
    end

end