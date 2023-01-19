function sort_ims_imlogs(im_directory, im_ext, imlog_directory, sort_directory)
% This function will sort images into subfolders based on their wavelength
% which is stored in the image log spreadsheets.
%
% INPUT 
% im_directory              Path to the directory with the images
% im_ext                    extension for the images
% imlog_directory:          Path to the directory with the image logs
% sort_directory:           where should the sorting folders be 
%
% OUTPUT
% there are no outputs, but the files will be moved 
%
% Written by R.A. Manzuk (based on sortFiles.m by Bolton)
% 10/07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % if there is no start put it in front of the image extension so we can
    % find all the proper files
    if im_ext(1) ~= '*'
        im_ext = ['*' im_ext];
    end
    
    % get all of the image names in one cell array
    im_files = dir(fullfile(im_directory, im_ext));
    im_files(strncmp({im_files.name}, '.', 1)) = []; %remove files in dir starting with '.'
    
    % and get a set of extensionless names for later
    [~,im_names] = fileparts({im_files.name});
    
    % load in all of the image logs
    imlog_dir = dir(fullfile(imlog_directory,'*.csv'));
    imlog_dir(strncmp({imlog_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'
    imlogs_cell = arrayfun(@(s)readtable(fullfile(imlog_directory, s.name)),imlog_dir,'uniform',false);
    
    % and combine all the tables into one
    imlogs = imlogs_cell{1};
    for i = 2:numel(imlogs_cell)
        imlogs = [imlogs; imlogs_cell{i}];
    end
    
    % now we can creat a set of unique wavelengths from the imlogs
    wavelengths = unique(imlogs.im_wavelength);
    
    % Make the directories for each wavelength
    for i = 1:length(wavelengths)
      
        % This will check if the directory exists, if not, it will make them 
        mkdir([sort_directory '/' num2str(wavelengths(i)) 'nm']);
    end
    
    % Make a waitbar to keep track of how long it's taking
    f = waitbar(0, 'Sorting the images');
    
    for im = 1:size(imlogs,1)
        
        % Get the name of the image (extensionless)
        [~,this_name] = fileparts(imlogs.file_name{im});
        
        % Now find that image in the Output folder
        check = all(ismember(im_names,this_name),1);
        idx = find(check,1);
        
        % ok, time to move the file only do this part if we actually found a
        % match
        if sum(check) == 1
            % put together the path of this image
            file_to_move = fullfile(im_directory, [this_name im_ext(2:end)]);
            % and the path of its destination
            destination = fullfile(sort_directory, [num2str(imlogs.im_wavelength(im)) 'nm']);
            % and move it
            movefile(file_to_move,destination);
        end
       
        % update the waitbar
        waitbar(im/size(imlogs,1), f, 'Still Sorting')
    end
    
    % Close the waitbar
    close(f)
end
