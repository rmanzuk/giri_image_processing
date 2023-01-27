function rename_dir(input_folder, ext, name_map)
% Working with the capture one given names of the images as opposed to
% sample names got annoying. This goes in and renames the files given a
% naming map.
%
% INPUT 
% input_folder: String with the full path to files to be renamed
% ext: String with the extension for the input files
% name_map: (optional) n_sample x 2 cell array with the file names of the
%
% OUTPUT
% there are no outputs, but the files will be renamed
%
% Written by R.A. Manzuk
% Thursday, January 26, 2023 at 5:42:05 PM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % if there is no star put it in front of the image extension so we can
    % find all the proper files
    if ext(1) ~= '*'
        ext = ['*' ext];
    end

    % make a directory of the input to start
    input_dir = dir(fullfile(input_folder, ext));
    input_dir(strncmp({input_dir.name}, '.', 1)) = []; %remove files in dir starting with '.'

    % loop through all files
    for i = 1:numel(input_dir)
        % get just the name
        [~,in_name,~] = fileparts(input_dir(i).name); 
        % figure out what we're searching to match
        [~,search_names,~] = cellfun(@fileparts,name_map(:,1),'UniformOutput',false);
        % get the index of the sample name
        map_ind = ismember(search_names,in_name);
    
        % now we can make the output name
        output_name = [name_map{map_ind,2} ext(2:end)];

        % use movefile to rename
        movefile([input_folder '/' in_name ext(2:end)], [input_folder '/' output_name]);
    end
end
