function subjectFolders = findSubjectFolders(rootPath, subjectPattern)
    % Recursively search for subdirectories matching the pattern

    subjectFolders = {}; % Initialize empty cell array
    allDirs = dir(rootPath);
    
    % Filter out non-directory entries and hidden folders
    allDirs = allDirs([allDirs.isdir] & ~startsWith({allDirs.name}, '.'));

    for j = 1:length(allDirs)
        folderName = allDirs(j).name;
        folderPath = fullfile(rootPath, folderName);
        
        % Only select folders that contain "sub-phy" in their name
        if contains(folderName, subjectPattern)
            subjectFolders{end+1} = folderName; % Add to list
        else
            %!! Die Idee war, dass noch weitere subfolder durchsucht
            %werden, allerdings findet es bei mir dann auch den "sub-phy"
            %in anderen (ungewollten) Pfaden) Müsste Optimiert werden

            %%Recursively search inside this folder
            %deeperSubjects = findSubjectFolders(folderPath, subjectPattern);
            %subjectFolders = [subjectFolders, deeperSubjects]; % Append results
        end
    end
end



