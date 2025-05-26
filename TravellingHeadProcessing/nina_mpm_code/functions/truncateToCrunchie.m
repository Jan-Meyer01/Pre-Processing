% This function helps to break down each (dynamic) path to crunchie and to later search for the raw data

function truncatedPath = truncateToCrunchie(fullPath)
    % Break the path into its parts
    pathParts = strsplit(fullPath, filesep);
    
    % Find the position of 'crunchie'
    idx = find(strcmp(pathParts, 'crunchie'), 1, 'last');
    
    % If 'crunchie' was found, rebuild the path to it
    if ~isempty(idx)
        truncatedPath = fullfile(pathParts{1:idx});
    else
        warning('Ordner "crunchie" nicht im Pfad gefunden.');
        truncatedPath = fullPath; % Falls nicht gefunden, bleibt der ursprüngliche Pfad erhalten
    end
end