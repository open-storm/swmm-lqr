% read_subcatch_outlet

%%
% Assumes storage curves have the same name as corresponding storage names
fid = fopen(inp);
A = fread(fid,'*char')';


% Create Maps to keep track of storage curves
subcatch_outlet = containers.Map('keytype','char','valuetype','any');
subcatch_area   = containers.Map('keytype','char','valuetype','any');

% Find and navigate to the line where '[CURVES]' appears
i_storage = strfind(A,'[SUBCATCHMENTS]');
fseek(fid,i_storage,'bof');

% Read the first line 'SUBCATCHMENTS]' 
tline = fgets(fid);
% and skip to the second line ';;Name       Rain Gage ...'
tline = fgets(fid);

% Keep reading each line until reading the next subsection (ie
% [SUBAREAS]), as marked by '['
while tline(1) ~= '['
    % Assume the format below
    % Name           Rain Gage        Outlet           Area     %Imperv  Width    %Slope   CurbLen  SnowPack 
    
    % Avoid parsing commented lines
    if tline(1) ~= ';'
        scan_results = textscan(tline,'%s %s %s %f %f %f %f %f');
        
        % Avoid parsing a blank line
        if ~isempty(scan_results{1})
            
            tmp = containers.Map(scan_results{1}{:},scan_results{3}{:},'uniformvalues',0);
            subcatch_outlet = [subcatch_outlet; tmp];
            
            tmp_area = containers.Map(scan_results{1}{:},scan_results{4},'uniformvalues',0);
            subcatch_area   = [subcatch_area; tmp_area];
            
        end
    end
    
    tline = fgets(fid);
end

fclose(fid);