% read_orifice_from_node

%%
% Assumes storage curves have the same name as corresponding storage names
fid = fopen(inp);
A = fread(fid,'*char')';


% Create Maps to keep track of storage curves
orifice_from_node = containers.Map('keytype','char','valuetype','any');

% Find and navigate to the line where '[CURVES]' appears
i_storage = strfind(A,'[ORIFICES]');
fseek(fid,i_storage,'bof');

% Read the first line 'ORIFICES]' 
tline = fgets(fid);
% and skip to the second line ';;Name       From Node ...'
tline = fgets(fid);

% Keep reading each line until reading the next subsection (ie
% [XSECTIONS]), as marked by '['
while tline(1) ~= '['
    % Assume the format below
    % Name           From Node        To Node          Type         Offset     Qcoeff     Gated    CloseTime 
    
    % Avoid parsing commented lines
    if tline(1) ~= ';'
        scan_results = textscan(tline,'%s %s %s %s %f %f %s %f');
        
        % Avoid parsing a blank line
        if ~isempty(scan_results{1})
            
            tmp = containers.Map(scan_results{2}{:},scan_results{1}{:},'uniformvalues',0);
            
            orifice_from_node = [orifice_from_node; tmp];
            
        end
    end
    
    tline = fgets(fid);
end

fclose(fid);

node_from_orifice = containers.Map(orifice_from_node.values,orifice_from_node.keys);