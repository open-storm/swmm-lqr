
%%
% Assumes storage curves have the same name as corresponding storage names
fid = fopen(inp);
A = fread(fid,'*char')';


% Create Maps to keep track of storage curves
storage_curves_h = containers.Map('keytype','char','valuetype','any');
storage_curves_A = containers.Map('keytype','char','valuetype','any');
storage_curves_As= containers.Map('keytype','char','valuetype','any');

% Find and navigate to the line where '[CURVES]' appears
i_storage = strfind(A,'[STORAGE]');
fseek(fid,i_storage,'bof');

% Read the first line
tline = fgets(fid);

% Keep reading each line until reading the next subsection (ie
% [TIMESERIES]), as marked by '['
while tline(1) ~= '['
    
    % Assume the below format, look for FUNCTIONAL Shapes (get the functional
    % curves below)
    % Name           Elev.    MaxDepth   InitDepth  Shape      Curve Name/Params (A,B,C)   Misc
    if strfind(tline,'FUNCTIONAL')
        %disp(tline)
        
        scan_results = textscan(tline,'%s %f %f %f %s %f %f %f');
        
        % Get the name of the storage curve
        s_node = scan_results{1}{:};
        
        % Parse storage curve parameters        
        h_max  = scan_results{3}(1);
        A      = scan_results{end-2}(1);
        B      = scan_results{end-1}(1);
        C      = scan_results{end  }(1);
        
        height = linspace(0, h_max, 11);
        area   = A*height.^B + C;        
        
        % Calculate average area for the storage curve
        As = trapz(height,area)/height(end);
        
        % Create a temporary Map for the heights, areas, and average_area
        tmp_h = containers.Map(s_node,height,'uniformvalues',0);
        tmp_A = containers.Map(s_node,area,'uniformvalues',0);          
        tmp_As= containers.Map(s_node,As,'uniformvalues',0);
        
        % Assume storage curve is complete once the latest line does not
        % contain the string 's_node'
        % --> Append to storage curve Maps
        storage_curves_h  = [storage_curves_h; tmp_h];
        storage_curves_A  = [storage_curves_A; tmp_A];
        storage_curves_As = [storage_curves_As;tmp_As];
                
    end
    
    tline = fgets(fid);
end

fclose(fid);
%%

fid = fopen(inp);
A = fread(fid,'*char')';

% Find and navigate to the line where '[CURVES]' appears
i_curves = strfind(A,'[CURVES]');
fseek(fid,i_curves,'bof');

% Read the first line
tline = fgets(fid);

% Keep reading each line until reading the next subsection (ie
% [TIMESERIES]), as marked by '['
while tline(1) ~= '['
    
    %disp(tline)
    
    % Assume the first line in each storage curve includes 'STORAGE'
    % and use that to mark when to save a specific storage curve
    if strfind(tline,'STORAGE')
        
        % Look for the format [NAME TYPE X-VALUE Y-VALUE] in the first line
        % of the storage curve
        scan_results = textscan(tline,'%s %s %f %f');
        
        % Parse the results
        s_node = scan_results{1}{:};
        height = scan_results{end-1};
        area   = scan_results{end};
        
        % Read each subsequent line 
        tline = fgets(fid);
        
        % .. and append results to the temporary arrays
        while strfind(tline,s_node)
            scan_results = textscan(tline,'%s %f %f');
            
            height = [height scan_results{end-1}];
            area   = [area scan_results{end}];
            
            tline = fgets(fid);
        end
        
        % Calculate average area for the storage curve
        As = trapz(height,area)/height(end);
        
        % Create a temporary Map for the heights, areas, and average_area
        tmp_h = containers.Map(s_node,height,'uniformvalues',0);
        tmp_A = containers.Map(s_node,area,'uniformvalues',0);          
        tmp_As= containers.Map(s_node,As,'uniformvalues',0);
        
        % Assume storage curve is complete once the latest line does not
        % contain the string 's_node'
        % --> Append to storage curve Maps
        storage_curves_h  = [storage_curves_h; tmp_h];
        storage_curves_A  = [storage_curves_A; tmp_A];
        storage_curves_As = [storage_curves_As;tmp_As];
        
    end
    
    tline = fgets(fid);
end

fclose(fid);

%%
% tline = fgets(fid);
% while ischar(tline)
%     disp(tline)
%     tline = fgets(fid);
% end
% 
% fclose(fid);