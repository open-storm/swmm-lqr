% 
% - outfall must connect to orifice
% - each storage node can have max 1 subcacthment

%%


A_mats = cell( length(nodes) );
G_mats = cell( length(nodes), length(lambdas) );
B_mats = cell( length(nodes), length(subcatch) );
C_mats_h = cell(length(lambdas), length(nodes) );
C_mats_Q = cell(length(lambdas), length(nodes) );
C_mats_OF= cell(length(outfalls),length(nodes) );

% Construct the block diagonals for the storage nodes and junctions
for m = 1:length(nodes)
    node = nodes{m};
    if ~ismember(node,{'','(null)'})
        if ismember(node,storages)
            if show_component
                disp(['storage' node])
            end
            
            %As = storage_curves_As(node);
            A_k = interp1(storage_curves_h(node), storage_curves_A(node), heights(node));
            A_k( heights(node) < 0.3048 ) = interp1(storage_curves_h(node), storage_curves_A(node), 0.3048); % if the area is 0, go up 1 foot
            
            % Read the storage curve parameters
            % swmm.get_from_input(inp,node,swmm.STORAGE_A) % Doesnt work
            A_mats{m,m} = create_storage_block_v2(T,A_k,lambdas(node));
            
        elseif ismember(node,junctions)
            if show_component
                disp(['junction' node])
            end
            
            i_tail   = find(ismember(nodes,node));
            %conduit = links{I(:,i_tail) == 1}
            conduit = links{I(i_tail,:) == 1};
            
            % the length of the conduit link (assumed in feet)
            l_conduit = swmm.get_from_input(inp,conduit,swmm.LENGTH);
            
            % Number of delay terms for flow 
            %  NOTE: 1000*(T/300) is a heuristic. Revisit using slope and length
            n_delay = ceil( l_conduit / ( 1000 *(T/300) ) );
            
            % Create the conduit block
            A_mats{m,m} = diag( ones(n_delay,1) ,-1);
            
        else
            % Then the node is either an outfall or divider
            % --> Assume it's an outfall
            % REVISIT and add a check to see if its in the "outfall" list
            if show_component
                disp([node])
            end
        end
    end
end

% Construct the off-diagonals for the orifices and conduits
for n = 1:length(links)
    link = links{n};
    if ~ismember(link,{'','(null)'})
        
        I_col   = find(ismember(links,link));
        
        % Get the names of the head and tailnodes for lookup
        head_node = nodes{I(:,I_col) == 1};
        tail_node = nodes{I(:,I_col) == -1};
        
        % Get the indices of the head and tail nodes
        i_head = find(ismember(nodes,head_node));
        i_tail = find(ismember(nodes,tail_node));
            
        % Create the off diagonal block
        A_mats{i_tail,i_head} = zeros( length(A_mats{i_tail,i_tail}), length(A_mats{i_head,i_head}) );        
        
        % Update the off diagonal block
        if ismember(link,orifices)

            % Assuming "head_node" is an orifice that drains a storage
            % node that is a key in the hash-map "lambdas"
            % .. and that "tail_node" is not an outfall node
            if ~isnan( lambdas(head_node) ) && ismember(tail_node,union(storages,junctions))
                A_mats{i_tail,i_head}(1,end) = lambdas(head_node);
            end
            
            if show_component
                disp(['orifice' link head_node tail_node])
            end
             
        elseif ismember(link,conduits)
            
            n_delay = size( A_mats{i_tail,i_head}, 2 ) - 1;
            
            if n_delay > 0
                A_mats{i_tail,i_head}(1, 1:end-1) = 1/n_delay;
            end
            
            if show_component
                disp(['conduit' link head_node tail_node])
            end
            
        else
            % Then the node is either an outfall or divider
            % --> Assume it's an outfall
            if show_component
                disp([link])
            end
        end
    end
end

% Fill in the rest of the zeros and create the A-matrix
% - Pre-allocate space for the A-matrix
% A = blkdiag(A_mats{find(eye(size(A_mats)))});

% Nevermind
A = [];

for m = 1:size(A_mats,1)
    for n = 1:size(A_mats,2)
        if isempty(A_mats{m,n})
            A_mats{m,n} = zeros( length(A_mats{m,m}), length(A_mats{n,n}) );
        end
    end

    if ~isempty(A_mats{m,m}) % 4/25: Check this hack. It's a workaround for orifice nodes bc A{m,m} is empty
        A = [A; A_mats{m,:}];
    end
end


%%


% B = zeros(size(A,1),2);
% B(1,1) = 1;
% B(3,2) = 1;

B   = [];
C_h = [];
C_Q = [];

% Determine the number of storage nodes, 
%  fill out the matrices with zeros,
n_storage = length(storages);

for m = 1:length(nodes) % or size(B_mats,1)
    for n = 1:length(subcatch)
        if ~isempty(A_mats{m,m})
            
            B_mats{m,n}        = A_mats{m,m}(:,1) * 0;
            
        end
    end
    
    for n = 1:n_storage % or size(B_mats,2)
        if ~isempty(A_mats{m,m})
            
            C_mats_h{n,m}      = A_mats{m,m}(1,:)' * 0;
            C_mats_Q{n,m}      = A_mats{m,m}(1,:)' * 0;
            
        end
    end
end

for n = 1:length(subcatch)
   sub = subcatch{n};
   if ~ismember(sub,{'','(null)'})
       
       % search through all nodes to find the index of the storage node
       % that corresponds with the subcatchment
       i_storage = find(ismember(nodes,subcatch_outlet(sub)));
       
       if i_storage > 0
           % Put in the 1 to assign a subcatchment to a storage node
           node = subcatch_outlet(sub);
           A_k = interp1(storage_curves_h(node), storage_curves_A(node), heights(node));
           if A_k < 0.3048 
               A_k = interp1(storage_curves_h(node), storage_curves_A(node), 0.3048); % if the area is 0, go up 1 foot           
           end
           
           B_mats{i_storage,n}(1) = T/A_k;
       end
       
   end
end

for m = 1:length(nodes)
    node = nodes{m};
    if ~ismember(node,{'','(null)'})
        
        % Find the storage node because that's where the runoff goes
        i_storage = find( ismember(storages,node) );
        
        if i_storage > 0
%             % Put in the 1's, lambda's where needed
%             B_mats{m,i_storage}(1) = 1;
                        
            C_mats_h{i_storage,m}(end) = 1;
            if ~(isnan(lambdas(node)))
                C_mats_Q{i_storage,m}(end) = lambdas(node);
            end
        end
%{
%         % in case of an outfall(s)
%         
%         i_outfall = find( ismember(outfalls,node) );
%         if i_outfall > 0
%             i_link = find( I(m,:) == -1 );
%             C_mats_OF{i_outfall,i_link} = A_mats{m,m}(1,:) * 0;
%             % Check if link is draining a controlled storage node
%             % or a junction.  Then insert lambda or 1 accordingly
%         end
%}
    end
    
    B = [B; B_mats{m,:}];
    C_h = [C_h; C_mats_h{:,m}];
    C_Q = [C_Q; C_mats_Q{:,m}];
    
end

C = [C_h'; C_Q'];

% C = zeros(3,size(A,2));
% C(1,2) = 1;
% C(2,4) = 1;
% C(3,6) = 1;


D = 0;

%%

for m = 1:length(nodes)
    for n = 1:length(orifices)
        if ~isempty( A_mats{m,m} )
            G_mats{m,n} = A_mats{m,m}(:,1)' * 0;
        end
    end 
end

% develop broad understanding with depth

G = [];

for n = 1:length(orifices)
    % Assign current link ID to "orifice"
    orifice = orifices{n};
    
    % Find which column the orifice belongs to in the incidence matrix, I,
    % so that we can find its corresponding head and tail node
    i_orifice = find( ismember(links,orifice) );
    
    % Find row corresponding to the head node
    i_head = find( I(:,i_orifice) == 1 );
    
    % Find row corresponding to the tail node
    i_tail = find( I(:,i_orifice) == -1 );
    
    head_node = nodes{i_head};
    tail_node = nodes{i_tail};
    
    %As = storage_curves_As(head_node);
    
    % Note: The lambda and As terms should correspond to the storage area of the
    % head node
    if ismember(head_node,junctions)                
        
        G_mats{i_head,n}(end) = I(i_head,i_orifice) * G_lambdas( head_node );
        
    elseif ismember(head_node,storages)
        
        A_k = interp1(storage_curves_h(head_node), storage_curves_A(head_node), heights(head_node));
        A_k( heights(head_node) < 0.3048 ) = interp1(storage_curves_h(head_node), storage_curves_A(head_node), 0.3048); % if the area is 0, go up 1 foot
        
        G_mats{i_head,n}(end) = I(i_head,i_orifice) * G_lambdas( head_node ) * T/A_k;
    else
        
        % the node is an outfall so G_mats{..}(end) returns an error
        
    end
    
    % Note: The lambda and As terms should correspond to the storage area of the
    % head node
    if ismember(tail_node,junctions)

        G_mats{i_tail,n}(end) = I(i_tail,i_orifice) * G_lambdas( head_node );
        
    elseif ismember(tail_node,storages)
        
        A_k = interp1(storage_curves_h(head_node), storage_curves_A(head_node), heights(head_node));
        A_k( heights(head_node) < 0.3048 ) = interp1(storage_curves_h(head_node), storage_curves_A(head_node), 0.3048); % if the area is 0, go up 1 foot
        
        G_mats{i_tail,n}(end) = I(i_tail,i_orifice) * G_lambdas( head_node ) * T/A_k;
    else
        
        % the node is an outfall so G_mats{..}(end) returns an error
        
    end
    
    G = [G; G_mats{:,n}];
end

G = G';
G(:,~isnan(lambda_array)) = [];
