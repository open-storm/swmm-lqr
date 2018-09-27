close all

max_col = 4;
plot_linear_model = 0;%1;

i_subcatch = [1:length(storages)];
n_selected = length(i_subcatch);
n_figures  = ceil( n_selected/max_col );

% Remove "yellow" colors that are too light
ccontrol = linspecer(n_selected*2);
if n_selected > 20 
    ccontrol(11:20,:)  = [];
elseif n_selected > 11
    ccontrol(11:end,:) = [];
end

% Set default last color to black
ccontrol(end+1,:)   = [0 0 0];

subcatch_outlet_keys = subcatch_outlet.keys;

% Extract the model results 
% - Assume model has been run and n_storage exists
% - Assume y_height and y_flow have the same indices for given storage node
%      - Make sure of this when running the simulation
y_height = y(:,1:n_storage);
y_flow   = y(:,1+n_storage:end);

if plot_linear_model
    y_height_linear = y_linear(:,1:n_storage);
    y_flow_linear   = y_linear(:,1+n_storage:end);
end

if size(y,1) == length(t)
    time = t;
elseif size(y,1) == length(t_swmm)
    time = t_swmm;
else
    warning('check time vector')
end

for m = 1:n_figures
    f = figure('name',num2str(m),'units','normalized','outerposition',[0 0 1 1]);%figure(m);
    
    % Make sure not to exceed max number of storage elements
    if m*max_col > length( i_subcatch )
        ind = i_subcatch( 1+(m-1)*max_col:length(i_subcatch) );
    else
        ind = i_subcatch( 1+(m-1)*max_col:m*max_col );
    end
    
    for n = 1:length(ind)
        
        %disp( (m-1)*max_col + n )
        n_color = (m-1)*max_col + n;
        
        %{
        % Notes:
          - Place storage name as title
          - Add legends
          - Skip light colors from linspecer
          - Eventually, plot both blinearized and matswmm results
        %}
        
        % Precip, runoff plots
        n_subplot     = (n-1)+1;
        ax(n_subplot) = subplot(4,length(ind), n_subplot );
        outlets = subcatch_outlet_keys(ismember(subcatch_outlet.values,storages{ind(n)}));
        plot(  t, runoff(:,ismember(subcatch,outlets)) ,'-','color',ccontrol(n_color,:),'linewidth',3)
        %plot(  t, runoff ,'-','color',ccontrol(n_color,:),'linewidth',3)
        title(storages{n_color},'fontsize',20)
        if n == 1
            ylabel('Rainfall [in]','fontsize',15)
        end
        
        % Storage height
        n_subplot     = (n-1)+length(ind)*1+1;
        ax(n_subplot) = subplot(4,length(ind), n_subplot );
        plot( t, depth(:,ind(n)),'--','color',ccontrol(end,:),'linewidth',2)
        hold on
        plot( time, y_height(:,ind(n)),'-','color',ccontrol(n_color,:),'linewidth',3)
        if plot_linear_model
            plot( t, y_height_linear(:,ind(n)),'-','color','red')
        end
        hold off
        if n ==  1
            ylabel('Pond Height [ft]','fontsize',15)
        end
        ylim([0 15])
        
        % Storage outflow
        n_subplot     = (n-1)+length(ind)*2+1;
        ax(n_subplot) = subplot(4,length(ind), n_subplot );
        % /!\ IMP: Remember to find the orifices corresponding to the
        % storage node
        i_orifice = find( ismember(orifices, orifice_from_node(storages{ind(n)}) ) );
        plot( t, flow(:,i_orifice),'--','color',ccontrol(end,:),'linewidth',2)
        hold on
        plot( time, y_flow(:,ind(n)),'-','color',ccontrol(n_color,:),'linewidth',3)
        if plot_linear_model
            plot( t, y_flow_linear(:,ind(n)),'-','color','red')
        end
        hold off
        if n == 1
            ylabel('Discharge [ft^{3} s^{-1}]','fontsize',15)
        end
        
        % Orifice percent opening
        n_subplot     = (n-1)+length(ind)*3+1;
        ax(n_subplot) = subplot(4,length(ind), n_subplot );
        %{
        - Find the orifices corresponding to the storage node
        - Use those indices to find their gate positions
        - Plot        
        %}
        plot( time, OR_setting(:,ind(n)),'-','color',ccontrol(n_color,:),'linewidth',3)
        hold on
        if plot_linear_model
            plot( t, OR_setting_linear(:,ind(n)),'-','color','red')
        end
        hold off
        if n == 1
            ylabel('Gate Position [%]','fontsize',15)
        end
        xlabel('Timestep','fontsize',15)
        ylim([0 1])
        
    end
    
    %hgexport(f,[figure_path num2str(m)])
    if save_plot
        export_fig([figure_path num2str(m) '.png'],'-nocrop','-transparent','-r150',f)
    end

end