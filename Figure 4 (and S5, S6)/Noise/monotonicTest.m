%Function to separate the noisy data for the single dose-responses by
%monotonicity.
%INPUTS:
    % p: 'struct' object that contains 
        % days: days of the experiment
        % cells: cell type of the experiment
    % cytokines: cytokine names
    % noise: previously generated noise data
%OUTPUTS:
    % total: vector containing the size of each "compartment"
        % first element is the size of the monotonic increasing noisy data
        % second element is the size of the monotonic decreasing noisy data
        % third element is the size of the multiphasic noisy data
    % monotonicInc: matrix containing the noisy data that is monotonic 
    % increasing
    %monotonicDec: matrix containing the noisy data that is monotonic
    %decreasing
    %multiphasic: matrix containing the noisy data that is multiphasic

function [total,monotonicInc,monotonicDec,multiphasic] = monotonicTest(p,cytokines,noise)
        
            y_bar_val = [];
            x_bar_val = cytokines;

            for k = 1:length(cytokines)
                 p.name = num2str(cell2mat(cytokines(k)));
    
                 Y = noise.(p.days).(p.cell).(p.name);
    
                 monotonicInc.(p.name) = [];
                 monotonicDec.(p.name) = [];
                 multiphasic.(p.name) = [];
        
                for i = 1:size(Y,1)
                     idx = find(Y(i,:)==max(Y(i,:)));
                     if (idx== length(Y(i,:)))
                         monotonicInc.(p.name) = [monotonicInc.(p.name);Y(i,:)];
                     elseif (idx== 1)
                         monotonicDec.(p.name) = [monotonicDec.(p.name);Y(i,:)];
                     else
                         multiphasic.(p.name)= [multiphasic.(p.name); Y(i,:)];
                     end 
                end 

                total.(p.name) = [size(monotonicInc.(p.name),1),size(monotonicDec.(p.name),1),size(multiphasic.(p.name),1)];

                y_bar_val = [y_bar_val;total.(p.name)];
            end 

    fig = figure('Visible','off');
    b = bar(x_bar_val,y_bar_val,"stacked",'FaceColor','flat');
    b(1).CData = [0.54 0.75 0.58]; 
    b(2).CData = [0.75 0.58 0.54]; 
    b(3).CData = [0.58 0.54 0.75];
 
    fontsize(20,"points")
    yticks([0 500 1000])

    file_name = strcat('Results/Figures_',p.days,'_',p.cell,'_Bar','.fig');
    savefig(fig, file_name); 
end 