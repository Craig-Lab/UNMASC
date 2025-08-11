%Function to generate data points and noisy data based on the polynomial
%coefficients in Michaels et al., 2022. 
%INPUTS:
    % k: number associated with the cytokine. Each cytokine was attributed a number from
        % 1 to 6.
    % p: 'struct' object containing the name of the cytokine, the days of
        % the experiment and the cell type.
    % y: the function handles
    % scaled_conc: scaled concentration of the cytokines
    % cytConc: concentrations for each cytokine (in ng/mL)
    % coeff: coefficients of the polynomials
    % stde: standard error of the polynomial coefficients

%OUTPUTS:
    %noise_data: generated noise data (each row is a new set of noisy data
    %and each column corresponds to a different cytokine concentration)

function [noise_data,data] = dataNoise_fun(k,p,y,scaled_conc,cytConc,coeff,stde)
    
%Initialize variables
    fun_handle = y.(p.days).(p.cell).(p.name);
    conc = scaled_conc;
    conc_ng = cytConc(:,k);
    poly_coeff = coeff.(p.days).(p.cell).(p.name);
    stde_val = stde.(p.days).(p.cell).(p.name);

%Generate data points
data = (fun_handle(conc)).^2./max((fun_handle(conc)).^2);

%Generate noise data
noise_data = noise_fun(fun_handle,conc,poly_coeff,stde_val);


% err_neg = data - min(noise_data);
% err_pos = max(noise_data) - data;

% fig = figure('Visible','off');
% plot(log10(conc_ng),noise_data,'Color','#d5d5d5','LineWidth',2)
% hold on
% scatter(log10(conc_ng),data,'Color',"#1468a2",'LineWidth',2)
% hold on
% errorbar(log10(conc_ng),data,err_neg,err_pos,'o','MarkerSize',5,'MarkerFaceColor',...
%     "#A2142F",'MarkerEdgeColor',"#A2142F",'LineWidth',2,'Color',"#1468a2")
% % title(name)
% fontsize(20,"points")
% % ylim([0 1.2])
% % yticks([])
% % xticks([])
% axis square
% 
%  file_name = strcat('Noise_figures/',p.days,'/',p.cell,'/',p.name,'.fig');
%  savefig(fig, file_name); % Saves the figure 