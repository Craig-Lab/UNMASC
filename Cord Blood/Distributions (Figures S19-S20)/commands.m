clear variables

%Load MCMC parameter distributions from single fits
load("MCMC_samples.mat")

%Load single fit parameters
file_name1 = "single_fits_CD7_cells_MCMC.xlsx";

singleFits.CD7pCD5n = readtable(file_name1,"Sheet","CD7pCD5n",'ReadRowNames',true);
singleFits.CD7pCD5p = readtable(file_name1,"Sheet","CD7pCD5p",'ReadRowNames',true);


%Load combination fit parameters
filename = "combination_fits_cordBloodData.xlsx";

combinationFits.CD7pCD5n = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "CD7pCD5n");

combinationFits.CD7pCD5p = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "CD7pCD5p");


% Generate all figures for new CD7 data

cellTypes = ["CD7pCD5n", "CD7pCD5p"];
vials = ["vial2", "vial3", "vial4"];
interactions = ["SCF_IL7", "IL3_IL7"];

cytokinesByInteraction.SCF_IL7 = ["SCF", "IL7"];
cytokinesByInteraction.IL3_IL7 = ["IL3", "IL7"];

for c = 1:numel(cellTypes)

    p.cell = cellTypes(c);

    for v = 1:numel(vials)

        p.vial = vials(v);

        for i = 1:numel(interactions)

            p.interaction = interactions(i);

            cytokines = cytokinesByInteraction.(p.interaction);

            p.drug1 = cytokines(1);
            p.drug2 = cytokines(2);

            p.int = replace(p.interaction, "_", "-");

            p.combCol = matlab.lang.makeValidName( ...
                string(p.interaction) + "_" + string(p.vial));

            drug1_samples = samples.(p.vial).(p.interaction).(p.cell).(p.drug1);
            drug2_samples = samples.(p.vial).(p.interaction).(p.cell).(p.drug2);

            % Non-monotonic:
            % E0, E1, EC50_1, h1, EC50_2, h2, sigma  
            
            % Monotonic:
            % E0, EC50_1, h1, E1, sigma             

            if size(drug1_samples, 2) == 7
                p.drug1_mono = 0;
            else
                p.drug1_mono = 1;
            end

            if size(drug2_samples, 2) == 7
                p.drug2_mono = 0;
            else
                p.drug2_mono = 1;
            end

            fprintf("Generating figure for %s %s %s\n", ...
                string(p.cell), string(p.vial), string(p.int));

            fun_Dist_SingleFits(p, samples, singleFits, combinationFits);

        end
    end
end

