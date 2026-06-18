clear variables

%Load MCMC parameter distributions from single fits
load("MCMC_samples.mat")

%Load single fit parameters
file_name1 = "single_fits_ipscs_MCMC.xlsx";

singleFits.D7.PROT = readtable(file_name1,"Sheet","D7 PROT",'ReadRowNames',true);
singleFits.D7.CD4ISP = readtable(file_name1,"Sheet","D7 CD4ISP",'ReadRowNames',true);
singleFits.D7.CD3N = readtable(file_name1,"Sheet", "D7 CD3N",'ReadRowNames',true);
singleFits.D14.CD3N = readtable(file_name1,"Sheet", "D14 CD3N",'ReadRowNames',true);
singleFits.D14.CD3P = readtable(file_name1,"Sheet", "D14 CD3P",'ReadRowNames',true);
singleFits.D14.SP = readtable(file_name1,"Sheet", "D14 SP",'ReadRowNames',true);


%Load combination fit parameters
filename = "combination_fits_ipscs.xlsx";

combinationFits.D7.PROT = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "D7 PROT");

combinationFits.D7.CD4ISP = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "D7 CD4ISP");

combinationFits.D7.CD3N = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "D7 CD3N");

combinationFits.D14.CD3N = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "D14 CD3N");

combinationFits.D14.CD3P = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "D14 CD3P");

combinationFits.D14.SP = readtable(filename, "ReadRowNames", true, ...
    "VariableNamingRule", "preserve", "Sheet", "D14 SP");


%% Generate all figures

sheetInfo = {
    "D7",  "PROT";
    "D7",  "CD4ISP";
    "D7",  "CD3N";
    "D14", "CD3N";
    "D14", "CD3P";
    "D14", "SP"
};

for s = 1:size(sheetInfo, 1)

    p.days = sheetInfo{s, 1};
    p.cell = sheetInfo{s, 2};

    comboNames = string(combinationFits.(p.days).(p.cell).Properties.VariableNames);

    for c = 1:length(comboNames)

        p.int = comboNames(c);

        drugs = split(p.int, "-");
        p.drug1 = drugs(1);
        p.drug2 = drugs(2);

        drug1_samples = samples.(p.days).(p.cell).(p.drug1);
        drug2_samples = samples.(p.days).(p.cell).(p.drug2);

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
            string(p.days), string(p.cell), string(p.int));

        fun_Dist_SingleFits(p, samples, singleFits, combinationFits);

        close(gcf)

    end
end


