%% analyze microbiome data

load('..\..\data\Microbiome_Interactios_Table.mat')

species = unique(Microbiome_Interactios_Table.Species);
species(839:844) = []; % they are not microbial species
metabolites = unique(Microbiome_Interactios_Table.Smallmoleculemetaboliteormacromolecule);
n = length(species);
m = length(metabolites);

export = zeros(n,m);
import = export;
degregation = import;

for i=1:size(Microbiome_Interactios_Table,1)
   I = strfind(species,Microbiome_Interactios_Table.Species(i));
   J = find(ismember(metabolites,Microbiome_Interactios_Table.Smallmoleculemetaboliteormacromolecule(i)));
   if contains(cellstr(Microbiome_Interactios_Table.Metabolicactivity(i)),'export') 
       export(I,J) = 1;
   elseif contains(cellstr(Microbiome_Interactios_Table.Metabolicactivity(i)),'import') 
       import(I,J) = 1;
   else
       degregation(I,J) = 1;
   end
end

norm_import = import./sum(import,2);
norm_import(isnan(norm_import)) = 0;

norm_export = export./sum(export,1);
norm_export(isnan(norm_export)) = 0;

competition = norm_import*norm_import'; % !!!
complementarity = norm_import*norm_export';

competition(1:n+1:end) = 0;
complementarity(1:n+1:end) = 0;


%%

save('..\..\data\MicrobiomeNetworks','species','metabolites','competition','complementarity','norm_import')
