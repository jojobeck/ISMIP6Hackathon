%% Run ANOVA ON ISMIP6 AIS2300 DATA
% The example here is a 3-way anova test

%% Define data to be explained
% This script is based on Helene Seroussi's anova test for the AIS 2300 ISMIP6 paper currently in review

varn = 'shelfmelt'; 

%% Path to data, use mounted drive here 

% local
%computed_name = '/home/ronja/projects/MeltSensitivity/10528582/ComputedScalars/ComputedScalars'
% mounted nci
%computed_name = '/mnt/ronja/nci/ismip6_hackathon/ComputedScalars'

computed_name = '/home/ronja/projects/MeltSensitivity/ComputedScalarsHelene'


%% Load metadata information

% FIXME update with Jojos new structure
%metadata = readtable("../../Metadata.txt");
%metadata = readtable("Metadata_icefront.txt"); %FIXME is this the latest??
%metadata = readtable("../../data_preprocessing/Metadata_IceFront_InitMethod.txt"); %FIXME is this the latest??
metadata = readtable("../../data_preprocessing/Metadata_IceFront_InitMethod_GIA_StressBalance.txt"); %FIXME is this the latest??


% %% Exapmle use
% % List all variables in the table
% metadata.Properties.VariableNames
% % access a column called 'ClimateForcing' by typing
% metadata(:,'ClimateForcing')
% metadata.ClimateForcing
% % Select a specific entry 
% i = (strcmp(metadata.fileID, 'DC_ISSM') & strcmp(metadata.Experiment,'expAE02')) ; 
% metadata(i,'MeltParameterisation')
% metadata.MeltParameterisation(i)

%%

% Determine Melt Sensitivity groups
ms = readtable("tables/best_ms_method_per_model.csv");


%%
yearlen = 360*24*60*60; %FIXME 360 day year OK?


%% Load data 

sector = "AIS" ;% "FRIS", %"RIS", "ASE", "Aurora", "Wilkes"

data_alltimeseries = NaN*ones(4*8,285);

expnum=0;
for iexp=1:length(metadata.Model)
    % skip if not correct experiment
    if strcmp(convertCharsToStrings(metadata.Experiment{iexp}),"ctrlAE")
        continue
    %elseif strcmp(convertCharsToStrings(metadata.MainSubmission{iexp}),"False")
    %    continue    % FIXME continue to exclude these?
    %elseif (strcmp(convertCharsToStrings(metadata.Group{iexp}),"IMAU") | strcmp(convertCharsToStrings(metadata.Group{iexp}),"DOE") )
    %    continue   % FIXME continue to exclude these?
    elseif (strcmp(convertCharsToStrings(metadata.Group{iexp}),"IMAU"))
        filename = append(computed_name,'/',metadata.Experiment{iexp},'/',varn,'/', ...
            'jb_computed_', varn, '_AIS_', metadata.fileID{iexp},'_',metadata.Experiment{iexp},'.nc');
    else
        filename = append(computed_name,'/',metadata.Experiment{iexp},'/',varn,'/', ...
            'computed_', varn, '_AIS_', metadata.fileID{iexp},'_',metadata.Experiment{iexp},'.nc');
    end
        if strcmp(sector, "ASE")
            data=ncread(filename,append(varn+"_sector_4")) ; 
        elseif strcmp(sector,"Aurora")
            data=ncread(filename,append(varn+"_sector_8")) ; 
        elseif strcmp(sector, "Wilkes")
            data=ncread(filename,append(varn+"_sector_7")) ; 
        elseif strcmp(sector, "FRIS")
            data1=ncread(filename,append(varn+"_sector_5")) ;
            data2=ncread(filename,append(varn+"_sector_11")) ; 
            data = data1+data2;
        elseif strcmp(sector, "RIS")
            data1=ncread(filename,append(varn+"_sector_2")) ;
            data2=ncread(filename,append(varn+"_sector_6")) ; 
            data = data1+data2;
        else
            %fprintf("No sector selected (correctly), loading whole AIS ")
            data=ncread(filename,append(varn)) ; 
       end;
        
        %time=2016:2016+length(data)-1;

        expnum = expnum+1;

		% FOR dynamic SLE
        %data_alltimeseries(expnum,:)= data(1:285)*1e12/1028/3.6e14;  % Gt convert to m SLE%-(data(1:285)-data(1))/(362.5*1000); %FIXME replace anomaly with ctrl?
		% FOR BMB
        %data_alltimeseries(expnum,:)=-(data(1:285)-data(1))/1e12*yearlen;
        data_alltimeseries(expnum,:)=cumsum(-(data(1:285)-data(1))/1e12*yearlen); %FIXME replace anomaly with ctrl?
end

% Save data 
%%
save("data/data_alltimeseries_"+sector+"_cumBMR.mat", 'data_alltimeseries');
%save("data/time_"+sector+".mat")

%% Load data if you want to skip the part before
sector = "AIS" ;% "FRIS", %"RIS", "ASE", "Aurora", "Wilkes"

data_alltimeseries = load("data/data_alltimeseries_"+sector+"_cumBMR.mat", 'data_alltimeseries').data_alltimeseries;

%%
% these are string that sort the data into classes, e.g., subgrid melt
% "yes" or "no", these groups will be tested in the ANOVA test
ice_alltimeseries = {};
climate_alltimeseries = {};
calving_alltimeseries = {};
meltsens_alltimeseries = {};
meltsens_alltimeseries_values = {};
meltparameterisation_alltimeseries = {};
meltparameters_alltimeseries = {};
glresolution_alltimeseries = {};
subglmelt_alltimeseries = {};
resolution_alltimeseries = {};
init_alltimeseries = {};
gia_alltimeseries = {};
stressbalance_alltimeseries = {};

% Calving groups:
%calvinggroup1 = {'VUW_PISM1', 'VUW_PISM1_s1', 'VUW_PISM1_s2', 'VUW_PISM1_s3', 'VUW_PISM1_s4' ,'VUW_PISM2' ,'VUW_PISM2_s1' ,'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4', 'PIK_PISM', 'LSCE_GRISLI2' ,'LSCE_GRISLI','UCM_Yelmo'}
%calvinggroup2 = {'DC_ISSM', 'ILTS_SICOPOLIS' ,'VUB_AISMPALEO','NCAR_CISM1' ,'NORCE_CISM3-MAR364-ERA-t1-nonlocal', 'NORCE_CISM4-MAR364-ERA-t1-nonlocal' ,'NORCE_CISM5-MAR364-ERA-t1-nonlocal','NCAR_CISM2', 'NORCE_CISM3-MAR364-ERA-t1-local' ,'NORCE_CISM4-MAR364-ERA-t1-local' ,'NORCE_CISM5-MAR364-ERA-t1-local', 'UNN_Ua','NORCE_CISM2-MAR364-ERA-t1', 'NORCE_CISM3-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-JRA-t1','NORCE_CISM5-MAR364-ERA-t1', 'ULB_fETISh-KoriBU2', 'IGE_ElmerIce' }
%%calvinggroup3 = {'ULB_fETISh-KoriBU1', 'UCSD_ISSM', 'UTAS_ElmerIce'}
%calvinggroup3 = {'ULB_fETISh-KoriBU1', 'UCSD_ISSM', 'DOE_MALI_4km', 'DOE_MALI_8km_Ant95', 'DOE_MALI_8km_AntMean','IMAU_UFEMISM1', 'IMAU_UFEMISM2','IMAU_UFEMISM3','IMAU_UFEMISM4',  'UTAS_ElmerIce'}
calvinggroup1 = {'VUW_PISM1', 'VUW_PISM1_s1', 'VUW_PISM1_s2', 'VUW_PISM1_s3','VUW_PISM1_s4' ,'VUW_PISM2','VUW_PISM2_s1', 'VUW_PISM2_s2', 'VUW_PISM2_s3', 'VUW_PISM2_s4', 'PIK_PISM', 'LSCE_GRISLI2', 'LSCE_GRISLI','UCM_Yelmo','IMAU_UFEMISM1', 'IMAU_UFEMISM2', 'IMAU_UFEMISM3' ,'IMAU_UFEMISM4'}
calvinggroup2 = {'DC_ISSM', 'ILTS_SICOPOLIS' ,'VUB_AISMPALEO','NCAR_CISM1',    'NORCE_CISM3-MAR364-ERA-t1-nonlocal','NORCE_CISM4-MAR364-ERA-t1-nonlocal','NORCE_CISM5-MAR364-ERA-t1-nonlocal','NCAR_CISM2', 'NORCE_CISM3-MAR364-ERA-t1-local','NORCE_CISM4-MAR364-ERA-t1-local' , 'NORCE_CISM5-MAR364-ERA-t1-local','UNN_Ua','NORCE_CISM2-MAR364-ERA-t1', 'NORCE_CISM3-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-ERA-t1', 'NORCE_CISM4-MAR364-JRA-t1','NORCE_CISM5-MAR364-ERA-t1', 'UTAS_ElmerIce','ULB_fETISh-KoriBU2','IGE_ElmerIce','DOE_MALI_4km','DOE_MALI_8km_Ant95' ,'DOE_MALI_8km_AntMean'}
calvinggroup3 = {'ULB_fETISh-KoriBU1', 'UCSD_ISSM'}


yearlen = 360*24*60*60; %FIXME 360 day year OK?
expnum=0;

for iexp=1:length(metadata.Model)
    if strcmp(convertCharsToStrings(metadata.Experiment{iexp}),"ctrlAE")
        continue
    %elseif strcmp(convertCharsToStrings(metadata.MainSubmission{iexp}),"False")
    %    continue % FIXME keep consistent with above loop
    %elseif (strcmp(convertCharsToStrings(metadata.Group{iexp}),"IMAU") | strcmp(convertCharsToStrings(metadata.Group{iexp}),"DOE") )
    %    continue % FIXME keep consistent with above loop
    else
        
        expnum = expnum+1;

        %  ice flow model
        ice_alltimeseries{end+1}= metadata.Model{iexp};
		
        % climate model
        if strcmp(convertCharsToStrings(metadata.Experiment{iexp}),"expAE02")
			climate_alltimeseries{end+1}='CCSM4'; %FIXME CHECK
		elseif strcmp(convertCharsToStrings(metadata.Experiment{iexp}),"expAE03")
			climate_alltimeseries{end+1}='HadGEM2';
	    elseif strcmp(convertCharsToStrings(metadata.Experiment{iexp}),"expAE04")
			climate_alltimeseries{end+1}='CESM2';
		elseif strcmp(convertCharsToStrings(metadata.Experiment{iexp}),"expAE05")
			climate_alltimeseries{end+1}='UKESM';
		else 
			error('exp not supported yet');
        end
                
        % Add melt sensitivity group based on melt sens       
        i = (strcmp(ms.Model, metadata.fileID{iexp}) ) ; 
        
        if strcmp(sector, "ASE")
            sens = ms.melt_sensetivity_Amundsen(i); 
            small_sens = mean(ms.melt_sensetivity_Amundsen) - std(ms.melt_sensetivity_Amundsen)/2; % smaller than mean - half a standard deviation
            mean_sens = mean(ms.melt_sensetivity_Amundsen) ; % 
            large_sens = mean(ms.melt_sensetivity_Amundsen) + std(ms.melt_sensetivity_Amundsen)/2;
        elseif strcmp(sector,"Aurora") % FIXME WHOLE EAIS
            sens = ms.melt_sensetivity_Aurora(i);
            small_sens = mean(ms.melt_sensetivity_Aurora) - std(ms.melt_sensetivity_Aurora)/2; % smaller than mean - half a standard deviation
            mean_sens = mean(ms.melt_sensetivity_Aurora) ; % 
            large_sens = mean(ms.melt_sensetivity_Aurora) + std(ms.melt_sensetivity_Aurora)/2;
        elseif strcmp(sector, "Wilkes") % FIXME WHOLE EAIS
            sens = ms.melt_sensetivity_Wilkes(i);
            small_sens = mean(ms.melt_sensetivity_Wilkes) - std(ms.melt_sensetivity_Wilkes)/2; % smaller than mean - half a standard deviation
            mean_sens = mean(ms.melt_sensetivity_Wilkes) ; % 
            large_sens = mean(ms.melt_sensetivity_Wilkes) + std(ms.melt_sensetivity_Wilkes)/2;
        elseif strcmp(sector, "FRIS" )
            sens = ms.melt_sensetivity_FilchnerRonne(i);
            small_sens = mean(ms.melt_sensetivity_FilchnerRonne) - std(ms.melt_sensetivity_FilchnerRonne)/2; % smaller than mean - half a standard deviation
            mean_sens = mean(ms.melt_sensetivity_FilchnerRonne) ; % 
            large_sens = mean(ms.melt_sensetivity_FilchnerRonne) + std(ms.melt_sensetivity_FilchnerRonne)/2;
        elseif strcmp(sector, "RIS")
            sens = ms.melt_sensetivity_Ross(i);
            small_sens = mean(ms.melt_sensetivity_Ross) - std(ms.melt_sensetivity_Ross)/2; % smaller than mean - half a standard deviation
            mean_sens = mean(ms.melt_sensetivity_Ross) ; % 
            large_sens = mean(ms.melt_sensetivity_Ross) + std(ms.melt_sensetivity_Ross)/2;
        else
            %fprintf("No sector selected (correctly), loading whole AIS ")
            sens = ms.melt_sensetivity_AIS(i); 
            very_small_sens = mean(ms.melt_sensetivity_AIS) - std(ms.melt_sensetivity_AIS); % smaller than mean - a standard deviation
            small_sens = mean(ms.melt_sensetivity_AIS) - std(ms.melt_sensetivity_AIS)/2; % smaller than mean - half a standard deviation
            mean_sens = mean(ms.melt_sensetivity_AIS) ; % 
            large_sens = mean(ms.melt_sensetivity_AIS) + std(ms.melt_sensetivity_AIS)/2;
            very_large_sens = mean(ms.melt_sensetivity_AIS) + std(ms.melt_sensetivity_AIS);
        end

        if sens <= very_small_sens
             sens_group='very low melt sensitivity';
        elseif sens <= small_sens
             sens_group='low melt sensitivity';
        elseif sens <= large_sens
            sens_group='medium melt sensivity';
        elseif sens <= very_large_sens
            sens_group='high melt sensivity';
        else 
             sens_group='very high melt sensitivity';
        end
%         if sens <= mean_sens
%             sens_group='low melt sensitivity'
%         else 
%             sens_group='high melt sensitivity';
%         end            
        meltsens_alltimeseries_values{end+1} = sens; 
        meltsens_alltimeseries{end+1}=sens_group;
                
        % Define calving groups based on calving criertia
        i1 = ( strcmp(convertCharsToStrings(metadata.fileID{iexp}), calvinggroup1)) ;
        i2 = ( strcmp(convertCharsToStrings(metadata.fileID{iexp}), calvinggroup2)) ;
        i3 = ( strcmp(convertCharsToStrings(metadata.fileID{iexp}), calvinggroup3)) ;
        if sum(i1)>0
            calving_alltimeseries{end+1}='Strong calving';
        elseif sum(i2)>0
            calving_alltimeseries{end+1}='Weak calving';
        elseif sum(i3)>0
            calving_alltimeseries{end+1}='No calving';
            %calving_alltimeseries{end+1}='Weak calving';
        else
            error('no calving group assigned');
        end;

        meltparameterisation_alltimeseries{end+1} = metadata.MeltParameterisation{iexp};
        glresolution_alltimeseries{end+1} = metadata.ResolutionGL{iexp}; 
        subglmelt_alltimeseries{end+1} = metadata.GLMelt{iexp}; 
         
        %resolution_alltimeseries{end+1} = metadata.Resolution{iexp} ; 
        % only 2 groups with resolution
        if strcmp(convertCharsToStrings(metadata.Resolution{iexp}), '16km')
            resolution_alltimeseries{end+1} = 'coarse';
        elseif strcmp(convertCharsToStrings(metadata.initMethod{iexp}), '32km') 
            resolution_alltimeseries{end+1} = 'coarse';
        %elseif strcmp(convertCharsToStrings(metadata.initMethod{iexp}), '8km') 
        %    resolution_alltimeseries{end+1} = 'coarse';
        %elseif strcmp(convertCharsToStrings(metadata.initMethod{iexp}), '4km') 
        %    resolution_alltimeseries{end+1} = 'coarse';
        else
            resolution_alltimeseries{end+1} = 'fine';
        end;        

        %init_alltimeseries{end+1} = metadata.initMethod{iexp};
        % only two groups for init method
        if strcmp(convertCharsToStrings(metadata.initMethod{iexp}), 'SP')
            init_alltimeseries{end+1} = 'SP';
        elseif strcmp(convertCharsToStrings(metadata.initMethod{iexp}), 'SP+') 
            init_alltimeseries{end+1} = 'SP';
        else
            init_alltimeseries{end+1} = 'DA';
        end;

        %gia_alltimeseries{end+1} = metadata.GIA{iexp};
        % only two groups for init method
        if strcmp(convertCharsToStrings(metadata.GIA{iexp}), 'VE')
            gia_alltimeseries{end+1} = 'GIA';
        elseif strcmp(convertCharsToStrings(metadata.GIA{iexp}), 'ELRA') 
            gia_alltimeseries{end+1} = 'GIA';
        else
            gia_alltimeseries{end+1} = 'no GIA';
        end;

        %stressbalance_alltimeseries{end+1} = metadata.stressBalance{iexp};
        % only two groups for init method
        if strcmp(convertCharsToStrings(metadata.stressBalance{iexp}), 'SIA+SSA')
            stressbalance_alltimeseries{end+1} = 'SSA';
        %elseif strcmp(convertCharsToStrings(metadata.stressBalance{iexp}), 'L1L2') 
        %    gia_alltimeseries{end+1} = 'GIA';
        else
            stressbalance_alltimeseries{end+1} = metadata.stressBalance{iexp};
        end;

    end
end


%%
% Print number of members in a bin?

to_test = meltsens_alltimeseries;
bin = 'high melt sensitivity';
%bin = 'medium melt sensivity';
%bin = 'low melt sensitivity';

to_test = calving_alltimeseries;
bin = 'Strong calving';
%bin = 'Weak calving';
%bin = 'No calving';

bin_size = 0;
for iexp=1:length(to_test)
    if strcmp(convertCharsToStrings(to_test{iexp}), bin)
        bin_size = bin_size +1;
    end;
end;

bin, bin_size


%%

% Print number of members in a bin?

to_test1 = meltsens_alltimeseries;
bin1 = 'high melt sensitivity'; % Removed bin1 = 'medium melt sensivity';
%bin1 = 'low melt sensitivity';

to_test2 = calving_alltimeseries;
bin2 = 'Strong calving';
bin2 = 'Weak calving';
bin2 = 'No calving';

bin_size = 0;
for iexp=1:length(to_test)
    if strcmp(convertCharsToStrings(to_test1{iexp}), bin1) & strcmp(convertCharsToStrings(to_test2{iexp}), bin2)
        bin_size = bin_size +1;
    end;
end;

bin1, bin2, bin_size

%%

% 3-way ANOVA

% select groups for 3-way anova
g1_alltimeseries = climate_alltimeseries;
%g1_alltimeseries = init_alltimeseries;
g2_alltimeseries = meltsens_alltimeseries;
%g2_alltimeseries = meltparameterisation_alltimeseries; 
g3_alltimeseries = calving_alltimeseries; 
%g3_alltimeseries = ice_alltimeseries; 
%g3_alltimeseries = resolution_alltimeseries; 
%g3_alltimeseries = glresolution_alltimeseries; 
%g3_alltimeseries = subglmelt_alltimeseries; 
%g3_alltimeseries = init_alltimeseries;
%g3_alltimeseries = gia_alltimeseries; 
%g3_alltimeseries = stressbalance_alltimeseries;


%meltparameters_alltimeseries = {};

%titlestring = 'var: , g1: init method, g2: melt sens, g3: stress balance'
%titlestring = 'var: , g1: init method, g2: melt sens, g3: resolution'
%titlestring = 'var: , cum BMR: climate, g2: melt sens, g3: calving'
%titlestring = 'var: cum BMR, g1: climate, g2: melt sens, g3: ice model'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: resolution'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: gl resolution'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: subglmelt'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: melt parameterisation'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: init method'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: gia'
%titlestring = 'var: , g1: climate, g2: melt sens, g3: stress balance'
%titlestring = 'var: cum BMR, g1: climate, g2: melt param, g3: calving'
titlestring = 'var: cum BMR, g1: climate, g2: melt param, g3: ice model'



%Calculate variance and relative variance
var_ensemble=var(data_alltimeseries,1);
var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_g2=zeros(285,1);
var_g3=zeros(285,1);
var_g1g2=zeros(285,1);
var_g1g3=zeros(285,1);
var_g2g3=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285,
		individual_results=data_alltimeseries(:,yeari);
		[p,tabl,stats]=anovan(individual_results,{g1_alltimeseries,g2_alltimeseries,g3_alltimeseries},...
			'model',2,'display','off', 'sstype', 1) ;
		var_g1(yeari)=tabl{2,2}/length(g1_alltimeseries);
		var_g2(yeari)=tabl{3,2}/length(g1_alltimeseries);
		var_g3(yeari)=tabl{4,2}/length(g1_alltimeseries);
		var_g1g2(yeari)=tabl{5,2}/length(g1_alltimeseries);
		var_g1g3(yeari)=tabl{6,2}/length(g1_alltimeseries);
		var_g2g3(yeari)=tabl{7,2}/length(g1_alltimeseries);
		var_error(yeari)=tabl{8,2}/length(g1_alltimeseries);
		var_ensemble_anova(yeari)=tabl{end,2}/length(g1_alltimeseries);
end

figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_g2(2:end)),'color', [102 204 238]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g3(2:end)),'color', [204 187 68]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g2(2:end)),'color', [170 51 119]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g3(2:end)),'color', [34 136 51]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g2g3(2:end)),'color', [68 19 170]/256,'linewidth',2)
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
%ylabel('\sigma (m SLE)','fontsize',13)
ylabel('\sigma (Gt)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

title(sector)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_g2(2:end),var_g3(2:end),var_g1g2(2:end),var_g1g3(2:end),var_g2g3(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;
a(2).FaceColor=[102 204 238]/256;
a(3).FaceColor=[204 187 68]/256;
a(4).FaceColor=[170 51 119]/256;
a(5).FaceColor=[34 136 51]/256;
a(6).FaceColor=[68 19 170]/256;
a(7).FaceColor=[187 187 187]/256;
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA_withcollapse.pdf');


title(titlestring)


%%
% figure; hold all;
% for i=1:length(meltsens_alltimeseries_values)
%     scatter(i,meltsens_alltimeseries_values{i})
%     text(i,meltsens_alltimeseries_values{i}, meltparameterisation_alltimeseries{i})
% end

figure; hold all;
for i=1:length(meltsens_alltimeseries_values)
    scatter(i,data_alltimeseries(i,85))
    text(i,data_alltimeseries(i,85), meltsens_alltimeseries{i})
end

% 
% figure; hold all;
% for i=1:length(meltsens_alltimeseries_values)
%     plot(data_alltimeseries(i,:))
%     %text(i,data_alltimeseries(i,85), meltsens_alltimeseries{i})
% end

%%

%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensCalving.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_cumBMR_'+sector+'_ClimateMeltSensCalving.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_cumBMR_'+sector+'_ClimateMeltParameterisationCalving.pdf');

print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_cumBMR_'+sector+'_ClimateMeltParameterisationIcemodel.pdf');


%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensGLResolution.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensSubglmelt.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensMeltParameterisation.pdf');

%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_DSLC_ClimateMeltSensInitAllFactorlevels.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_DSLC_ClimateMeltSensInitSPvsDA.pdf');

%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_DSLC_ClimateMeltSensGIA.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_DSLC_ClimateMeltSensStressbalance.pdf');

%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_DSLC_InitiMethodMeltSensStressbalance.pdf');


%% test 1 way anova

g1_data = meltsens_alltimeseries;
%g1_data = calving_alltimeseries;
%g1_data = climate_alltimeseries;
%g1_data = ice_alltimeseries;

titlestring = 'var: sle, g1: melt sensitivity'
%titlestring = 'var: sle, g1: calving group'
%titlestring = 'var: sle, g1: climate model'
%titlestring = 'var: sle, g1: ice model'


var_ensemble=var(data_alltimeseries,1);

var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285
    individual_results=data_alltimeseries(:,yeari);
    
    [p,tabl]= anovan(individual_results, {g1_data},'display','off');
    var_g1(yeari)=tabl{2,2}/length(g1_data);
    var_error(yeari)=tabl{end-1,2}/length(g1_data);
    var_ensemble_anova(yeari)=tabl{end,2}/length(g1_data);
end

%var_g1/var_total*100

%
figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','error','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
ylabel('\sigma (m SLE)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;
a(2).FaceColor=[187 187 187]/256;

legend('g1','error','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');


title(titlestring)

%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA1_MeltSens.pdf');


%% Try 2 way anova

g1_data = climate_alltimeseries;
%g2_data = calving_alltimeseries;
g2_data = meltsens_alltimeseries;

%titlestring='var: sle, g1: climate, g2: calving';
titlestring='var: sle, g1: climate, g2: melt sens';


var_ensemble=var(data_alltimeseries,1);

var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_g2=zeros(285,1);
%var_g1g2=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285
    individual_results=data_alltimeseries(:,yeari);
    [p,tabl]= anovan(individual_results, {g1_data,g2_data},'display','off');
    var_g1(yeari)=tabl{2,2}/length(g1_data);
    var_g2(yeari)=tabl{3,2}/length(g1_data);
    %var_g1g2(yeari)=tabl{4,2}/length(g1_data);
    var_error(yeari)=tabl{end-1,2}/length(g1_data);
    var_ensemble_anova(yeari)=tabl{end,2}/length(g1_data);
end



% PLOT 2 way anova

figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_g2(2:end)),'color', [102 204 238]/256,'linewidth',2)
%plot(years(2:end), sqrt(var_g3(2:end)),'color', [204 187 68]/256,'linewidth',2)
%plot(years(2:end), sqrt(var_g1g2(2:end)),'color', [170 51 119]/256,'linewidth',2)
%plot(years(2:end), sqrt(var_g1g3(2:end)),'color', [34 136 51]/256,'linewidth',2)
%plot(years(2:end), sqrt(var_g2g3(2:end)),'color', [68 19 170]/256,'linewidth',2)
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','g2','2-way interaction','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
ylabel('\sigma (m SLE)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_g2(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;
a(2).FaceColor=[102 204 238]/256;
a(3).FaceColor=[187 187 187]/256;
legend('g1','g2','2-way interaction','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');

title(titlestring)

%%
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA2_SLE_ClimateCalving.pdf');

print(gcf, '-dpdf', '-painters', 'Figures/ANOVA2_SLE_ClimateMeltSens.pdf');



%% 3-way ANOVA

% select groups for 3-way anova
g1_alltimeseries = climate_alltimeseries;
g2_alltimeseries = meltsens_alltimeseries;
g3_alltimeseries = calving_alltimeseries; 
%g3_alltimeseries = ice_alltimeseries; 
%g3_alltimeseries = resolution_alltimeseries; 
%g3_alltimeseries = glresolution_alltimeseries; 
%g3_alltimeseries = subglmelt_alltimeseries; 
%g3_alltimeseries = meltparameterisation_alltimeseries; 

%meltparameters_alltimeseries = {};


titlestring = 'var: sle, g1: climate, g2: melt sens, g3: calving'
%titlestring = 'var: sle, g1: climate, g2: melt sens, g3: ice model'
%titlestring = 'var: sle, g1: climate, g2: melt sens, g3: resolution'
%titlestring = 'var: sle, g1: climate, g2: melt sens, g3: gl resolution'
%titlestring = 'var: sle, g1: climate, g2: melt sens, g3: subglmelt'
%titlestring = 'var: sle, g1: climate, g2: melt sens, g3: melt parameterisation'

%Calculate variance and relative variance
var_ensemble=var(data_alltimeseries,1);
var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_g2=zeros(285,1);
var_g3=zeros(285,1);
var_g1g2=zeros(285,1);
var_g1g3=zeros(285,1);
var_g2g3=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285,
		individual_results=data_alltimeseries(:,yeari);
		[p,tabl]=anovan(individual_results,{g1_alltimeseries,g2_alltimeseries,g3_alltimeseries},...
			'model',2,'display','off', 'sstype', 3) ;
		var_g1(yeari)=tabl{2,2}/length(g1_alltimeseries);
		var_g2(yeari)=tabl{3,2}/length(g1_alltimeseries);
		var_g3(yeari)=tabl{4,2}/length(g1_alltimeseries);
		var_g1g2(yeari)=tabl{5,2}/length(g1_alltimeseries);
		var_g1g3(yeari)=tabl{6,2}/length(g1_alltimeseries);
		var_g2g3(yeari)=tabl{7,2}/length(g1_alltimeseries);
		var_error(yeari)=tabl{8,2}/length(g1_alltimeseries);
		var_ensemble_anova(yeari)=tabl{end,2}/length(g1_alltimeseries);
end

figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_g2(2:end)),'color', [102 204 238]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g3(2:end)),'color', [204 187 68]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g2(2:end)),'color', [170 51 119]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g3(2:end)),'color', [34 136 51]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g2g3(2:end)),'color', [68 19 170]/256,'linewidth',2)
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
ylabel('\sigma (m SLE)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_g2(2:end),var_g3(2:end),var_g1g2(2:end),var_g1g3(2:end),var_g2g3(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;
a(2).FaceColor=[102 204 238]/256;
a(3).FaceColor=[204 187 68]/256;
a(4).FaceColor=[170 51 119]/256;
a(5).FaceColor=[34 136 51]/256;
a(6).FaceColor=[68 19 170]/256;
a(7).FaceColor=[187 187 187]/256;
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA_withcollapse.pdf');


title(titlestring)


%%

%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensCalving.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensIceModel.pdf');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensGLResolution.pdf');
print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensSubglmelt.pdf');
print(gcf, '-dpdf', '-painters', 'Figures/ANOVA3_SLE_ClimateMeltSensMeltParameterisation.pdf');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD


%% test 1 way anova

g1_data = meltsens_alltimeseries;
%g1_data = calving_alltimeseries;
g1_data = climate_alltimeseries;


var_ensemble=var(data_alltimeseries,1);

var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285
    individual_results=data_alltimeseries(:,yeari);
    
    [p,tabl]= anovan(individual_results, {g1_data},'display','off');
    var_g1(yeari)=tabl{2,2}/length(g1_data);
    var_error(yeari)=tabl{end-1,2}/length(g1_data);
    var_ensemble_anova(yeari)=tabl{end,2}/length(g1_data);
end

%var_g1/var_total*100

%%
figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','error','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
%ylabel('\sigma (m SLE)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;

legend('g1','error','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA_withcollapse.pdf');




%%
% select groups for 3-way anova
g1_alltimeseries = climate_alltimeseries;
g2_alltimeseries = calving_alltimeseries;
g3_alltimeseries = meltsens_alltimeseries; 

titlestring = 'var: melt, g1: climate, g2: calving type, g3: melt sens'

%Calculate variance and relative variance
var_ensemble=var(data_alltimeseries,1);
var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_g2=zeros(285,1);
var_g3=zeros(285,1);
var_g1g2=zeros(285,1);
var_g1g3=zeros(285,1);
var_g2g3=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285,
		individual_results=data_alltimeseries(:,yeari);
		[p,tabl]=anovan(individual_results,{g1_alltimeseries,g2_alltimeseries,g3_alltimeseries},...
			'model',2,'display','off') ;
		var_g1(yeari)=tabl{2,2}/length(g1_alltimeseries);
		var_g2(yeari)=tabl{3,2}/length(g1_alltimeseries);
		var_g3(yeari)=tabl{4,2}/length(g1_alltimeseries);
		var_g1g2(yeari)=tabl{5,2}/length(g1_alltimeseries);
		var_g1g3(yeari)=tabl{6,2}/length(g1_alltimeseries);
		var_g2g3(yeari)=tabl{7,2}/length(g1_alltimeseries);
		var_error(yeari)=tabl{8,2}/length(g1_alltimeseries);
		var_ensemble_anova(yeari)=tabl{end,2}/length(g1_alltimeseries);
end

figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_g2(2:end)),'color', [102 204 238]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g3(2:end)),'color', [204 187 68]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g2(2:end)),'color', [170 51 119]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g3(2:end)),'color', [34 136 51]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g2g3(2:end)),'color', [68 19 170]/256,'linewidth',2)
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
%ylabel('\sigma (m SLE)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_g2(2:end),var_g3(2:end),var_g1g2(2:end),var_g1g3(2:end),var_g2g3(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;
a(2).FaceColor=[102 204 238]/256;
a(3).FaceColor=[204 187 68]/256;
a(4).FaceColor=[170 51 119]/256;
a(5).FaceColor=[34 136 51]/256;
a(6).FaceColor=[68 19 170]/256;
a(7).FaceColor=[187 187 187]/256;
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA_withcollapse.pdf');


title(titlestring)




























%% OLD
%%
%Go over computed scalars and prepare the data
expnum=0;
for iexp=1:length(expnames),
		file_exp=dir([computed_name '/' expnames{iexp} '/' varn]);
		for ifile=length(file_exp):-1:1,
			if length(file_exp(ifile).name)<3,
				file_exp(ifile)=[];
			end
		end

		for ifile=1:length(file_exp),
			filename=[computed_name '/' expnames{iexp} '/' varn '/' file_exp(ifile).name];
			data=ncread(filename,varn);
			time=2016:2016+length(data)-1;
            % This makes sure to use only one run per group
            % FIXME: Do we want that?
			if (count(file_exp(ifile).name,'NORCE') & ~count(file_exp(ifile).name,'NORCE_CISM3-MAR364-ERA-t1-non')) | ...
					(count(file_exp(ifile).name,'VUW') & ~count(file_exp(ifile).name,'VUW_PISM1_e')) | ...
					(count(file_exp(ifile).name,'DOE') & ~count(file_exp(ifile).name,'DOE_MALI_4km')) | ...
					(count(file_exp(ifile).name,'NCAR') & ~count(file_exp(ifile).name,'NCAR_CISM2')) | ...
					(count(file_exp(ifile).name,'LSCE') & ~count(file_exp(ifile).name,'LSCE_GRISLI_e')) | ...
					(count(file_exp(ifile).name,'ULB') & ~count(file_exp(ifile).name,'ULB_fETISh-KoriBU1')) | ...
					(count(file_exp(ifile).name,'IMAU') & ~count(file_exp(ifile).name,'IMAU_UFEMISM1')) | ...
					(count(file_exp(ifile).name,'IGE_Elmer')) | ...
					(count(file_exp(ifile).name,'NCAR')) | ...
					(count(file_exp(ifile).name,'NORCE')) | ...
					(count(file_exp(ifile).name,'PIK')) | ...
					(count(file_exp(ifile).name,'UTAS')) | ...
					(count(file_exp(ifile).name,'VUB')) | ...
					(count(file_exp(ifile).name,'VUW')) | ...
					(count(file_exp(ifile).name,'DC_ISSM'))
				%do nothing not to bias ensemble
			else
				expnum = expnum+1;
				data_alltimeseries(expnum,:)=-(data(1:285)-data(1))/(362.5*1000);
				ice_alltimeseries{end+1}=file_exp(ifile).name(19:end-11);
				coll_alltimeseries{end+1}='no';
				if iexp==1,
					climate_alltimeseries{end+1}='CCSM4';
				elseif iexp==2;
					climate_alltimeseries{end+1}='HadGEM2';
				elseif iexp==3;
					climate_alltimeseries{end+1}='CESM2';
				elseif iexp==4;
					climate_alltimeseries{end+1}='UKESM';
				else 
					error('exp not supported yet');
                end
                
                % Define melt sensitivity group based on melt sens
                % load melt sensitivity
                ms = readtable("tables/linear_meltsensetivity.csv");
                i = (strcmp(ms.Model, file_exp(ifile).name) & strcmp(ms.Experiment,expnames(iexp))) ; 
                sens = ms.AIS(i);
                
                % FIXME set these better
                small_sens = 5;
                %large_sens = 5;

                if sens <= small_sens
                    sens_group='low melt sensitivity';
                %elseif sens <= large_sens
                %    sens_group='large melt sensivity';
                else 
                    sens_group='high melt sensitivity';
                end;
                meltsens_alltimeseries{end+1}=sens_group;
                
                
                % Define calving groups based on calving criertia

                i = (strcmp(metadata.Model, file_exp(ifile).name) & strcmp(metadata.Experiment,expnames(iexp))) ; 
                
                %calving_alltimeseries


%                 % Examples on how add group qualifiers based on the metadata 
%                 s = erase(erase(file_exp(ifile).name,"computed_"+varn+"_AIS_"), strcat("_",expnames(iexp),".nc" ));
%                 i = (strcmp(metadata.fileID, s) & strcmp(metadata.Experiment,expnames(iexp))) ; 
%                 s = metadata.MeltParameterisation(i);
%                 meltparameterisation_alltimeseries{end+1} = s{:} ;
%                 s = metadata.MeltParameters(i);
%                 meltparameters_alltimeseries{end+1} = s{:} ;
%                 s = metadata.ResolutionGL(i);
%                 glresolution_alltimeseries{end+1} = s{:}; 
%                 s = metadata.GLMelt(i);
%                 subglmelt_alltimeseries{end+1} = s{:}; 
%                 s = metadata.Resolution(i);
%                 resolution_alltimeseries{end+1} = s{:}; 
			end
		end
end

% %%
% %Add collapse experiments
% expnames={'expAE11','expAE12','expAE13','expAE14'};
% 
% %Go over computed scalars and prepare the data
% for iexp=1:length(expnames),
% 		file_exp=dir([computed_name '/' expnames{iexp} '/shelfmelt']);
% 		for ifile=length(file_exp):-1:1,
% 			if length(file_exp(ifile).name)<3,
% 				file_exp(ifile)=[];
% 			end
% 		end
% 
% 		for ifile=1:length(file_exp),
% 			filename=[computed_name '/' expnames{iexp} '/shelfmelt/' file_exp(ifile).name];
% 			data=ncread(filename,'shelfmelt');
% 			time=2016:2016+length(data)-1;
% 			if (count(file_exp(ifile).name,'NORCE') & ~count(file_exp(ifile).name,'NORCE_CISM3-MAR364-ERA-t1-non')) | ...
% 					(count(file_exp(ifile).name,'VUW') & ~count(file_exp(ifile).name,'VUW_PISM1_e')) | ...
% 					(count(file_exp(ifile).name,'DOE') & ~count(file_exp(ifile).name,'DOE_MALI_4km')) | ...
% 					(count(file_exp(ifile).name,'NCAR') & ~count(file_exp(ifile).name,'NCAR_CISM2')) | ...
% 					(count(file_exp(ifile).name,'LSCE') & ~count(file_exp(ifile).name,'LSCE_GRISLI_e')) | ...
% 					(count(file_exp(ifile).name,'ULB') & ~count(file_exp(ifile).name,'ULB_fETISh-KoriBU1')) | ...
% 					(count(file_exp(ifile).name,'IMAU') & ~count(file_exp(ifile).name,'IMAU_UFEMISM1')) | ...
% 					(count(file_exp(ifile).name,'IGE_Elmer')) | ...
% 					(count(file_exp(ifile).name,'NCAR')) | ...
% 					(count(file_exp(ifile).name,'NORCE')) | ...
% 					(count(file_exp(ifile).name,'PIK')) | ...
% 					(count(file_exp(ifile).name,'UTAS')) | ...
% 					(count(file_exp(ifile).name,'VUB')) | ...
% 					(count(file_exp(ifile).name,'VUW')) | ...
% 					(count(file_exp(ifile).name,'DC_ISSM'))
% 				%do nothing not to bias ensemble
% 			else
% 				expnum = expnum+1;
% 				data_alltimeseries(end+1,:)=-(data(1:285)-data(1))/(362.5*1000);
% 				ice_alltimeseries{end+1}=file_exp(ifile).name(19:end-11);
% 				coll_alltimeseries{end+1}='yes';
% 				if iexp==1,
% 					climate_alltimeseries{end+1}='CCSM4';
% 				elseif iexp==2;
% 					climate_alltimeseries{end+1}='HadGEM2';
% 				elseif iexp==3;
% 					climate_alltimeseries{end+1}='CESM2';
% 				elseif iexp==4;
% 					climate_alltimeseries{end+1}='UKESM';
% 				else 
% 					error('exp not supported yet');
%                 end
%                 s = erase(erase(file_exp(ifile).name,"computed_shelfmelt_AIS_"), strcat("_",expnames(iexp),".nc" ));
%                 i = (strcmp(metadata.fileID, s)); % & strcmp(metadata.Experiment,expnames(iexp))) ; 
%                 s = metadata.MeltParameterisation(i);
%                 meltparametersation_alltimeseries{end+1} = s{1} ;
% 
% 			end
% 		end
% end


%% test 1 way anova

for yeari=50:50 %285
    individual_results=data_alltimeseries(:,yeari);
    group = climate_alltimeseries;
    [p,tabl]= anova1(individual_results, group, 'display', 'off')
    var_g1=tabl{2,2}/length(meltsens_alltimeseries);
    var_total=tabl{end,2}/length(meltsens_alltimeseries);
end

var_g1/var_total*100

%%
for yeari=50:50 %285
    individual_results=data_alltimeseries(:,yeari);
    group = meltsens_alltimeseries;
    [p,tabl]= anova1(individual_results, group, 'display', 'off')
    var_g1=tabl{2,2}/length(meltsens_alltimeseries);
    var_total=tabl{end,2}/length(meltsens_alltimeseries);
end

var_g1/var_total*100

%% FIXME: Try two way anova with both, but really, we want to add the calving group!


%%
% select groups for 3-way anova
g1_alltimeseries = climate_alltimeseries;
g2_alltimeseries = meltsens_alltimeseries;
g3_alltimeseries = ice_alltimeseries;

titlestring = 'var: melt, g1: climate, g2: melt sens, g3: ice model'



%Calculate variance and relative variance
var_ensemble=var(data_alltimeseries,1);
var_ensemble_anova=zeros(285,1);
var_g1=zeros(285,1);
var_g2=zeros(285,1);
var_g3=zeros(285,1);
var_g1g2=zeros(285,1);
var_g1g3=zeros(285,1);
var_g2g3=zeros(285,1);
var_error=zeros(285,1);

for yeari=2:285,
		individual_results=data_alltimeseries(:,yeari);
		[p,tabl]=anovan(individual_results,{g1_alltimeseries,g2_alltimeseries,g3_alltimeseries},...
			'model',2,'display','off') ;
		var_g1(yeari)=tabl{2,2}/length(g1_alltimeseries);
		var_g2(yeari)=tabl{3,2}/length(g1_alltimeseries);
		var_g3(yeari)=tabl{4,2}/length(g1_alltimeseries);
		var_g1g2(yeari)=tabl{5,2}/length(g1_alltimeseries);
		var_g1g3(yeari)=tabl{6,2}/length(g1_alltimeseries);
		var_g2g3(yeari)=tabl{7,2}/length(g1_alltimeseries);
		var_error(yeari)=tabl{8,2}/length(g1_alltimeseries);
		var_ensemble_anova(yeari)=tabl{end,2}/length(g1_alltimeseries);
end

figure(1); clf; set(gcf,'color','w'); hold on;
set(gcf,'position',[100 200 1200 920]);


subplot(2,1,1)
years=(1:285)+2015;
plot(years(2:end), sqrt(var_g1(2:end)),'color', [238 102 119]/256,'linewidth',2)
hold on
plot(years(2:end), sqrt(var_g2(2:end)),'color', [102 204 238]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g3(2:end)),'color', [204 187 68]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g2(2:end)),'color', [170 51 119]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g1g3(2:end)),'color', [34 136 51]/256,'linewidth',2)
plot(years(2:end), sqrt(var_g2g3(2:end)),'color', [68 19 170]/256,'linewidth',2)
plot(years(2:end), sqrt(var_error(2:end)),'color', [187 187 187]/256,'linewidth',2)
plot(years(2:end), sqrt(var_ensemble_anova(2:end)),'color', [0 0 0]/256,'linewidth',2)
hold off
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
%ylim([0 1.8])
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','total','location','West','fontsize',13)
xlabel('Year','fontsize',13)
%ylabel('\sigma (m SLE)','fontsize',13)
text(2005,0,'a','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)
%
subplot(2,1,2)
years=(1:285)+2015;
a=area(years(2:end),([var_g1(2:end),var_g2(2:end),var_g3(2:end),var_g1g2(2:end),var_g1g3(2:end),var_g2g3(2:end),var_error(2:end)]./var_ensemble_anova(2:end)*100));
a(1).FaceColor=[238 102 119]/256;
a(2).FaceColor=[102 204 238]/256;
a(3).FaceColor=[204 187 68]/256;
a(4).FaceColor=[170 51 119]/256;
a(5).FaceColor=[34 136 51]/256;
a(6).FaceColor=[68 19 170]/256;
a(7).FaceColor=[187 187 187]/256;
legend('g1','g2','g3','g1-g2','g1-g3','g2-g3','3-way interaction','location','West','fontsize',13)
grid on
ax = gca;
ax.Layer = 'top';
xlim([2017 2300])
ylim([1 100])
xlabel('Year','fontsize',13)
ylabel('Percentage of variance','fontsize',13)
text(2005,0,'b','VerticalAlignment','middle','HorizontalAlignment','right','fontsize',16,'fontweight','b');
set(gca,'fontsize',13)

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
pos = get(h,'Position');
%print(gcf, '-dpdf', '-painters', 'Figures/ANOVA_withcollapse.pdf');


title(titlestring)



