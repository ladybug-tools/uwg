% Read Energy Plus Output File (based on entire year)

clear;
close all;
file = {'RefBldgFullServiceRestaurantPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgHospitalPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgLargeHotelPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgLargeOfficePost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgMediumOfficePost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgMidriseApartmentPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgOutPatientPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgPrimarySchoolPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgQuickServiceRestaurantPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgSecondarySchoolPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgSmallHotelPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgSmallOfficePost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgStand-aloneRetailPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgStripMallPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgSuperMarketPost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx';
    'RefBldgWarehousePost1980_v1.4_7.2_5A_USA_IL_CHICAGO-OHARE.xlsx'};
FLarea = [511,22422,11345,46320,4982,3135,3089,6871,232,19592,4014,511,2294,2090,4181,4835];
ymax = [200,100,100,100,100,100,100,100,200,100,100,100,100,100,100,100];
BLD = {'FullServiceRestaurant';
    'Hospital';
    'LargeHotel';
    'LargeOffice';
    'MedOffice';
    'MidRiseApartment';
    'OutPatient';
    'PrimarySchool';
    'QuickServiceRestaurant';
    'SecondarySchool';
    'SmallHotel';
    'SmallOffice';
    'StandAloneRetail';
    'StripMall';
    'SuperMarket';
    'WareHouse'};

hour24 = 1:1:24;
month = {' (Jan)',' (Jul)'};
rowstart = [1,4345];
rowend = [744,5088];
NRGcool = zeros(16,1);
NRGheat = zeros(16,1);

% Building no.
for i = 16:16
%    input = strcat('Results\EPEnergyUseBostonRural\',file{i});
    input = strcat('Results\EPEnergyUseBostonUrbanMidRise\',file{i});

    [~,txt,~] = xlsread(input,1,'B1:R1');
    [num,~,~] = xlsread(input,1);   
    num = num/3600;
    
    % Jan (k = 1) or July (k = 2)
    for k = 1:1
        % change number from J/hr to W
        ElecTotal = zeros(24,1);
        ElecLight = zeros(24,1);
        ElecEquip = zeros(24,1);
        ElecCool = zeros(24,1);
        ElecFan = zeros(24,1);
        ElecPump = zeros(24,1);
        ElecRefg = zeros(24,1);
        ElecHum = zeros(24,1);
        GasTotal = zeros(24,1);
        GasHeat = zeros(24,1);
        GasEquip = zeros(24,1);
        GasWater = zeros(24,1);

        for col = 1:numel(txt)
            if strcmp(txt(col),'Electricity:Facility [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end
                    ElecTotal(hour) = ElecTotal(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'Cooling:Electricity [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end
                    ElecCool(hour) = ElecCool(hour) + num(j,col)/31;
                end
                NRGcool(i) = sum(num(:,col))/1000; 
            end
            if strcmp(txt(col),'InteriorLights:Electricity [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    ElecLight(hour) = ElecLight(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'InteriorEquipment:Electricity [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    ElecEquip(hour) = ElecEquip(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'Fans:Electricity [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    ElecFan(hour) = ElecFan(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'Pumps:Electricity [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    ElecPump(hour) = ElecPump(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'Refrigeration:Electricity [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    ElecRefg(hour) = ElecRefg(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'Gas:Facility [J](Hourly)')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    GasTotal(hour) = GasTotal(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'Heating:Gas [J](Hourly)') || strcmp(txt(col),'Heating:Gas [J](Hourly) ')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    GasHeat(hour) = GasHeat(hour) + num(j,col)/31;
                end
                NRGheat(i) = sum(num(:,col))/1000; 
            end
            if strcmp(txt(col),'InteriorEquipment:Gas [J](Hourly)') || strcmp(txt(col),'InteriorEquipment:Gas [J](Hourly) ')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    GasEquip(hour) = GasEquip(hour) + num(j,col)/31;
                end
            end
            if strcmp(txt(col),'WaterSystems:Gas [J](Hourly)') || strcmp(txt(col),'WaterSystems:Gas [J](Hourly) ')
                for j = rowstart(k):rowend(k)
                    hour = mod(j,24);
                    if hour == 0
                        hour = 24;
                    end                    
                    GasWater(hour) = GasWater(hour) + num(j,col)/31;
                end
            end
        end
%         figure;
%         Elec = [ElecCool,ElecLight,ElecEquip,ElecTotal-ElecCool-ElecLight-ElecEquip]/FLarea(i);
%         area(hour24,Elec);
%         ax = gca;
%         ax.XTick = [0 6 12 18 24];
%         
%         grid;
%         grid minor;
%         title (strcat(BLD(i),': Electricity Usage ',month(k),' - EnergyPlus'));
%         ylim = [0 20];
%         xlabel ('Hour'); 
%         ylabel ('W/m^2');
%         legend ('Cool','Light','Equip','Misc','Location','NorthWest');
%         
        figure;
        Gas = [GasHeat, GasWater, GasEquip]/FLarea(i);
        area(hour24,Gas);
        ax = gca;
        ax.XTick = [0 6 12 18 24];
        grid;
        grid minor;
        title (strcat(BLD(i),': Gas Usage ',month(k),' - Urban'));
        xlabel ('Hour');
        ylabel ('W/m^2');
        legend ('Heat','Water','Equip','Location','NorthWest');
    end
end