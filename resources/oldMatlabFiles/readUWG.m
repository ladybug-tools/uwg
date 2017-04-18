% Read Energy Plus Output File (based on entire year)

clear;
file = {'Results\Midrise July - UWGoutput - v5.xlsx'};
%file = {'Results\Midrise Jan - UWGoutput - v6.xlsx'};

%file = {'UWGoutput.xlsx'};

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

hour = 1:1:24;

% Building no.
%input = strcat('output\',file{1});
input = file{1};
month = {' (Jan)',' (Jul)'};
k =2;


for i = 15:15

   [~,txt,~] = xlsread(input,3,'B1:R1');
    [num,~,~] = xlsread(input,2+i,'A2:V25');   

    ElecLight = num(:,10);
    ElecEquip = num(:,9);
    ElecCool = num(:,19);
    Gas = num(:,8);

    figure;
    Elec = [ElecCool,ElecLight,ElecEquip,ElecCool*0];
    area(hour,Elec);
    ax = gca;
    ax.XTick = [0 6 12 18 24];

    grid;
    grid minor;
    title (strcat(BLD(i),': Electricity Usage ',month(k),' - UWG'));
    xlabel ('Hour'); 
    ylabel ('W/m^2');
    legend ('Cool','Light','Equip','Location','NorthWest');
        
%         figure;
%         area(hour,Gas);
%         ax = gca;
%         ax.XTick = [0 6 12 18 24];
%         grid;
%         grid minor;
%         title (strcat(BLD(i),': Gas Usage ',month(k),' - UWG'));
%         xlabel ('Hour');
%         ylabel ('W/m^2');

end