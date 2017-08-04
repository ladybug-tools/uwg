
"""
Translated from: https://github.com/hansukyang/UWG_Matlab/blob/master/readDOE.m
Translated to Python by Saeran Vasanthakumar (saeranv@gmail.com) - April, 2017
"""

import csv

class UWG_Unit_Test:
    def __init__(self):
        self.fail = 0
        self.success = 0
        self.total = self._get_total()
        self.test_history = "Test Results:\n"
    def __repr__(self):
        return "{a} successful and {b} failed tests".format(a=self.success,b=self.fail)
    def test_results(self):
        return self.test_history + self.__repr__()
    def _get_total(self):
        return self.fail + self.success
    def test_equality(self,a,b,toggle=True):
        if a == b:
            s = "test_equality: {y} == {z} success\n".format(y=a,z=b)
            self.success += 1
        else:
            s = "test_equality: {y} != {z} fail\n".format(y=a,z=b)
            self.fail += 1
        if toggle:
            self.test_history += s
    def test_in_string(self,a,b,toggle=True):
        if type(a)!=type("") or type(b)!=type(""):
            s = "test_in_string: {y} or {z} not a string\n".format(y=b,z=a)
            self.fail += 1
        else:
            if b in a:
                s = "test_in_string:: {y} in {z} success\n".format(y=b,z=a)
                self.success += 1
            else:
                s = "test_in_string: {y} in {z} fail\n".format(y=b,z=a)
                self.fail += 1
        if toggle:
            self.test_history += s


def readDOE():
    """
    Read Excel files of DOE buildings
    Sheet 1 = BuildingSummary
    Sheet 2 = ZoneSummary
    Sheet 3 = LocationSummary
    Sheet 4 = Schedules
    Note BLD8 & 10 = school

    args:
        ...
    returns:
        ...
    """

    #Make a test object
    test_readDOE = UWG_Unit_Test()
    test_readDOE.toggle = True

    # DOE Building Types
    bldType = [
        'FullServiceRestaurant',
        'Hospital',
        'LargeHotel',
        'LargeOffice',
        'MedOffice',
        'MidRiseApartment',
        'OutPatient',
        'PrimarySchool',
        'QuickServiceRestaurant',
        'SecondarySchool',
        'SmallHotel',
        'SmallOffice',
        'StandAloneRetail',
        'StripMall',
        'SuperMarket',
        'WareHouse']

    zoneType = [
        '1A (Miami)',
        '2A (Houston)',
        '2B (Phoenix)',
        '3A (Atlanta)',
        '3B-CA (Los Angeles)',
        '3B (Las Vegas)',
        '3C (San Francisco)',
        '4A (Baltimore)',
        '4B (Albuquerque)',
        '4C (Seattle)',
        '5A (Chicago)',
        '5B (Boulder)',
        '6A (Minneapolis)',
        '6B (Helena)',
        '7 (Duluth)',
        '8 (Fairbanks)']

    builtEra = ['Pre80',
        'Pst80',
        'New']


    # Building (Type,Era,Zone), Type = 1-16, Era = 1-3, Zone = 1-16
    #refDOE (16,3,16) = Building;
    #Schedule (16,3,16) = SchDef;
    #refBEM (16,3,16) = BEMDef;

    def read_doe_csv(file_doe_name_):
        file_doe = open(file_doe_name_,"r")
        gen_doe = csv.reader(file_doe, delimiter=",")
        list_doe = map(lambda r: r,gen_doe)
        file_doe.close()
        return list_doe

    #Purpose: Loop through every DOE reference csv and extract building data
    #Nested loop = 16 types, 3 era, 16 zones
    #Therefore time complexity O(n*m*k) = 768
    dir_doe_name = "DOERefBuildings"
    for i in xrange(1):#16
        print "BLD", i+1
        # Read building summary (Sheet 1)
        file_doe_name_bld = "{x}\\BLD{y}\\BLD{y}_BuildingSummary.csv".format(x=dir_doe_name,y=i+1)
        list_doe1 = read_doe_csv(file_doe_name_bld)

        nFloor      = map(lambda n: float(n), list_doe1[3][3:])      # Number of Floors
        glazing     = map(lambda n: float(n), list_doe1[4][3:])      # [?] Total
        hCeiling    = map(lambda n: float(n), list_doe1[5][3:])      # [m] Ceiling height
        ver2hor     = map(lambda n: float(n), list_doe1[7][3:])      # Wall to Skin Ratio
        AreaRoof    = map(lambda n: float(n), list_doe1[8][3:])      # [m2] Gross Dimensions - Total area

        #Tests for sheet 1
        test_readDOE.test_in_string(list_doe1[0][1],"Building Summary")

        # Read zone summary (Sheet 2)
        file_doe_name_zone = "{x}\\BLD{y}\\BLD{y}_ZoneSummary.csv".format(x=dir_doe_name,y=i+1)
        list_doe2 = read_doe_csv(file_doe_name_zone)

        AreaFloor   = map(lambda n: float(n), [list_doe2[2][5],list_doe2[3][5],list_doe2[4][5]])       # [m2]
        Volume      = map(lambda n: float(n), [list_doe2[2][6],list_doe2[3][6],list_doe2[4][6]])       # [m3]
        AreaWall    = map(lambda n: float(n), [list_doe2[2][8],list_doe2[3][8],list_doe2[4][8]])       # [m2]
        AreaWindow  = map(lambda n: float(n), [list_doe2[2][9],list_doe2[3][9],list_doe2[4][9]])       # [m2]
        Occupant    = map(lambda n: float(n), [list_doe2[2][11],list_doe2[3][11],list_doe2[4][11]])    # Number of People
        Lights      = map(lambda n: float(n), [list_doe2[2][12],list_doe2[3][12],list_doe2[4][12]])    # [W/m2]
        Elec        = map(lambda n: float(n), [list_doe2[2][13],list_doe2[3][13],list_doe2[4][13]])    # [W/m2] Electric Plug and Process
        Gas         = map(lambda n: float(n), [list_doe2[2][14],list_doe2[3][14],list_doe2[4][14]])    # [W/m2] Gas Plug and Process
        SHW         = map(lambda n: float(n), [list_doe2[2][15],list_doe2[3][15],list_doe2[4][15]])    # [Litres/hr] Peak Service Hot Water
        Vent        = map(lambda n: float(n), [list_doe2[2][17],list_doe2[3][17],list_doe2[4][17]])    # [L/s/m2] Ventilation
        Infil       = map(lambda n: float(n), [list_doe2[2][20],list_doe2[3][20],list_doe2[4][20]])    # Air Changes Per Hour (ACH) Infiltration

        #Tests sheet 2
        test_readDOE.test_equality(list_doe2[0][2],"Zone Summary")

        # Read location summary (Sheet 3)
        file_doe_name_location = "{x}\\BLD{y}\\BLD{y}_LocationSummary.csv".format(x=dir_doe_name,y=i+1)
        list_doe3 = read_doe_csv(file_doe_name_location)


        if not test_readDOE.toggle:
            for i,row in enumerate(list_doe3):
                print 'ROW {:d}:'.format(i),
                for j,col in enumerate(row):
                    print '[{:d}]'.format(j), col, ",",
                print ''
            print '--'


        #TypeWall =  [3][]
        """
        % Read location summary (Sheet 3)
        [num, text, ~] = xlsread(file,3);
        TypeWall = [text(4,5:20); text(15,5:20); text(26,5:20)];
        RvalWall = [num(5,5:20); num(16,5:20); num(27,5:20)];
        TypeRoof = [text(6,5:20); text(17,5:20); text(28,5:20)];
        RvalRoof = [num(7,5:20); num(18,5:20); num(29,5:20)];
        Uwindow = [num(8,5:20); num(19,5:20); num(30,5:20)];
        SHGC = [num(9,5:20); num(20,5:20); num(31,5:20)];
        HVAC = [num(10,5:20); num(21,5:20); num(32,5:20)];
        HEAT = [num(11,5:20); num(22,5:20); num(33,5:20)];
        COP = [num(12,5:20); num(23,5:20); num(34,5:20)];
        EffHeat = [num(13,5:20); num(24,5:20); num(35,5:20)];
        FanFlow = [num(14,5:20); num(25,5:20); num(36,5:20)];
        """
        #Test sheet 3
        test_readDOE.test_equality(list_doe3[0][2],"Location Summary")

        """
        % Read schedule (Sheet 4)
        [num, text, ~] = xlsread(file,4);
        SchEquip = num(2:4,7:30);
        SchLight = num(5:7,7:30);
        SchOcc = num(8:10,7:30);
        SetCool = num(11:13,7:30);
        SetHeat = num(14:16,7:30);
        SchGas = num(17:19,7:30);
        SchSWH = num(20:22,7:30);

        for j = 1:3
            for k = 1:16
                refDOE (i,j,k) = Building(hCeiling(j),...  % floorHeight
                    1,...                           % intHeatNight
                    1,...                           % intHeatDay
                    0.1,...                         % intHeatFRad
                    0.1,...                         % intHeatFLat
                    Infil(j),...                    % infil (ACH)
                    Vent(j)/1000,...                % vent (m^3/s/m^2)
                    glazing(j),...                  % glazing ratio
                    Uwindow(j,k),...                % uValue
                    SHGC(j,k),...                   % shgc
                    'AIR',...                       % a/c type
                    COP(j,k),...                    % cop
                    297,...                         % coolSetpointDay
                    297,...                         % coolSetpointNight
                    293,...                         % heatSetpointDay
                    293,...                         % heatSetpointNight
                    HVAC(j,k)*1000/AreaFloor(j),... % coolCap (refDOE in kW)
                    EffHeat(j,k),...                % heatEff
                    293);                           % initialTemp

                refDOE(i,j,k).heatCap = HEAT(j,k)*1000/AreaFloor(j);
                refDOE(i,j,k).Type = bldType(i);
                refDOE(i,j,k).Era = builtEra(j);
                refDOE(i,j,k).Zone = zoneType(k);

                % Define wall, roof, road, and mass
                % Material (thermalCond, volHeat = specific heat * density);
                Stucco = Material (0.6918, 837.0 * 1858.0);
                Concrete = Material (1.311, 836.8 * 2240);
                Insulation = Material (0.049, 836.8 * 265.0);
                Gypsum = Material (0.16, 830.0 * 784.9);
                Wood = Material (0.11,1210.0*544.62);

                % Wall (1 in stucco, concrete, insulation, gypsum)
                if strcmp(TypeWall(j,k),'MassWall')
                    % 1" stucco, 8" concrete, tbd insulation, 1/2" gypsum
                    Rbase = 0.271087; % based on stucco, concrete, gypsum
                    Rins = RvalWall(j,k) - Rbase;
                    D_ins = Rins * Insulation.thermalCond;
                    if D_ins > 0.01
                        thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;D_ins;0.0127];
                        layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
                    else
                        thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;0.0127];
                        layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Gypsum];
                    end

                    wall = Element(0.08,0.92,thickness,layers,0,293,0);

                    % If mass wall, assume mass foor (4" concrete)
                    % Mass (assume 4" concrete);
                    alb = 0.2;
                    emis = 0.9;
                    thickness = [0.054;0.054];
                    concrete = Material (1.31,2240.0*836.8);
                    mass = Element(alb,emis,thickness,[concrete;concrete],0,293,1);

                elseif strcmp(TypeWall(j,k),'WoodFrame')
                    % 0.01m wood siding, tbd insulation, 1/2" gypsum
                    Rbase = 0.170284091;    % based on wood siding, gypsum
                    Rins = RvalWall(j,k) - Rbase;
                    D_ins = Rins * Insulation.thermalCond;

                    if D_ins > 0.01
                        thickness = [0.01;D_ins;0.0127];
                        layers = [Wood;Insulation;Gypsum];
                    else
                        thickness = [0.01;0.0127];
                        layers = [Wood;Gypsum];
                    end
                    wall = Element(0.22,0.92,thickness,layers,0,293,0);

                    % If wood frame wall, assume wooden floor
                    alb = 0.2;
                    emis = 0.9;
                    thickness = 0.05 * ones (2,1);
                    wood = Material (1.31,2240.0*836.8);
                    mass = Element(alb,emis,thickness,[wood;wood],0,293,1);

                elseif strcmp(TypeWall(j,k),'SteelFrame')
                    % 1" stucco, 8" concrete, tbd insulation, 1/2" gypsum
                    Rbase = 0.271087; % based on stucco, concrete, gypsum
                    Rins = RvalWall(j,k) - Rbase;
                    D_ins = Rins * Insulation.thermalCond;
                    if D_ins > 0.01
                        thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;D_ins;0.0127];
                        layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
                    else    % If insulation is too thin, assume no insulation
                        thickness = [0.0254;0.0508;0.0508;0.0508;0.0508;0.0127];
                        layers = [Stucco;Concrete;Concrete;Concrete;Concrete;Gypsum];
                    end
                    wall = Element(0.15,0.92,thickness,layers,0,293,0);

    %                 wall = Element(0.2,0.92,thickness,layers,0,293,0); % Singapore case

    %                 % TOULOUS case
    %                 thickness = [0.05;0.05;0.05;0.05;0.05;0.05;0.03;0.0127];
    %                 layers = [Concrete;Concrete;Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
    %                 wall = Element(0.25,0.93,thickness,layers,0,293,0);
    %
    %                 % BUBBLE case
    %                 thickness = [0.05;0.05;0.05;0.05;0.03;0.0127];
    %                 layers = [Concrete;Concrete;Concrete;Concrete;Insulation;Gypsum];
    %                 wall = Element(0.15,0.93,thickness,layers,0,293,0);

                    % If mass wall, assume mass foor
                    % Mass (assume 4" concrete);
                    alb = 0.2;      % Adjusted for Bubble/Toulouse case (0.08 per Energy Plus)
                    emis = 0.93;
                    thickness = 0.05 * ones (2,1);
    %                thickness = 0.102*ones(2,1);       % for Bubble
                    mass = Element(alb,emis,thickness,[Concrete;Concrete],0,293,1);

                elseif strcmp(TypeWall(j,k),'MetalWall')

                    % metal siding, insulation, 1/2" gypsum
                    alb = 0.2;
                    emis = 0.9;
                    D_ins = max((RvalWall(j,k) * Insulation.thermalCond)/2,0.01);
                    wall = Element(alb,emis,[D_ins;D_ins;0.0127],[Insulation;Insulation;Gypsum],0,293,0);

                    % Mass (assume 4" concrete);
                    alb = 0.2;
                    emis = 0.9;
                    thickness = 0.05 * ones (2,1);
                    concrete = Material (1.31,2240.0*836.8);
                    mass = Element(alb,emis,thickness,[concrete;concrete],0,293,1);

                end

                % Roof
                if strcmp(TypeRoof(j,k),'IEAD')
                    % IEAD-> membrane, insulation, decking
                     alb = 0.2;
                     emis = 0.93;
                     D_ins = max(RvalRoof(j,k) * Insulation.thermalCond / 2,0.01);
                     roof = Element(alb,emis,[D_ins;D_ins],[Insulation;Insulation],0,293,0);

    %                 % TOULOUS/BUBBLE Case
    %                 alb = 0.15;
    %                 emis = 0.93;
    %                 thickness = [0.06;0.05;0.05;0.05;0.05;0.03];
    %                 roof = Element(alb,emis,thickness,[Concrete;Wood;Wood;Wood;Wood;Insulation],0,293,0);

                elseif strcmp(TypeRoof(j,k),'Attic')
                    % IEAD-> membrane, insulation, decking
                    alb = 0.2;
                    emis = 0.9;
                    D_ins = max(RvalRoof(j,k) * Insulation.thermalCond/2,0.01);
                    roof = Element(alb,emis,[D_ins;D_ins],[Insulation;Insulation],0,293,0);

                elseif strcmp(TypeRoof(j,k),'MetalRoof')
                    % IEAD-> membrane, insulation, decking
                    alb = 0.2;
                    emis = 0.9;
                    D_ins = max(RvalRoof(j,k) * Insulation.thermalCond/2,0.01);
                    roof = Element(alb,emis,[D_ins;D_ins],[Insulation;Insulation],0,293,0);
                end

                % Define bulding energy model, set fraction to zero
                refBEM(i,j,k) = BEMDef(refDOE(i,j,k),mass,wall,roof,0);
                refBEM(i,j,k).building.FanMax = FanFlow(j,k);

                Schedule(i,j,k).Elec = SchEquip;   % 3x24 matrix of schedule for electricity (WD,Sat,Sun)
                Schedule(i,j,k).Light = SchLight;  % 3x24 matrix of schedule for light (WD,Sat,Sun)
                Schedule(i,j,k).Gas = SchGas;      % 3x24 matrix of schedule for gas (WD,Sat,Sun)
                Schedule(i,j,k).Occ = SchOcc;      % 3x24 matrix of schedule for occupancy (WD,Sat,Sun)
                Schedule(i,j,k).Cool = SetCool;    % 3x24 matrix of schedule for cooling temp (WD,Sat,Sun) % off for BUBBLE case
                Schedule(i,j,k).Heat = SetHeat;    % 3x24 matrix of schedule for heating temp (WD,Sat,Sun)
                Schedule(i,j,k).SWH = SchSWH;      % 3x24 matrix of schedule for SWH (WD,Sat,Sun)

                Schedule(i,j,k).Qelec = Elec(j);         % W/m^2 (max) for electrical plug process
                Schedule(i,j,k).Qlight = Light(j);       % W/m^2 (max) for light
                Schedule(i,j,k).Nocc = Occupant(j)/AreaFloor(j); % Person/m^2
                Schedule(i,j,k).Qgas = Gas(j);           % W/m^2 (max) for gas
                Schedule(i,j,k).Vent = Vent(j)/1000;     % m^3/m^2 per person
                Schedule(i,j,k).Vswh = SHW(j)/AreaFloor(j);    % litres per hour per m^2 of floor

            end
        end
    end

    % % BUBBLE/TOULOUSE adjustment Case
    % refBEM(6,2,5).building.glazingRatio = 0.3;

    % % Singapore adjustment Case
    % refBEM(6,2,1).building.glazingRatio = 0.3;

    save ('RefDOE.mat','refDOE','refBEM','Schedule');

    % wall & roof definition based on material
    % Material,
    %     1/2IN Gypsum,            !- Name
    %     Smooth,                  !- Roughness
    %     0.0127,                  !- Thickness {m}
    %     0.1600,                  !- Conductivity {W/m-K}
    %     784.9000,                !- Density {kg/m3}
    %     830.0000,                !- Specific Heat {J/kg-K}
    %     0.9000,                  !- Thermal Absorptance
    %     0.9200,                  !- Solar Absorptance
    %     0.9200;                  !- Visible Absorptance
    %
    % Material,
    %     1IN Stucco,              !- Name
    %     Smooth,                  !- Roughness
    %     0.0253,                  !- Thickness
    %     0.6918,                  !- Conductivity
    %     1858.0000,               !- Density
    %     837.0000,                !- Specific Heat
    %     0.9000,                  !- Thermal Absorptance
    %     0.9200,                  !- Solar Absorptance
    %     0.9200;                  !- Visible Absorptance
    %
    % Material,
    %     8IN CONCRETE HW,  !- Name
    %     Rough,                   !- Roughness
    %     0.2032,                  !- Thickness {m}
    %     1.3110,                  !- Conductivity {W/m-K}
    %     2240.0000,               !- Density {kg/m3}
    %     836.8000,                !- Specific Heat {J/kg-K}
    %     0.9000,                  !- Thermal Absorptance
    %     0.7000,                  !- Solar Absorptance
    %     0.7000;                  !- Visible Absorptance
    %
    % Material,
    %     Mass NonRes Wall Insulation, !- Name
    %     MediumRough,             !- Roughness
    %     0.0484268844343858,      !- Thickness {m}
    %     0.049,                   !- Conductivity {W/m-K}
    %     265.0000,                !- Density {kg/m3}
    %     836.8000,                !- Specific Heat {J/kg-K}
    %     0.9000,                  !- Thermal Absorptance
    %     0.7000,                  !- Solar Absorptance
    %     0.7000;                  !- Visible Absorptance
    %
    % Material,
    %     Std Wood 6inch,          !- Name
    %     MediumSmooth,            !- Roughness
    %     0.15,                    !- Thickness {m}
    %     0.12,                    !- Conductivity {W/m-K}
    %     540.0000,                !- Density {kg/m3}
    %     1210,                    !- Specific Heat {J/kg-K}
    %     0.9000000,               !- Thermal Absorptance
    %     0.7000000,               !- Solar Absorptance
    %     0.7000000;               !- Visible Absorptance! Common Materials
    %
    % Material,
    %     Wood Siding,             !- Name
    %     MediumSmooth,            !- Roughness
    %     0.0100,                  !- Thickness {m}
    %     0.1100,                  !- Conductivity {W/m-K}
    %     544.6200,                !- Density {kg/m3}
    %     1210.0000,               !- Specific Heat {J/kg-K}
    %     0.9000,                  !- Thermal Absorptance
    %     0.7800,                  !- Solar Absorptance
    %     0.7800;                  !- Visible Absorptance
    """

    print test_readDOE.test_results()


if __name__ == "__main__":
    readDOE()
