from utilities import read_csv

class Weather(object):
    """
    Weather
    Read a epw file

    properties
        staTemp   % air temperature (C)
        staRhum   % air relative humidity (%)
        staPres   % air pressure (Pa)
        staInfra  % horizontal Infrared Radiation Intensity (W m-2)
        staHor    % horizontal radiation
        staDir    % normal solar direct radiation (W m-2)
        staDif    % horizontal solar diffuse radiation (W m-2)
        staUdir   % wind direction ()
        staUmod   % wind speed (m s-1)
        staRobs   % Precipitation (mm h-1)
        staHum    % specific humidty (kg kg-1)
    """

    DIR_EPW_NAME = "data\\epw\\"
    def __init__(self,climate_file,HI,HF):
        #print self.DIR_EPW_NAME + climate_file
        self.climate_data = read_csv(self.DIR_EPW_NAME + climate_file)
        self.location = self.climate_data[0][1]
        print 'HI: ', HI #"July 30th, timestep = hourly"
        print 'HF:', HF #7X24 = 168
        #H1 and HF define the row we want
        #HI = Initial sensor data index for row
        #H2 = final sensor data index for row
        #M = csvread(filename,R1,C1,[R1 C1 R2 C2]) reads only the range bounded by row offsets R1 and R2 and column offsets C1 and C2.
        #http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
        #Date,HH:MM,Datasource,DryBulb {C},DewPoint {C},RelHum {%},Atmos Pressure {Pa},
        #ExtHorzRad {Wh/m2},ExtDirRad {Wh/m2},HorzIRSky {Wh/m2},GloHorzRad {Wh/m2},
        #DirNormRad {Wh/m2},DifHorzRad {Wh/m2},GloHorzIllum {lux},DirNormIllum {lux},DifHorzIllum {lux},ZenLum {Cd/m2},WindDir {deg},WindSpd {m/s},TotSkyCvr {.1},OpaqSkyCvr {.1},Visibility {km},Ceiling Hgt {m},PresWeathObs,PresWeathCodes,Precip Wtr {mm},Aerosol Opt Depth {.001},SnowDepth {cm},Days Last Snow,Albedo {.01},Rain {mm},Rain Quantity {hr}
        day = map(lambda r: r[2], self.climate_data[HI:HF+1])
        hr = map(lambda r: r[3], self.climate_data[HI:HF+1])
        self.staTemp = map(lambda r: r[6], self.climate_data[HI:HF+1])
        self.staRhum = map(lambda r: r[8], self.climate_data[HI:HF+1])


        #for i in xrange(len(day)):
        #    print 'date', day[i], hr[i]
        #    print tdb[i]


        """
        self.staTemp = csvread(climate_data,HI,6,[HI,6,HF,6])
        self.staRhum = csvread(climate_data,HI,8,[HI,8,HF,8])
        self.staPres = csvread(climate_data,HI,9,[HI,9,HF,9])
        self.staInfra = csvread(climate_data,HI,12,[HI,12,HF,12])
        self.staHor = csvread(climate_data,HI,13,[HI,13,HF,13])
        self.staDir = csvread(climate_data,HI,14,[HI,14,HF,14])
        self.staDif = csvread(climate_data,HI,15,[HI,15,HF,15])
        self.staUdir = csvread(climate_data,HI,20,[HI,20,HF,20])
        self.staUmod = csvread(climate_data,HI,21,[HI,21,HF,21])
        self.staRobs = csvread(climate_data,HI,33,[HI,33,HF,33])
        self.staHum = zeros(size(self.staTemp,1),1)
        for i=1:size(self.staTemp,1)
          self.staHum(i) = HumFromRHumTemp(self.staRhum(i),...
              self.staTemp(i),self.staPres(i))
        end
        self.staTemp = self.staTemp + 273.15
        """
    def __repr__(self):
        return "Weather @ {a}".format(a=self.location)
"""
function W = HumFromRHumTemp(RH,T,P)

    % Saturation vapour pressure from ASHRAE
    C8=-5.8002206E3 C9=1.3914993 C10=-4.8640239E-2 C11=4.1764768E-5
    C12=-1.4452093E-8 C13=6.5459673
    T=T+273.15
    PWS = exp(C8/T+C9+C10*T+C11*T^2+C12*T^3+C13*log(T))

    PW=RH*PWS/100          % Vapour pressure
    W = 0.62198*PW/(P-PW)  % 4. Specific humidity

 """
