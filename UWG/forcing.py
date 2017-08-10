class Forcing (self):
    """
    FORCING
    Force Restore method to estimate ground soil temperature?? -sv

    properties
        deepTemp;   % deep soil temperature (K)
        infra;      %
        wind;
        uDir;
        hum;
        pres;       % Pressure (Pa)
        temp;
        rHum;
        prec;       % Precipitation (m s-1)
        dir;
        dif;
        waterTemp;  % ground water temp, set to temp at 2m

    """

    def __init__(self,staTemp,weather):
        self.deepTemp = mean(staTemp)
        self.waterTemp = mean(staTemp)
        self.infra = [weather.staInfra]
        self.uDir = [weather.staUdir]
        self.hum  = [weather.staHum]
        self.pres = [weather.staPres]
        self.temp = [weather.staTemp]
        self.rHum = [weather.staRhum]
        self.dir = [weather.staDir]
        self.dif = [weather.staDif]
        self.prec = [weather.staRobs]/3.6e6
        self.wind = [weather.staUmod]
