class SchDef(object):
    """
    Schedule class\

    Elec        % Schedule for electricity (WD,Sat,Sun)
    Gas         % Schedule for gas (WD,Sat,Sun)
    Light       % Schedule for light (WD,Sat,Sun)
    Occ         % Schedule for occupant (WD,Sat,Sun)
    Cool        % Temperature schedule for cooling (WD,Sat,Sun)
    Heat        % Temperature schedule for heating (WD,Sat,Sun)
    SWH         % Hot water schedule (tbd)

    % Internal Heat Load from DOE
    Qelec       % W/m^2 (max) for electrical plug process
    Qgas        % W/m^2 (max) for gas process
    Qlight      % W/m^2 (max) for light process
    Nocc        % #/m^2 (max) for occupancy
    Vswh        % Hot water vol/hr (max)
    Vent        % litres/s/person for ventilation
    """

    def __init__(self,Elec=None,Gas=None,Light=None,Occ=None,Cool=None,Heat=None,SWH=None):
        self.Elec = Elec
        self.Gas = Gas
        self.Light = Light
        self.Occ = Occ
        self.Cool = Cool
        self.Heat = Heat
        self.SWH = SWH

    def __repr__(self):
        return "Schedule:(weekday from 8:00 - 18:00) \n Heating: {a}\n Cooling: {b}".format(
            a= "Null" if self.Heat is None else self.Heat[0][7:17],
            b= "Null" if self.Cool is None else self.Cool[0][7:17]
            )
