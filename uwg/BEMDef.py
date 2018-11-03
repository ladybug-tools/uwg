from __future__ import print_function


class BEMDef(object):
    """
    Building Energy Model (BEM) class definition


    attributes:
        building;       # building class type
        mass;           # mass element class
        wall;           # wall element class
        roof;           # roof element class

        frac;           # fraction of the urban floor space of this typology
        fl_area;        # Total building floor area in the urban area

        ElecTotal;      # Actual total electricity (W/m^2)
        Elec;           # Actual electricity consumption(W/m^2)
        Gas;            # Actual gas consumption(W/m^2)
        Light;          # Actual light (W/m^2)
        Qocc;           # Actual heat load from occupant (W/m^2)
        SWH;            # Actual hot water usage
        Nocc;           # Actual number of people

        T_wallex;       # Wall surface temp (ext)
        T_wallin;       # Wall surface temp (int)
        T_roofex;       # Roof surface temp
        T_roofin;       #
        Q_wallex;       # Wall surface flux (ext)
        Q_wallin;       # Wall surface flux (int)
        Q_roofex;
        Q_roofin;
    """

    def __init__(self,building,mass,wall,roof,frac):
        self.building = building    # building class type
        self.mass = mass            # mass element class
        self.wall = wall            # wall element class
        self.roof = roof            # roof element class
        self.frac = frac            # fraction of the urban floor space of this typology (~ bld_footprint/urba_footprint ?)

    def __repr__(self):
        return "BEMDef: Type = {}, Zone = {}, Era = {}, Construction = {}".format(
            self.building.Type,
            self.building.Zone,
            self.building.Era,
            self.wall._name
            )
