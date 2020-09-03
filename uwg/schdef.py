"""Class for Schedule Definition."""


class SchDef(object):
    """Schedule class

    Args:
        Elec: Matrix of numbers for schedule for electricity (WD,Sat,Sun)
        Gas: Matrix of numbers for schedule for gas (WD,Sat,Sun)
        Light: Matrix of numbers for schedule for light (WD,Sat,Sun)
        Occ: Matrix of numbers for chedule for occupant (WD,Sat,Sun)
        Cool: Matrix of numbers for temperature schedule for cooling (WD,Sat,Sun)
        Heat: Matrix of numbers for temperature schedule for heating (WD,Sat,Sun)
        SWH: Matrix of numbers for hot water schedule

    Properties:
        * Elec
        * Gas
        * Light
        * Occ
        * Cool
        * Heat
        * SWH
    """

    def __init__(self, Elec=None, Gas=None, Light=None, Occ=None, Cool=None, Heat=None,
                 SWH=None):
        self.Elec = Elec
        self.Gas = Gas
        self.Light = Light
        self.Occ = Occ
        self.Cool = Cool
        self.Heat = Heat
        self.SWH = SWH

    def __repr__(self):
        return 'Schedule:(weekday from 8:00 - 18:00) \n Heating: {}\n ' \
            'Cooling: {}'.format(
                'Null' if self.Heat is None else self.Heat[0][7:17],
                'Null' if self.Cool is None else self.Cool[0][7:17])
