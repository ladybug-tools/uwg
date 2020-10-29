"""Class for Schedule Definition."""

from .utilities import float_positive, REF_BUILTERA_SET, REF_BUILTERA


try:
    str = basestring
except NameError:
    pass


class SchDef(object):
    """Schedule definition class.

    The internal load weekly schedules consists of three lists of 24 values representing hours in
    weekday, Saturday, and Sunday.

    Args:
        elec: Weekly schedule of fractional electricity plug process loads.
        gas: Weekly schedule of fractional gas process loads.
        light: Weekly schedule of fractional light process loads.
        occ: Weekly schedule of fractional occupant number.
        cool: Weekly schedule of cooling temperatures.
        heat: Weekly schedule of heating temperatures.
        q_elec: Maximum electrical plug process load [W/m2].
        q_gas: Maximum gas process load per unit area [W/m2].
        q_light: Maximum light process load per unit area [W/m2].
        n_occ: Maximum number of occupants per unit area [person/m2].
        vent: Maximum ventilation rate per unit area [m3/s/m2].
        swh: Optional property for weekly schedule of fractional hot water rate. This
            property will be a weekly schedule of zero values, as default.
        v_swh: Optional property for maximum volumetric hot water rate per unit area
            [L/hr/m2]. (Default: 0).
        bldtype: Text referring to a building type. By default, 16 building types are
            defined in the UWG according to models from the Department of Energy (DOE).
            Custom building types can also be defined with a new name. Note that this
            value along with the BEMDef builtera must exactly match the identifiers in
            the UWG bld list in order to specify the fraction of total built stock the
            building occupies in the UWG simulation. Choose from the following to
            reference or overwrite a BEM associated with a DOE reference building type:
            'fullservicerestaurant', 'hospital', 'largehotel', 'largeoffice',
            'mediumoffice', 'midriseapartment', 'outpatient', 'primaryschool',
            'quickservicerestaurant', 'secondaryschool', 'smallhotel', 'smalloffice',
            'standaloneretail', 'stripmall', 'supermarket', or 'warehouse'.
        builtera: Text defining building built era. Must be one of the following:
            "pre80" (pre-1980s), "pst80" (post-1980s), or "new" (new construction).
            This value along with the bldtype must exactly match the identifiers in
            the bld array in order to specify the fraction of total built stock the
            building occupies in the UWG simulation.

    Properties:
        * elec
        * gas
        * light
        * occ
        * cool
        * heat
        * swh
        * q_elec
        * q_gas
        * q_light
        * n_occ
        * vent
        * v_swh
        * bldtype
        * builtera
        * zonetype
    """
    DEFAULT_SWH = [[0 for j in range(24)] for i in range(3)]

    def __init__(self, elec, gas, light, occ, cool, heat, q_elec, q_gas, q_light,
                 n_occ, vent, bldtype, builtera, swh=DEFAULT_SWH, v_swh=0):
        self.elec = elec
        self.gas = gas
        self.light = light
        self.occ = occ
        self.cool = cool
        self.heat = heat
        self.q_elec = q_elec
        self.q_gas = q_gas
        self.q_light = q_light
        self.n_occ = n_occ
        self.vent = vent
        self.v_swh = v_swh
        self.swh = swh
        # Properties to be set in readDOE
        self.bldtype = bldtype  # DOE reference building type
        self.builtera = builtera  # pre80, pst80, new
        self.zonetype = None  # climate zone number (only used in testing).

    @property
    def elec(self):
        """Get or set weekly schedule of fractional electricity plug process loads.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._elec

    @elec.setter
    def elec(self, value):
        self._elec = SchDef.check_week_validity(value, 'elec')

    @property
    def gas(self):
        """Get or set weekly schedule of fractional gas process loads.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._gas

    @gas.setter
    def gas(self, value):
        self._gas = SchDef.check_week_validity(value, 'gas')

    @property
    def light(self):
        """Get or set weekly schedule of fractional light process loads.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._light

    @light.setter
    def light(self, value):
        self._light = SchDef.check_week_validity(value, 'light')

    @property
    def occ(self):
        """Get or set weekly schedule of occupant number.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._occ

    @occ.setter
    def occ(self, value):
        self._occ = SchDef.check_week_validity(value, 'occ')

    @property
    def cool(self):
        """Get or set weekly schedule of cooling temperatures.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._cool

    @cool.setter
    def cool(self, value):
        self._cool = SchDef.check_week_validity(value, 'cool')

    @property
    def heat(self):
        """Get or set weekly schedule of heating temperatures.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._heat

    @heat.setter
    def heat(self, value):
        self._heat = SchDef.check_week_validity(value, 'heat')

    @property
    def swh(self):
        """Get or set weekly schedule of fractional hot water rate.

        Weekly schedule consists of three lists of 24 values representing hours in
        weekday, Saturday, and Sunday.
        """
        return self._swh

    @swh.setter
    def swh(self, value):
        self._swh = SchDef.check_week_validity(value, 'swh')

    @property
    def q_elec(self):
        """Get or set maximum electrical plug process load [W/m2]."""
        return self._q_elec

    @q_elec.setter
    def q_elec(self, value):
        self._q_elec = float_positive(value, 'q_elec')

    @property
    def q_gas(self):
        """Get or set maximum gas process load per unit area [W/m2]."""
        return self._q_gas

    @q_gas.setter
    def q_gas(self, value):
        self._q_gas = float_positive(value, 'q_gas')

    @property
    def q_light(self):
        """Get or set maximum light process load per unit area [W/m2]."""
        return self._q_light

    @q_light.setter
    def q_light(self, value):
        self._q_light = float_positive(value, 'q_light')

    @property
    def n_occ(self):
        """Get or set maximum number of occupants per unit area [person/m2]."""
        return self._n_occ

    @n_occ.setter
    def n_occ(self, value):
        self._n_occ = float_positive(value, 'n_occ')

    @property
    def vent(self):
        """Get or set maximum ventilation rate per unit area [m3/s/m2]."""
        return self._vent

    @vent.setter
    def vent(self, value):
        self._vent = float_positive(value, 'vent')

    @property
    def v_swh(self):
        """Get or set maximum volumetric hot water rate per unit area [L/hr/m2]."""
        return self._v_swh

    @v_swh.setter
    def v_swh(self, value):
        self._v_swh = float_positive(value, 'v_swh')

    @property
    def bldtype(self):
        """Get or set text for building type.

        By default, 16 building types are defined in the UWG according to models from
        the Department of Energy (DOE). Choose from the following to reference or
        overwrite a schedule associated with a DOE reference building type:

        * 'fullservicerestaurant'
        * 'hospital'
        * 'largehotel'
        * 'largeoffice'
        * 'medoffice'
        * 'midriseapartment'
        * 'outpatient'
        * 'primaryschool'
        * 'quickservicerestaurant'
        * 'secondaryschool'
        * 'smallhotel'
        * 'smalloffice'
        * 'standaloneretail'
        * 'stripmall'
        * 'supermarket'
        * 'warehouse'

        Custom building types can also be defined with a new name. If a custom SchDef is
        defined with the same name as a reference DOE building type from the list above,
        the reference SchDef will be overwritten by the custom SchDef. Note that this
        value along with the SchDef builtera must exactly match the identifiers in the
        UWG bld list in order to specify the fraction of total built stock the building
        occupies in the UWG simulation.
        """
        return self._bldtype

    @bldtype.setter
    def bldtype(self, value):
        assert isinstance(value, str), 'The bldtype must be a string. ' \
            'Got: {}.'.format(value.lower())
        self._bldtype = value

    @property
    def builtera(self):
        """Get or set text for built era.

        Must be one of the following:

        * 'pre80' - pre-1980s
        * 'pst80' - post-1980s
        * 'new' - new construction

        This value along with the bldtype must exactly match the identifiers in
        the bld array in order to specify the fraction of total built stock the
        building occupies in the UWG simulation.
        """
        return self._builtera

    @builtera.setter
    def builtera(self, value):
        assert isinstance(value, str) and value in REF_BUILTERA_SET, \
            'The builtera must be one of {}.Got: {}.'.format(
                REF_BUILTERA, value.lower())
        self._builtera = value

    @classmethod
    def from_dict(cls, data):
        """Create a SchDef object from a dictionary.

        Args:
            data: A SchDef dictionary following the format below.

        .. code-block:: python

            _example_week = [
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.7, 0.9, 0.9, 0.6, 0.6, 0.6, 0.6, 0.6,
                 0.7, 0.8, 0.9, 0.9, 0.8, 0.8, 0.7, 0.3, 0.2, 0.2],  # Weekday
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                 0.6, 0.7, 0.7, 0.7, 0.7, 0.5, 0.4, 0.3, 0.2, 0.2],  # Saturday
                [0.2, 0.2, 0.2, 0.2, 0.2, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
                 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.2, 0.2]]  # Sunday

            {
            "elec": _example_week,
            "gas": _example_week,
            "light": _example_week,
            "occ": _example_week,
            "cool": _example_week,
            "heat": _example_week,
            "swh": _example_week,
            "q_elec" = q_elec,
            "q_gas" = q_gas,
            "q_light" = q_light,
            "n_occ" = n_occ,
            "vent" = vent,
            "v_swh" = v_swh,
            "bldtype": 0,
            "builtera": 1
            }
        """
        assert data['type'] == 'SchDef', 'Expected ' \
            'SchDef dictionary. Got {}.'.format(data['type'])

        return cls(elec=data['elec'], gas=data['gas'], light=data['light'],
                   occ=data['occ'], cool=data['cool'], heat=data['heat'],
                   swh=data['swh'], q_elec=data['q_elec'], q_gas=data['q_gas'],
                   q_light=data['q_light'], n_occ=data['n_occ'], vent=data['vent'],
                   v_swh=data['v_swh'], bldtype=data['bldtype'],
                   builtera=data['builtera'])

    def to_dict(self):
        """SchDef dictionary representation."""
        base = {'type': 'SchDef'}
        base['elec'] = self.elec
        base['gas'] = self.gas
        base['light'] = self.light
        base['occ'] = self.occ
        base['cool'] = self.cool
        base['heat'] = self.heat
        base['swh'] = self.swh
        base['q_elec'] = self.q_elec
        base['q_gas'] = self.q_gas
        base['q_light'] = self.q_light
        base['n_occ'] = self.n_occ
        base['vent'] = self.vent
        base['v_swh'] = self.v_swh
        base['bldtype'] = self.bldtype
        base['builtera'] = self.builtera
        return base

    @staticmethod
    def check_week_validity(week, name):
        assert isinstance(week, (list, tuple)), 'The {} property must be a ' \
            'list or tuple. Got {}.'.format(name, week)

        assert len(week) == 3, 'The {} property must be a 3 x 24 matrix. Got ' \
            '{} rows.'.format(name, len(week))

        for i, day in enumerate(week):
            assert len(day) == 24, 'The {} property must be a 3 x 24 ' \
                'matrix. Got {} columns for row {}.'.format(name, len(day), i)
            for val in day:
                assert isinstance(val, (float, int)), 'The {} property ' \
                    'must contain 3 lists of numbers. Got : {}.'.format(
                        name, val)
        return week

    def __repr__(self):
        return 'Schedule, bldtype: {}\n builtera: {}\n q_elec: {}\n q_gas: {}\n ' \
            'q_light: {}\n n_occ: {}\n vent: {}\n v_swh: {}\n elec: {}\n ' \
            'gas: {}\n light: {}\n occ: {}\n cool: {}\n heat: {}\n swh: {}\n'.format(
                self.bldtype, self.builtera, self.q_elec, self.q_gas, self.q_light,
                self.n_occ, self.vent, self.v_swh, self.elec, self.gas, self.light,
                self.occ, self.cool, self.heat, self.swh)
