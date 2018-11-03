from __future__ import division

try:
    range = xrange
except NameError:
    pass

import math


class Element(object):
    """
    uwg Element

    # Note: In matlab not all instance variables are instantiated. They are assumed to be a 0-by-0 empty matrix
    # https://www.mathworks.com/help/matlab/matlab_oop/specifying-properties.html

    Attributes:
        albedo;          % outer surface albedo
        emissivity;      % outer surface emissivity
        layerThickness;  % vector of layer thicknesses (m)
        layerThermalCond;% vector of layer thermal conductivities (W m-1 K-1)
        layerVolHeat;    % vector of layer volumetric heat (J m-3 K-1)
        vegCoverage;     % surface vegetation coverage
        layerTemp;       % vector of layer temperatures (K)
        waterStorage;    % thickness of water film (m) (only for horizontal surfaces)
        horizontal;      % 1-horizontal, 0-vertical
        solRec;          % solar radiation received (W m-2)
        infra;           % net longwave radiation (W m-2)
        lat;             % surface latent heat flux (W m-2)
        sens;            % surface sensible heat flux (W m-2)
        solAbs;          % solar radiation absorbed (W m-2)
        aeroCond;        % convective heat transfer
        T_ext;           % external surface temperature
        T_int;           % internal surface temperature
        flux;            % external surface heat flux
    """

    THICKNESSLST_EQ_MATERIALLST_MSG = \
    "-----------------------------------------\n" +\
    "ERROR: the number of layer thickness must\n" +\
    "match the number of layer materials\n"
    "-----------------------------------------"
    CONDUCTION_INPUT_MSG = 'ERROR: check input parameters in the Conduction routine'

    def __init__(self, alb, emis, thicknessLst, materialLst, vegCoverage, T_init, horizontal,name=None):
        if len(thicknessLst) != len(materialLst):
            raise Exception(self.THICKNESSLST_EQ_MATERIALLST_MSG)
        else:
            self._name = name                                       # purely for internal process
            self.albedo = alb                                       # outer surface albedo
            self.emissivity = emis                                  # outer surface emissivity
            self.layerThickness = thicknessLst                      # vector of layer thicnesses (m)
            self.layerThermalCond = [0. for i in materialLst]       # vector of layer thermal conductivity (W m-1 K-1)
            self.layerVolHeat = [0. for i in materialLst]           # vector of layer volumetric heat (J m-3 K-1)

            # Create list of layer k and (Cp*density) from materialLst properties
            for i in range(len(materialLst)):
              self.layerThermalCond[i] = materialLst[i].thermalCond
              self.layerVolHeat[i] = materialLst[i].volHeat

            self.vegCoverage = vegCoverage                          # surface vegetation coverage
            self.layerTemp = [T_init] * len(thicknessLst)           # vector of layer temperatures (K)
            self.waterStorage = 0.                                  # thickness of water film (m) for horizontal surfaces only
            self.infra = 0.                                         # net longwave radiation (W m-2)
            self.horizontal = horizontal                            # 1-horizontal, 0-vertical
            self.sens = 0.                                          # surface sensible heat flux (W m-2)

            # B/c we have to explicity define this in python. Set as None
            self.solRec = None                                       # solar radiation received (W m-2)
            self.lat = None                                          # surface latent heat flux (W m-2)
            self.solAbs = None                                       # solar radiation absorbed (W m-2)
            self.aeroCond = None                                     # convective heat transfer
            self.T_ext = None                                        # external surface temperature
            self.T_int = None                                        # internal surface temperature
            self.flux = None                                         # external surface heat flux

    def __repr__(self):
        # Returns some representative wall properties
        s1 = "Element: {a}\n\tlayerNum={b}, totaldepth={c}\n\t".format(
            a=self._name,
            b=len(self.layerThickness),
            c=sum(self.layerThickness)
            )
        s2 = "e={d}, a={e}\n\tr_val={f}, Cp*dens_avg={g}\n\tlayerTemp: {h}".format(
            d=self.emissivity,
            e=self.albedo,
            f=round(sum(self.layerThermalCond)/2.,2),
            g=round(sum(self.layerVolHeat)/2.,2),
            h=self.layerTemp
            )
        return s1 + s2

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def SurfFlux(self,forc,parameter,simTime,humRef,tempRef,windRef,boundCond,intFlux):
        """ Calculate net heat flux, and update element layer temperatures
        """

        # Calculated per unit area (m^2)
        dens = forc.pres/(1000*0.287042*tempRef*(1.+1.607858*humRef)) # air density (kgd m-3)
        self.aeroCond = 5.8 + 3.7 * windRef         # Convection coef (ref: uwg, eq. 12))

        if (self.horizontal):     # For roof, mass, road
            # Evaporation (m s-1), Film water & soil latent heat
            if not self.is_near_zero(self.waterStorage) and self.waterStorage > 0.0:
                # N.B In the current uwg code, latent heat from evapotranspiration, stagnant water,
                # or anthropogenic sources is not modelled due to the difficulty of validation, and
                # lack of reliability of precipitation data from EPW files.Therefore this condition
                # is never run because all elements have had their waterStorage hardcoded to 0.
                qtsat = self.qsat([self.layerTemp[0]],[forc.pres],parameter)[0]
                eg = self.aeroCond*parameter.colburn*dens*(qtsat-humRef)/parameter.waterDens/parameter.cp
                self.waterStorage = min(self.waterStorage + simTime.dt*(forc.prec-eg),parameter.wgmax)
                self.waterStorage = max(self.waterStorage,0.) # (m)
            else:
                eg = 0.
            soilLat = eg*parameter.waterDens*parameter.lv

            # Winter, no veg
            if simTime.month < parameter.vegStart and simTime.month > parameter.vegEnd:
                self.solAbs = (1.-self.albedo)*self.solRec # (W m-2)
                vegLat = 0.
                vegSens = 0.
            else:    # Summer, veg
                self.solAbs = ((1.-self.vegCoverage)*(1.-self.albedo)+self.vegCoverage*(1.-parameter.vegAlbedo))*self.solRec
                vegLat = self.vegCoverage*parameter.grassFLat*(1.-parameter.vegAlbedo)*self.solRec
                vegSens = self.vegCoverage*(1.-parameter.grassFLat)*(1.-parameter.vegAlbedo)*self.solRec

            self.lat = soilLat + vegLat
            # Sensible & net heat flux
            self.sens = vegSens + self.aeroCond*(self.layerTemp[0]-tempRef)
            self.flux = -self.sens + self.solAbs + self.infra - self.lat # (W m-2)

        else:               # For vertical surfaces (wall)
            self.solAbs = (1.-self.albedo)*self.solRec
            self.lat = 0.

            # Sensible & net heat flux
            self.sens = self.aeroCond*(self.layerTemp[0]-tempRef)
            self.flux = -self.sens + self.solAbs + self.infra - self.lat # (W m-2)

        self.layerTemp = self.Conduction(simTime.dt, self.flux, boundCond, forc.deepTemp, intFlux)
        self.T_ext = self.layerTemp[0]
        self.T_int = self.layerTemp[-1]

    def Conduction(self, dt, flx1, bc, temp2, flx2):
        """
        Solve the conductance of heat based on of the element layers.
        arg:
            flx1  : net heat flux on surface
            bc    : boundary condition parameter (1 or 2)
            temp2 : deep soil temperature (ave of air temperature)
            flx2  : surface flux (sum of absorbed, emitted, etc.)

        key prop:
            za = [[ x00, x01, x02 ... x0w ]
                  [ x10, x11, x12 ... x1w ]
                            ...
                  [ xh0, xh1, xh2 ... xhw ]]

            where h = matrix row index    = element layer number
                  w = matrix column index = 3

        """
        t = self.layerTemp          # vector of layer temperatures (K)
        hc = self.layerVolHeat      # vector of layer volumetric heat (J m-3 K-1)
        tc = self.layerThermalCond  # vector of layer thermal conductivities (W m-1 K-1)
        d = self.layerThickness     # vector of layer thicknesses (m)

        # flx1                      : net heat flux on surface
        # bc                        : boundary condition parameter (1 or 2)
        # temp2                     : deep soil temperature (avg of air temperature)
        # flx2                      : surface flux (sum of absorbed, emitted, etc.)

        fimp = 0.5                  # implicit coefficient
        fexp = 0.5                  # explicit coefficient
        num = len(t)                # number of layers

        # Mean thermal conductivity over distance between 2 layers (W/mK)
        tcp = [0 for x in range(num)]
        # Thermal capacity times layer depth (J/m2K)
        hcp = [0 for x in range(num)]
        # lower, main, and upper diagonals
        za = [[0 for y in range(3)] for x in range(num)]
        # RHS
        zy = [0 for x in range(num)]

        #--------------------------------------------------------------------------
        # Define the column vectors for heat capactiy and conductivity
        hcp[0] = hc[0] * d[0]
        for j in range(1,num):
            tcp[j] = 2. / (d[j-1] / tc[j-1] + d[j] / tc[j])
            hcp[j] = hc[j] * d[j]

        #--------------------------------------------------------------------------
        # Define the first row of za matrix, and RHS column vector
        za[0][0] = 0.
        za[0][1] = hcp[0]/dt + fimp*tcp[1]
        za[0][2] = -fimp*tcp[1]
        zy[0] = hcp[0]/dt*t[0] - fexp*tcp[1]*(t[0]-t[1]) + flx1

        #--------------------------------------------------------------------------
        # Define other rows
        for j in range(1,num-1):
          za[j][0] = fimp*(-tcp[j])
          za[j][1] = hcp[j]/dt + fimp*(tcp[j]+tcp[j+1])
          za[j][2] = fimp*(-tcp[j+1])
          zy[j] = hcp[j]/dt * t[j] + fexp * \
            (tcp[j]*t[j-1] - tcp[j]*t[j] - tcp[j+1]*t[j] + tcp[j+1]*t[j+1])

        #--------------------------------------------------------------------------
        # Boundary conditions
        if self.is_near_zero(bc-1.): # heat flux
            za[num-1][0] = fimp * (-tcp[num-1])
            za[num-1][1] = hcp[num-1]/dt + fimp*tcp[num-1]
            za[num-1][2] = 0.
            zy[num-1] = hcp[num-1]/dt*t[num-1] + fexp*tcp[num-1]*(t[num-2]-t[num-1]) + flx2
        elif self.is_near_zero(bc-2.): # deep-temperature
            za[num-1][0] = 0.
            za[num-1][1] = 1.
            za[num-1][2] = 0.
            zy[num-1] = temp2
        else:
            raise Exception(self.CONDUCTION_INPUT_MSG)

        #--------------------------------------------------------------------------

        zx = self.invert(num,za,zy)
        #t(:) = zx(:);
        return zx # return zx as 1d vector of templayers

    def qsat(self,temp,pres,parameter):
        """
        Calculate (qsat_lst) vector of saturation humidity from:
            temp = vector of element layer temperatures
            pres = pressure (at current timestep).
        """
        gamw = (parameter.cl - parameter.cpv) / parameter.rv
        betaw = (parameter.lvtt/parameter.rv) + (gamw * parameter.tt)
        alpw = math.log(parameter.estt) + (betaw /parameter.tt) + (gamw * math.log(parameter.tt))
        work2 = parameter.r/parameter.rv
        foes_lst = [0 for i in range(len(temp))]
        work1_lst = [0 for i in range(len(temp))]
        qsat_lst = [0 for i in range(len(temp))]

        for i in range(len(temp)):
          # saturation vapor pressure
          foes_lst[i] = math.exp( alpw - betaw/temp[i] - gamw*math.log(temp[i])  )
          work1_lst[i] = foes_lst[i]/pres[i]
          # saturation humidity
          qsat_lst[i] = work2*work1_lst[i] / (1. + (work2-1.) * work1_lst[i])

        return qsat_lst


    def invert(self,nz,A,C):
        """
        Inversion and resolution of a tridiagonal matrix
                 A X = C
        Input:
         nz number of layers
         a(*,1) lower diagonal (Ai,i-1)
         a(*,2) principal diagonal (Ai,i)
         a(*,3) upper diagonal (Ai,i+1)
         c
        Output
         x     results
        """

        X = [0 for i in range(nz)]

        for i in reversed(range(nz-1)):
            C[i] = C[i] - A[i][2] * C[i+1]/A[i+1][1]
            A[i][1] = A[i][1] - A[i][2] * A[i+1][0]/A[i+1][1]

        for i in  range(1,nz,1):
            C[i] = C[i] - A[i][0] * C[i-1]/A[i-1][1]

        for i in range(nz):
            X[i] = C[i]/A[i][1]

        return X
