import os
import math
import pprint
ppr = lambda x: pprint.pprint(x)

class RSMDef(object):
    """
    % Rural Site & Vertical Diffusion Model (VDM)
    % Calculates the vertical profiles of air temperature above the weather
    % station per 'The UWG' (2012) Eq. 4, 5, 6.

    properties
        lat;           % latitude (deg)
        lon;           % longitude (deg)
        GMT;           % GMT hour correction
        height         % average obstacle height (m)
        z0r;           % rural roughness length (m)
        disp;          % rural displacement length (m)
        z;             % vertical height (m)
        dz;            % vertical discretization (m)
        nz0;           % layer number at zmt (m)
        nzref;         % layer number at zref (m)
        nzfor;         % layer number at zfor (m)
        nz10;          % layer number at zmu (m)
        nzi;           % layer number at zi_d (m)
        tempProf;      % potential temperature profile at the rural site (K)
        presProf;      % pressure profile at the rural site (Pa)
        tempRealProf;  % real temperature profile at the rural site (K)
        densityProfC;  % density profile at the center of layers (kg m-3)
        densityProfS;  % density profile at the sides of layers (kg m-3)
        windProf;      % wind profile at the rural site (m s-1)
        ublPres;       % Average pressure at UBL (Pa)
    end
    """

    Z_MESO_FILE_NAME = "z_meso.txt"

    def __init__(self,lat,lon,GMT,height,T_init,P_init,parameter,z_meso_path):

        # defines self.z_meso property
        self.load_z_meso(z_meso_path)   # TODO:look up in thesis and define z_meso

        self.lat = lat                  # latitude (deg)
        self.lon = lon                  # longitude (deg)
        self.GMT = GMT                  # GMT hour correction
        self.height = height            # average obstacle height (m)
        self.z0r = 0.1 * height         # rural roughness length (m)
        self.disp = 0.5 * height        # rural displacement lenght (m)

        # vertical grid at the rural site
        self.z  = [0 for x in range(len(self.z_meso)-1)] # Midht btwn each distance interval
        self.dz = [0 for x in range(len(self.z_meso)-1)] # Distance betweeen each interval

        for zi in xrange(len(self.z_meso)-1):
            self.z[zi] = 0.5 * (self.z_meso[zi] + self.z_meso[zi+1])
            self.dz[zi] = self.z_meso[zi+1] - self.z_meso[zi]

        # Define initial booleans
        ll = True
        mm = True
        nn = True
        oo = True
        pp = True

        # Define self.nz0, self.nzref, self.nzfor, self.nz10, self.nzi
        for iz in xrange(len(self.z_meso)-1):
            # self.nz0: self.z index >= reference height for weather station
            eq_th = self.is_near_zero(self.z[iz] - parameter.tempHeight)
            if (eq_th == True or self.z[iz] > parameter.tempHeight) and ll==True:
                self.nz0 = iz   # layer number at zmt (m)
                ll = False

            # self.nzref: self.z index >= reference inversion height
            eq_rh = self.is_near_zero(self.z[iz] - parameter.refHeight)
            if (eq_rh == True or self.z[iz] > parameter.refHeight) and mm==True:
              self.nzref = iz   # layer number at zref (m)
              mm = False

            # self.nzfor: self.z index >= nighttime boundary layer height
            eq_nh = self.is_near_zero(self.z[iz] - parameter.nightBLHeight)
            if (eq_nh == True or self.z[iz] > parameter.nightBLHeight) and nn==True:
              self.nzfor = iz   # layer number at zfor (m)
              nn = False

            # self.nz10: self.z index >= wind height
            eq_wh = self.is_near_zero(self.z[iz] - parameter.windHeight)
            if (eq_wh == True or self.z[iz] > parameter.windHeight) and oo==True:
              self.nz10 = iz    # layer number at zmu (m)
              oo = False

            eq_dh = self.is_near_zero(self.z[iz] - parameter.dayBLHeight)
            if (eq_dh == True or self.z[iz] > parameter.dayBLHeight) and pp==True:
              self.nzi = iz     # layer number at zi_d (m)
              pp = False


        # Define temperature, pressure and density vertical profiles
        self.tempProf = [T_init for x in range(self.nzref+1)]
        self.presProf = [P_init for x in range(self.nzref+1)]

        for iz in xrange(1,self.nzref+1):
            self.presProf[iz] = (self.presProf[iz-1]**(parameter.r/parameter.cp) -\
               parameter.g/parameter.cp * (P_init**(parameter.r/parameter.cp)) * (1./self.tempProf[iz] +\
               1./self.tempProf[iz-1]) * 0.5 * self.dz[iz])**(1./(parameter.r/parameter.cp))

        self.tempRealProf = [T_init for x in range(self.nzref+1)]
        for iz in xrange(self.nzref+1):
           self.tempRealProf[iz] = self.tempProf[iz] * (self.presProf[iz] / P_init)**(parameter.r/parameter.cp)

        self.densityProfC = [None for x in range(self.nzref+1)]
        for iz in xrange(self.nzref+1):
           self.densityProfC[iz] = self.presProf[iz] / parameter.r / self.tempRealProf[iz]

        self.densityProfS = [self.densityProfC[0] for x in range(self.nzref+2)]
        for iz in xrange(1,self.nzref+1):
           self.densityProfS[iz] = (self.densityProfC[iz] * self.dz[iz-1] +\
               self.densityProfC[iz-1] * self.dz[iz]) / (self.dz[iz-1]+self.dz[iz])

        self.densityProfS[self.nzref+1] = self.densityProfC[self.nzref]
        self.windProf = [1 for x in range(self.nzref+1)]

    def __repr__(self):
        return "RSM: obstacle ht={a}".format(
            a=self.height
            )
    def is_near_zero(self,num,eps=1e-16):
        return abs(float(num)) < eps

    def load_z_meso(self,z_meso_path):
        """ Open the z_meso.txt file and return heights as list """

        self.z_meso = []
        z_meso_file_path = os.path.join(z_meso_path, self.Z_MESO_FILE_NAME)

        # Check if exists
        if not os.path.exists(z_meso_file_path):
            raise Exception("z_meso.txt file: '{}' does not exist.".format(uwg_param_file))

        f = open(z_meso_file_path,'r')
        for txtline in f:
            z_ = float("".join(txtline.split())) # Strip all white spaces and change to float
            self.z_meso.append(z_)
        f.close()


    # Ref: The UWG (2012), Eq. (4)
    def VDM(self,forc,rural,parameter,simTime,debug):

        self.tempProf[0] = forc.temp    # Lower boundary condition

        # compute pressure profile
        for iz in reversed(range(self.nzref)[1:]):
           self.presProf[iz-1] = (math.pow(self.presProf[iz],parameter.r/parameter.cp) + \
               parameter.g/parameter.cp*(math.pow(forc.pres,parameter.r/parameter.cp)) * \
               (1./self.tempProf[iz] + 1./self.tempProf[iz-1]) * \
               0.5 * self.dz[iz])**(1./(parameter.r/parameter.cp))

        # compute the real temperature profile
        for iz in xrange(self.nzref):
            self.tempRealProf[iz]= self.tempProf[iz] * \
            (self.presProf[iz]/forc.pres)**(parameter.r/parameter.cp)

        # compute the density profile
        for iz in xrange(self.nzref):
           self.densityProfC[iz] = self.presProf[iz]/parameter.r/self.tempRealProf[iz]

        self.densityProfS[0] = self.densityProfC[0]

        #for iz=2:obj.nzref
        for iz in xrange(1,self.nzref):
           self.densityProfS[iz] = (self.densityProfC[iz] * self.dz[iz-1] + \
               self.densityProfC[iz-1] * self.dz[iz])/(self.dz[iz-1] + self.dz[iz])

        self.densityProfS[self.nzref+1] = self.densityProfC[self.nzref]

        # Ref: The UWG (2012), Eq. (5)
        # compute diffusion coefficient
        cd,ustarRur = self.DiffusionCoefficient(self.densityProfC[0], \
            self.z, self.dz, self.z0r, self.disp, \
            self.tempProf[0], rural.sens, self.nzref, forc.wind, \
            self.tempProf, parameter)

        """
        % solve diffusion equation
        self.tempProf = DiffusionEquation(self.nzref,simTime.dt,...
            self.tempProf,self.densityProfC,self.densityProfS,cd,self.dz);
        % compute wind profile
        for iz=1:self.nzref
            self.windProf(iz) = ustarRur/parameter.vk*...
                log((self.z(iz)-self.disp)/self.z0r);
        end
        % Average pressure
        self.ublPres = 0;
        for iz=1:self.nzfor
            self.ublPres = self.ublPres +...
                self.presProf(iz)*self.dz(iz)/...
                (self.z(self.nzref)+self.dz(self.nzref)/2);
        """

    def DiffusionEquation(self,nz,dt,co,da,daz,cd,dz):
        # Reference?
        """
        cddz = zeros(nz+2,1);
        a = zeros(nz,3);
        c = zeros(nz,1);

        %--------------------------------------------------------------------------
        cddz(1)= daz(1)*cd(1)/dz(1);
        for iz=2:nz
           cddz(iz) = 2.*daz(iz)*cd(iz)/(dz(iz)+dz(iz-1));
        end
        cddz(nz+1) = daz(nz+1)*cd(nz+1)/dz(nz);
        %--------------------------------------------------------------------------
        a(1,1)=0.;
        a(1,2)=1.;
        a(1,3)=0.;
        c(1)=co(1);
        for iz=2:nz-1
           dzv=dz(iz);
           a(iz,1)=-cddz(iz)*dt/dzv/da(iz);
           a(iz,2)=1+dt*(cddz(iz)+cddz(iz+1))/dzv/da(iz);
           a(iz,3)=-cddz(iz+1)*dt/dzv/da(iz);
           c(iz)  =co(iz);
        end
        a(nz,1)=-1.;
        a(nz,2)=1.;
        a(nz,3)=0.;
        c(nz)  =0;
        %--------------------------------------------------------------------------
        co = Invert (nz,a,c);
        """
        return co



    def DiffusionCoefficient(self,rho,z,dz,z0,disp,tempRur,heatRur,nz,uref,th,parameter):

        # Initialization
        Kt = [0 for x in xrange(nz+1)]
        ws = [0 for x in xrange(nz)]
        te = [0 for x in xrange(nz)]
        # Friction velocity (Louis 1979)
        ustar = parameter.vk * uref/math.log((10.-disp)/z0)


        # Monin-Obukhov length
        lengthRur = max(-rho*parameter.cp*ustar**3*tempRur/parameter.vk/parameter.g/heatRur,-50.)

        # Unstable conditions
        if heatRur > 1e-2:
            # Convective velocity scale
            wstar = (parameter.g*heatRur*parameter.dayBLHeight/rho/parameter.cp/tempRur)**(1/3.)
            # Wind profile function
            phi_m = (1-8.*0.1*parameter.dayBLHeight/lengthRur)**(-1./3.)

            for iz in xrange(nz):
                # Mixed-layer velocity scale
                ws[iz] = (ustar**3 + phi_m*parameter.vk*wstar**3*z[iz]/parameter.dayBLHeight)**(1/3.)
                # TKE approximation
                te[iz] = max(ws[iz]**2., 0.01)

        else: # Stable and neutral conditions
            for iz in xrange(nz):
                # TKE approximation
                te[iz] = max(ustar**2.,0.01)

        # lenght scales (l_up, l_down, l_k, l_eps)

        dlu,dld = self.DissipationBougeault(parameter.g,nz,z,dz,te,th)
        print dlu, dld
        """
        [dld,dls,dlk]= LengthBougeault(nz,dld,dlu,z);
        % Boundary-layer diffusion coefficient
        for iz=1:nz
           Kt(iz) = 0.4*dlk(iz)*sqrt(te(iz));
        end
        Kt(nz+1) = Kt(nz);
        """
        return Kt, ustar



    def DissipationBougeault(self,g,nz,z,dz,te,pt):

        dlu = [0 for x in xrange(nz)]
        dld = [0 for x in xrange(nz)]

        for iz in xrange(nz):
            zup=0.
            dlu[iz] = z[nz+1] - z[iz] - dz[iz]/2.
            zzz=0.
            zup_inf=0.
            beta=g/pt[iz]

            for izz in xrange(nz-1):
                dzt=(dz[izz+1]+dz[izz])/2.
                zup=zup-beta*pt[iz]*dzt
                zup=zup+beta*(pt[izz+1]+pt[izz])*dzt/2.
                zzz=zzz+dzt

                if (te[iz]<zup) and ((te[iz]>zup_inf) or self.is_near_zero(te[iz]-zup_inf)):
                    bbb=(pt[izz+1]-pt[izz])/dzt
                    if not self.is_near_zero(bbb-0.):
                        tl=(-beta*(pt[izz]-pt[iz])+ \
                        math.sqrt( max(0.,(beta*(pt[izz]-pt[iz]))**2.+ \
                        2.*bbb*beta*(te[iz]-zup_inf))))/bbb/beta
                    else:
                        tl=(te[iz]-zup_inf)/(beta*(pt[izz]-pt[iz]))
                    dlu[iz]=max(1.,zzz-dzt+tl)
                zup_inf=zup
            zdo=0.
            zdo_sup=0.
            dld[iz]=z[iz]+dz[iz]/2.
            zzz=0.
            for izz in reversed(range(iz)[1:]):
                dzt=(dz[izz-1]+dz[izz])/2.
                zdo=zdo+beta*pt[iz]*dzt
                zdo=zdo-beta*(pt[izz-1]+pt[izz])*dzt/2.
                zzz=zzz+dzt
                if (te[iz]<zdo) and ((te[iz]>zdo_sup) or self.is_near_zero(te[iz]-zdo_sup)):
                    bbb=(pt[izz]-pt[izz-1])/dzt
                    if not self.is_near_zero(bbb-0.):
                        tl=(beta*(pt[izz]-pt[iz])+ \
                            math.sqrt( max(0.,(beta*(pt[izz]-pt[iz]))**2.+ \
                            2.*bbb*beta*(te[iz]-zdo_sup))))/bbb/beta
                    else:
                        tl=(te[iz]-zdo_sup)/(beta*(pt[izz]-pt[iz]))
                    dld[iz]=max(1.,zzz-dzt+tl)
                zdo_sup=zdo

        return dlu,dld

    """
    function [dld,dls,dlk] = LengthBougeault(nz,dld,dlu,z)

        dlg = zeros(nz,1);
        dls = zeros(nz,1);
        dlk = zeros(nz,1);
        for iz=1:nz
            dlg(iz)=(z(iz)+z(iz+1))/2.;
        end

        for iz=1:nz
            dld(iz)=min(dld(iz),dlg(iz));
            dls(iz)=sqrt(dlu(iz)*dld(iz));
            dlk(iz)=min(dlu(iz),dld(iz));
        end

    end

    """

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

        X = [0 for i in xrange(nz)]

        for i in reversed(xrange(nz-1)):
            C[i] = C[i] - A[i][2] * C[i+1]/A[i+1][1]
            A[i][1] = A[i][1] - A[i][2] * A[i+1][0]/A[i+1][1]

        for i in  xrange(1,nz,1):
            C[i] = C[i] - A[i][0] * C[i-1]/A[i-1][1]

        for i in xrange(nz):
            X[i] = C[i]/A[i][1]

        return X
