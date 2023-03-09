import astropy.constants as c
import astropy.units as u
import numpy as np

from .binary_generic_modif_s0_univint import PSR_BINARY_DM

# PK
"""
Context: Interaction between the scalar DM and ordinary matter
Adm1 = amplitude of always present binary period modulation;        we want to find this from fitting
Adm2 = effective description of resonance effect                    we want to find this from fitting, too
Bdm  = local DM phase, in deg units                                 this we consider to be given (with a prior = continuous uniform distribution)
mdm  = DM mass (or rather corresponding frequency)                  given, too

omegab = 2*pi/Pb = binary's angular frequency
E = eccentric anomaly, in deg units
E/(360*u.deg)*2*np.pi = E in rad units
"""

# PK
"""
Function representing effect of DM on a binary system
It has two terms:
1st one - always present signal (propto ADM1)
2nd one - effective description of resonance effects
"""
def Rdm(Adm1, Adm2, Bdm, mdm, omegab, E):
    k = mdm/omegab  # array
    return_array = []


    #for j in range(len(k)):
    
    #    if ( k[j] <= 0.1 ):
    #        array_item =  Adm1 * ( k[j] * E[j] * np.sin(Bdm) + 0.5 * k[j] * k[j] * E[j] * E[j] * np.cos(Bdm) )  + Adm2 * E[j] / ( 360 * u.deg ) * 2 * np.pi
    #        return_array.append(array_item)

    #    else:
    #        array_item = Adm1 * ( np.cos(Bdm)  -  np.cos(E[j] * mdm / omegab + Bdm) ) + Adm2 * E[j] / ( 360 * u.deg ) * 2 * np.pi
    #        return_array.append(array_item)

    #    return_array = np.array(return_array)

    #    return return_array



    # PK: This works!!!:
    #return  Adm1 * ( np.cos(Bdm)  -  np.cos(E * mdm / omegab + Bdm) ) + Adm2 * E / ( 360 * u.deg ) * 2 * np.pi

    # PK: Scenerios with small masses
    return Adm1 * ( mdm / omegab * E / ( 360 * u.deg ) * 2 * np.pi * np.sin(Bdm) + 0.5 * (mdm/omegab)**2*np.cos(Bdm) * (E / ( 360 * u.deg ) * 2 * np.pi)**2 ) + Adm2 * E / ( 360 * u.deg ) * 2 * np.pi


# PK
"""
Derivative of Rdm w.r.t. ADM1 and ADM2
"""
def d_Rdm_d_ADM1(Bdm, mdm, omegab, E):
    #k = mdm/omegab
    #if ( k <= 0.1 ):
    #    return k * E * np.sin(Bdm) + 0.5 * k * k * E * E * np.cos(Bdm) 
    #else:
    #    return np.cos(Bdm)  -  np.cos(E * mdm / omegab + Bdm) 
    
    # return np.cos(Bdm)  -  np.cos(E * mdm / omegab + Bdm)
    return mdm / omegab * E / ( 360 * u.deg ) * 2 * np.pi * np.sin(Bdm) + 0.5 * (mdm/omegab)**2*np.cos(Bdm) * (E / ( 360 * u.deg ) * 2 * np.pi)**2 

def d_Rdm_d_ADM2(E):
    return E / ( 360 * u.deg ) * 2 * np.pi



# PK
"""
BTmodel_modif_s0_univint = class for computing binary time delays
                           a subclass of PSR_binary in the binary_generic_modif_s0_univint module in the same directory
                           to interact with PINT, it needs a pulsar binary wrapper - pint/models/binary_bt.py
"""


class BTmodel_modif_s0_univint(PSR_BINARY_DM):
    """This is a class independent from PINT platform for pulsar BT binary model.
    It is a subclass of PSR_BINARY class defined in file binary_generic.py in
    the same dierectory. This class is desined for PINT platform but can be used
    as an independent module for binary delay calculation.
    To interact with PINT platform, a pulsar_binary wrapper is needed.
    See the source file pint/models/binary_bt.py
    Refence
    ---------
    The 'BT' binary model for the pulse period. Model as in:
    W.M. Smart, (1962), "Spherical Astronomy", p35
    Blandford & Teukolsky (1976), ApJ, 205, 580-591

    Return
    ----------
    A bt binary model class with paramters, delay calculations and derivatives.
    Example
    ----------
    >>> import numpy
    >>> t = numpy.linspace(54200.0,55000.0,800)
    >>> binary_model = BTmodel()
    >>> paramters_dict = {'A0':0.5,'ECC':0.01}
    >>> binary_model.update_input(t, paramters_dict)
    Here the binary has time input and parameters input, the delay can be
    calculated.

    @param PB:          Binary period [days]
    @param ECC:         Eccentricity
    @param A1:          Projected semi-major axis (lt-sec)
    @param A1DOT:       Time-derivative of A1 (lt-sec/sec)
    @param T0:          Time of periastron passage (barycentric MJD)
    @param OM:          Omega (longitude of periastron) [deg]
    @param EDOT:        Time-derivative of ECC [0.0]
    @param PBDOT:       Time-derivative of PB [0.0]
    @param OMDOT:       Time-derivative of OMEGA [0.0]
    """

    # PK
    """
    Check:  param_default_value ... should be a directory
            set_param_values()  ... should be a method
            BTdelay_modif_s0_univint
            d_BTdelay_d_par
    """



    def __init__(self, t=None, input_params=None):
        super().__init__()
        self.binary_name = "BT_modif_s0_univint"
        self.binary_params = list(self.param_default_value.keys())
        self.set_param_values()  # Set parameters to default values.
        self.binary_delay_funcs = [self.BTdelay_modif_s0_univint]  # BTdelay is defined later in this page
        self.d_binarydelay_d_par_funcs = [self.d_BTdelay_d_par] # d_BTdelay_d_par is also defined later in this page
        if t is not None:
            self.t = t
        if input_params is not None:
            self.update_input(param_dict=input_params)


    def delayL1(self):
        """First term of Blandford & Teukolsky (1976), ApJ, 205,
        580-591, eq 2.33/ First left-hand term of W.M. Smart, (1962),
        "Spherical Astronomy", p35 delay equation.

        alpha * (cosE-e)
        alpha = a1*sin(omega)
        Here a1 is in the unit of light second as distance
        """
        # return self.a1() / c.c * np.sin(self.omega()) * (np.cos(self.E()) - self.ecc())


        # PK
        """
        Modifications of delayL1 - it includes a contribution coming from the DM

        (?) It seems to me that binary parameters are determinated in time T0
            i.e. t1 = T0 and E(t_1) = 0 
        """
         # waringing amd vs AMD ... check that!
        rdm = Rdm(self.ADM1, self.ADM2, self.BDM, self.MDM, 2*np.pi/self.pb(), self.E())  # I should check unit conversion.

        x = self.a1() / c.c
        alpha_t1 = x * np.sin(self.omega())
        alpha_t  = alpha_t1 * ( 1 + 2/3*rdm )

        return alpha_t * (np.cos(self.E()) - self.ecc())

    def delayL2(self):
        """Second term of Blandford & Teukolsky (1976), ApJ, 205,
        580-591, eq 2.33/ / Second left-hand term of W.M. Smart, (1962),
        "Spherical Astronomy", p35 delay equation.

        (beta + gamma) * sinE
        beta = (1-e^2)*a1*cos(omega)
        Here a1 is in the unit of light second as distance
        """
        #a1 = self.a1() / c.c
        #return (
        #    a1 * np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2) + self.GAMMA
        #) * np.sin(self.E())


        # what follows is a modification
        rdm = Rdm(self.ADM1, self.ADM2, self.BDM, self.MDM, 2*np.pi/self.pb(), self.E())

        x = self.a1() / c.c
        beta_t1 = x * np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2)
        beta_t = beta_t1 * ( 1 + 2/3*rdm ) 

        return (beta_t + self.GAMMA)*np.sin(self.E())

        

    def delayR(self):
        """Third term of Blandford & Teukolsky (1976), ApJ, 205,
        580-591, eq 2.33 / Right-hand term of W.M. Smart, (1962),
        "Spherical Astronomy", p35 delay equation.

        (alpha*(cosE-e)+(beta+gamma)*sinE)*(1-alpha*sinE - beta*sinE)/
        (pb*(1-e*coeE))
        (alpha*(cosE-e)+(beta+gamma)*sinE) is define in delayL1
        and delayL2
        """
        #a1 = self.a1() / c.c
        #omega = self.omega()
        #ecc = self.ecc()
        #E = self.E()
        #num = a1 * np.cos(omega) * np.sqrt(1 - ecc**2) * np.cos(E) - a1 * np.sin(
        #    omega
        #) * np.sin(E)
        #den = 1.0 - ecc * np.cos(E)

        # In BTmodel.C, they do not use pbprime here, just pb...
        # Is it not more appropriate to include the effects of PBDOT?
        # return 1.0 - 2*np.pi*num / (den * self.pbprime())
        #return 1.0 - 2 * np.pi * num / (den * self.pb().to(u.second))


        # what follows is a modification
        rdm = Rdm(self.ADM1, self.ADM2, self.BDM, self.MDM, 2*np.pi/self.pb(), self.E())

        x = self.a1() / c.c
        
        alpha_t1 = x * np.sin(self.omega())
        alpha_t  = alpha_t1*( 1 + 2/3*rdm )
        beta_t1 = x * np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2)
        beta_t = beta_t1 * ( 1 + 2/3*rdm ) 

        num = alpha_t * np.sin(self.E()) - (beta_t+self.GAMMA)*np.cos(self.E())
        den = ( 1.0 - self.ecc() * np.cos(self.E()) ) * self.pb().to(u.second)

        return 1.0 + 2*np.pi*num/den


    def BTdelay_modif_s0_univint(self):
        """Full (but modified) BT model delay"""
        return (self.delayL1() + self.delayL2()) * self.delayR()

    # NOTE: Below, OMEGA is supposed to be in RADIANS!
    # TODO: Fix UNITS!!!
    def d_delayL1_d_E(self):
        a1 = self.a1() / c.c
        return -a1 * np.sin(self.omega()) * np.sin(self.E())

    def d_delayL2_d_E(self):
        a1 = self.a1() / c.c
        return (
            a1 * np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2) + self.GAMMA
        ) * np.cos(self.E())

    def d_delayL1_d_A1(self):
        return np.sin(self.omega()) * (np.cos(self.E()) - self.ecc()) / c.c

    def d_delayL1_d_A1DOT(self):
        return self.tt0 * self.d_delayL1_d_A1()

    def d_delayL2_d_A1(self):
        return (
            np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2) * np.sin(self.E()) / c.c
        )

    def d_delayL2_d_A1DOT(self):
        return self.tt0 * self.d_delayL2_d_A1()

    def d_delayL1_d_OM(self):
        a1 = self.a1() / c.c
        return a1 * np.cos(self.omega()) * (np.cos(self.E()) - self.ecc())

    def d_delayL1_d_OMDOT(self):
        return self.tt0 * self.d_delayL1_d_OM()

    def d_delayL2_d_OM(self):
        a1 = self.a1() / c.c
        return (
            -a1 * np.sin(self.omega()) * np.sqrt(1 - self.ecc() ** 2) * np.sin(self.E())
        )

    def d_delayL2_d_OMDOT(self):
        return self.tt0 * self.d_delayL2_d_OM()

    def d_delayL1_d_ECC(self):
        a1 = self.a1() / c.c
        return a1 * np.sin(self.omega()) + self.d_delayL1_d_E() * self.d_E_d_ECC()

    def d_delayL1_d_EDOT(self):
        return self.tt0 * self.d_delayL1_d_ECC()

    def d_delayL2_d_ECC(self):
        a1 = self.a1() / c.c
        num = -a1 * np.cos(self.omega()) * self.ecc() * np.sin(self.E())
        den = np.sqrt(1 - self.ecc() ** 2)
        return num / den + self.d_delayL2_d_E() * self.d_E_d_ECC()

    def d_delayL2_d_EDOT(self):
        return self.tt0 * self.d_delayL2_d_ECC()

    def d_delayL1_d_GAMMA(self):
        return np.zeros(len(self.t)) * u.second / u.second

    def d_delayL2_d_GAMMA(self):
        return np.sin(self.E())

    def d_delayL1_d_T0(self):
        return self.d_delayL1_d_E() * self.d_E_d_T0()

    def d_delayL2_d_T0(self):
        return self.d_delayL2_d_E() * self.d_E_d_T0()

    #PK new derivatives of delayL1 with respect to ADM1 and ADM2
    # Note: I ignore terms with ADMX in derivatives, since their are not expected to be dominant
    def d_delayL1_d_ADM1(self):

        omegab = 2*np.pi/(self.pb())

        alpha_0 = self.a1()/c.c * np.sin(self.omega())

        return 2/3 * alpha_0 * (  np.cos(self.E()) - self.ecc()  ) * d_Rdm_d_ADM1(self.BDM, self.MDM, omegab, self.E() )        

    def d_delayL1_d_ADM2(self):

        alpha_0 = self.a1()/c.c * np.sin(self.omega())

        return 2/3 * alpha_0 * (  np.cos(self.E()) - self.ecc()  ) * d_Rdm_d_ADM2(self.E())   


    def d_delayL2_d_ADM1(self):

        omegab = 2*np.pi/(self.pb())

        beta_0 = self.a1()/c.c * np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2)

        return 2/3 * beta_0 *  np.sin(self.E()) * d_Rdm_d_ADM1(self.BDM, self.MDM, omegab, self.E() )        

    def d_delayL2_d_ADM2(self):

        beta_0 = self.a1()/c.c * np.cos(self.omega()) * np.sqrt(1 - self.ecc() ** 2)

        return 2/3 * beta_0 * np.sin(self.E()) * d_Rdm_d_ADM2(self.E())  






    def d_delayL1_d_par(self, par):
        if par not in self.binary_params:
            errorMesg = par + " is not in binary parameter list."
            raise ValueError(errorMesg)

        par_obj = getattr(self, par)
        if hasattr(self, "d_delayL1_d_" + par):
            func = getattr(self, "d_delayL1_d_" + par)
            return func()
        else:
            if par in self.orbits_cls.orbit_params:
                return self.d_delayL1_d_E() * self.d_E_d_par(par)
            else:
                return np.zeros(len(self.t)) * u.second / par_obj.unit

    def d_delayL2_d_par(self, par):
        if par not in self.binary_params:
            errorMesg = par + " is not in binary parameter list."
            raise ValueError(errorMesg)

        par_obj = getattr(self, par)
        if hasattr(self, "d_delayL2_d_" + par):
            func = getattr(self, "d_delayL2_d_" + par)
            return func()
        else:
            if par in self.orbits_cls.orbit_params:
                return self.d_delayL2_d_E() * self.d_E_d_par(par)
            else:
                return np.zeros(len(self.t)) * u.second / par_obj.unit

    def d_BTdelay_d_par(self, par):
        return self.delayR() * (self.d_delayL1_d_par(par) + self.d_delayL2_d_par(par))
