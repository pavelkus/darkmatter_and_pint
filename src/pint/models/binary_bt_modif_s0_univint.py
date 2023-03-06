""" The BT (Blandford & Teukolsky) model. """
""" PK: Modifications to include the spin-0 DM! """

from pint.models.parameter import floatParameter
from pint.models.pulsar_binary import PulsarBinary  # DO I NEED TO MODIFY STH THERE???
from pint.models.stand_alone_psr_binaries.BT_model_modif_s0_univint import BTmodel_modif_s0_univint # this has been changed, too; notice: "BT_model_modif_s0_univint" vs "BTmodel_modif_s0_univint" ... the difference lies in the placeholder
from pint.models.timing_model import MissingParameter

# changed the name of the class from BinaryBT
class BinaryBT_modif_s0_univint(PulsarBinary):
    """The binary model of Blandford and Teukolsky, which is extended to include the effect of DM on the pulsar period.

    The actual calculations for this are done in
    :class:`pint.models.stand_alone_psr_binaries.BT_model_modif_s0_univint.BTmodel_modif_s0_univint`.

    See Blandford & Teukolsky 1976, ApJ, 205, 580.
    """

    register = True

    def __init__(self):
        super().__init__()
        self.binary_model_name = "BT_modif_s0_univint"
        self.binary_model_class = BTmodel_modif_s0_univint

        # "GAMMA is added, but with zero value and is not fitted - it is not effectively present anywhere"
        self.add_param(
            floatParameter(
                name="GAMMA",   
                value=0.0,
                units="second",
                description="Time dilation & gravitational redshift",
            )
        )

        # 1st DM parameter is added
        # will be searched by optimizing
        self.add_param(
            floatParameter(
                name="ADM1",
                value=0.0,
                units="",
                description="Scalar DM, univer. int. - amplitude of always present signal",
            )
        )

        # 2nd DM parameter is added
        # will be searched by optimizing
        self.add_param(
            floatParameter(
                name="ADM2",
                value=0.0,
                units="",
                description="Scalar DM, univer. int. - amplitude of the effective resonant component",
            )
        )

        # 3rd DM parameter is added: DM angular frequency (closely related to the DM scale) 
        # will be "given" - a particular value from an interval
        # Note: DM mass scale: [1e-22, 1e-18] eV
        self.add_param(
            floatParameter(
                name="MDM",
                value=0.0,
                units="1/day",
                description="DM - angular oscillation frequency",
            )
        )

        # 4th DM parameter is added: DM phase
        # will be "given" -  a particular value from the interval [0,2pi]
        self.add_param(
            floatParameter(
                name="BDM",
                value=0.0,
                units="deg",
                description="DM - phase",
            )
        )

        self.remove_param("M2")
        self.remove_param("SINI")

    def validate(self):
        """Validate BT model parameters"""
        super().validate()
        for p in ("T0", "A1"):
            if getattr(self, p).value is None:
                raise MissingParameter("BT_modif_s0_univint", p, "%s is required for BT_modif_s0_univint" % p) # BT_modif_s0_univint instead of BT

        # If any *DOT is set, we need T0
        for p in ("PBDOT", "OMDOT", "EDOT", "A1DOT"):
            if getattr(self, p).value is None:
                getattr(self, p).value = "0"
                getattr(self, p).frozen = True

        if self.GAMMA.value is None:
            self.GAMMA.value = "0"
            self.GAMMA.frozen = True
        
        # we follow gamma's example and set ADM1 to 0 if the original value is None
        if self.ADM1.value is None:
            self.ADM1.value = "0"
            self.ADM1.frozen = True
        
        # here is the similar situation for the ADM2
        if self.ADM2.value is None:
            self.ADM2.value = "0"
            self.ADM2.frozen = True
        
        # here is the similar situation for the MDM
        if self.MDM.value is None:
            self.MDM.value = "0"
            self.MDM.frozen = True

        # here is the similar situation for the BDM
        if self.BDM.value is None:
            self.BDM.value = "0"
            self.BDM.frozen = True    
