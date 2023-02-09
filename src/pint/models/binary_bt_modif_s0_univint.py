"""The BT (Blandford & Teukolsky) model."""

from pint.models.parameter import floatParameter
from pint.models.pulsar_binary import PulsarBinary
from pint.models.stand_alone_psr_binaries.BT_model_modif_s0_univint import BTmodel_modif_s0_univint # this has been changed, too
from pint.models.timing_model import MissingParameter

# changed the name of the class from BinaryBT
class BinaryBT_modif_s0_univint(PulsarBinary):
    """Blandford and Teukolsky binary model.

    This binary model is described in Blandford and Teukolshy 1976. It is
    a relatively simple parametrized post-Keplerian model that does not
    support Shapiro delay calculations.

    The actual calculations for this are done in
    :class:`pint.models.stand_alone_psr_binaries.BT_model.BTmodel`.

    Parameters supported:

    .. paramtable::
        :class: pint.models.binary_bt.BinaryBT

    Notes
    -----
    Because PINT's binary models all support specification of multiple orbital
    frequency derivatives FBn, this is capable of behaving like the model called
    BTX in tempo2. The model called BTX in tempo instead supports multiple
    (non-interacting) companions, and that is not supported here. Neither can
    PINT accept "BTX" as an alias for this model.

    See Blandford & Teukolsky 1976, ApJ, 205, 580.
    """

    register = True

    def __init__(self):
        super().__init__()
        self.binary_model_name = "BT_modif_s0_univint"
        self.binary_model_class = BTmodel_modif_s0_univint

        self.add_param(
            floatParameter(
                name="GAMMA",
                value=0.0,
                units="second",
                description="Time dilation & gravitational redshift",
            )
        )

        # this parameter is added
        # will be found by optimalization
        self.add_param(
            floatParameter(
                name="ADM1",
                value=0.0,
                units="",
                description="Scalar DM, univer. int. - amplitude of always present signal",
            )
        )

        # this parameter is also added
        # will be found by optimalization
        self.add_param(
            floatParameter(
                name="ADM2",
                value=0.0,
                units="",
                description="Scalar DM, univer. int. - amplitude of the effective resonant component",
            )
        )

        # mass scale - added
        # will be fixed - particular value from some interval
        self.add_param(
            floatParameter(
                name="MDM",
                value=0.0,
                units="1/day",
                description="DM - mass scale",
            )
        )

        # mass scale - added
        # will be fixed - particular value from some interval
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
        
        # we follow gamma's example and set Adm1 to 0 if the original value is None
        if self.ADM1.value is None:
            self.ADM1.value = "0"
            self.ADM1.frozen = True
        
        # here is the similar situation for the Adm2
        if self.ADM2.value is None:
            self.ADM2.value = "0"
            self.ADM2.frozen = True
        
        # here is the similar situation for the mdm
        if self.MDM.value is None:
            self.MDM.value = "0"
            self.MDM.frozen = True

        # here is the similar situation for the mdm
        if self.BDM.value is None:
            self.BDM.value = "0"
            self.BDM.frozen = True    
