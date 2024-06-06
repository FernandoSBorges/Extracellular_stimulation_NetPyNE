from dataclasses import dataclass
from enum import Enum, auto


class AxonModificationMode(Enum):
    KEEP_AXON = auto()
    MYELINATED_AXON = auto()
    REPLACE_WITH_STRAIGHT_AXON = auto()
    REPLACE_WITH_MYELINATED_STRAIGHT_AXON = auto()


@dataclass
class CellModificationParameters:
    """
    A dataclass with parameters for the NeuronCell class.
    """
    axon_modification_mode: AxonModificationMode = AxonModificationMode.KEEP_AXON
    soma_area_scaling_factor: float = 1
    axon_diameter_scaling_factor: float = 1
    main_axon_diameter_scaling_factor: float = 1
    apic_diameter_scaling_factor: float = 1
    dend_diameter_scaling_factor: float = 1
    dend_length_scaling_factor: float = 1

    min_myelin_diameter: float = 0.2
    min_myelin_length: float = 20
    max_myelin_order: int = 0

@dataclass
class MaxhModificationParameters(CellModificationParameters):
    """
    Parameters to create a neuron of an adult, human with myelinated axon using L2/3 PC rat BB: human Eyal 2018
    """

    axon_modification_mode: AxonModificationMode = AxonModificationMode.MYELINATED_AXON
    soma_area_scaling_factor: float = 2.453
    axon_diameter_scaling_factor: float = 2.453
    main_axon_diameter_scaling_factor: float = 1
    apic_diameter_scaling_factor: float = 1.876
    dend_diameter_scaling_factor: float = 1.946
    dend_length_scaling_factor: float = 1.17
    
@dataclass
class UmaxHModificationParameters(CellModificationParameters):
    """
    Parameters to create a neuron of an adult, human with unmyelinated axon using L2/3 PC rat BB: human Eyal 2018
    """

    axon_modification_mode: AxonModificationMode = AxonModificationMode.KEEP_AXON
    soma_area_scaling_factor: float = 2.453
    axon_diameter_scaling_factor: float = 2.453
    main_axon_diameter_scaling_factor: float = 1
    apic_diameter_scaling_factor: float = 1.876
    dend_diameter_scaling_factor: float = 1.946
    dend_length_scaling_factor: float = 1.17
