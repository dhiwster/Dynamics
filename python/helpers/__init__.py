from .constants import (
    ELEMENTARY_CHARGE_C,
    HBAR_J_S,
    MICROELECTRONVOLT_TO_JOULE,
    MICROSECOND_TO_SECOND,
    NANOSECOND_TO_SECOND,
    RAD_PER_NS_PER_MICROELECTRONVOLT,
    RAD_PER_US_PER_MICROELECTRONVOLT,
    TWO_PI,
)
from .hamiltonian import (
    SingleDQDSpectrum,
    single_dqd_eigensystem,
    single_dqd_hamiltonian,
    single_dqd_numpy_hamiltonian,
    single_dqd_qubit_splitting,
    single_dqd_spectrum,
    single_dqd_tau_z_matrix_element,
)
from .operators import (
    SingleDQDOperatorsNumpy,
    single_dqd_numpy_operators,
)
from .pulses import Pulse, constant, piecewise, square

try:
    from .gates import (
        build_H,
        drift_H,
        edsr1_op,
        edsr2_op,
        eps1_op,
        eps2_op,
        initial_full_state,
        iswap_H,
        vacuum_eigenstates,
        xrot_H,
        zrot_H,
        zrot_delta_freq_GHz,
        zrot_time_ns,
    )
    from .hamiltonian import single_dqd_numpy_hamiltonian, single_dqd_qutip_hamiltonian
    from .operators import SingleDQDOperatorsQutip, single_dqd_qutip_operators
    from .sequence import DQDSequenceCompiler
    from .solver import (
        dqd_populations,
        dqd_qubit_unitary,
        run_state,
        run_unitary,
        two_qubit_unitary,
        vacuum_subspace_unitary,
    )
    from .system import DQDsystem
    from .viz import plot_pulse_schedule, print_summary
except ModuleNotFoundError as exc:
    if exc.name != "qutip":
        raise

__all__ = [
    "ELEMENTARY_CHARGE_C",
    "HBAR_J_S",
    "MICROELECTRONVOLT_TO_JOULE",
    "MICROSECOND_TO_SECOND",
    "NANOSECOND_TO_SECOND",
    "RAD_PER_NS_PER_MICROELECTRONVOLT",
    "RAD_PER_US_PER_MICROELECTRONVOLT",
    "TWO_PI",
    "build_H",
    "constant",
    "drift_H",
    "edsr1_op",
    "edsr2_op",
    "eps1_op",
    "eps2_op",
    "dqd_populations",
    "dqd_qubit_unitary",
    "DQDSequenceCompiler",
    "DQDsystem",
    "initial_full_state",
    "iswap_H",
    "piecewise",
    "plot_pulse_schedule",
    "print_summary",
    "Pulse",
    "run_state",
    "run_unitary",
    "SingleDQDOperatorsNumpy",
    "SingleDQDOperatorsQutip",
    "SingleDQDSpectrum",
    "single_dqd_eigensystem",
    "single_dqd_hamiltonian",
    "single_dqd_numpy_hamiltonian",
    "single_dqd_numpy_operators",
    "single_dqd_qutip_hamiltonian",
    "single_dqd_qutip_operators",
    "single_dqd_qubit_splitting",
    "single_dqd_spectrum",
    "single_dqd_tau_z_matrix_element",
    "square",
    "two_qubit_unitary",
    "vacuum_eigenstates",
    "vacuum_subspace_unitary",
    "xrot_H",
    "zrot_H",
    "zrot_delta_freq_GHz",
    "zrot_time_ns",
]
