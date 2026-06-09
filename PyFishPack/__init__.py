"""
PyFishPack - Python wrapper for FISHPACK Fortran library

FISHPACK is a collection of Fortran routines for solving separable
elliptic partial differential equations.
"""

__version__ = "0.1.0"

try:
    # Import the compiled fishpack extension
    from . import fishpack
    from .apps import (
        invert_3DOcean,
        invert_BrethertonHaidvogel,
        invert_Eliassen,
        invert_Fofonoff,
        invert_GillMatsuno,
        invert_GillMatsuno_test,
        invert_MultiGrid,
        invert_PV2D,
        invert_Poisson,
        invert_RefState,
        invert_RefStateSWM,
        invert_Stommel,
        invert_StommelArons,
        invert_StommelMunk,
        invert_Stommel_test,
        invert_geostrophic,
        invert_omega,
    )

    # Make commonly used functions available at package level
    # These will be available as PyFishPack.genbun, etc.
    from .fishpack import *

    __all__ = [name for name in dir(fishpack) if not name.startswith('_')]
    __all__.extend([
        "invert_3DOcean",
        "invert_BrethertonHaidvogel",
        "invert_Eliassen",
        "invert_Fofonoff",
        "invert_GillMatsuno",
        "invert_GillMatsuno_test",
        "invert_MultiGrid",
        "invert_PV2D",
        "invert_Poisson",
        "invert_RefState",
        "invert_RefStateSWM",
        "invert_Stommel",
        "invert_StommelArons",
        "invert_StommelMunk",
        "invert_Stommel_test",
        "invert_geostrophic",
        "invert_omega",
    ])

except ImportError as e:
    import warnings
    print(f"Warning: Could not import fishpack extension: {e}")
    print("Make sure the extension is compiled correctly.")
    warnings.warn(f"Could not import fishpack extension: {e}. "
                  f"Make sure the extension is compiled correctly.",
                  ImportWarning)

    # Fallback - define empty __all__
    __all__ = []

# Package metadata
__author__ = "Qianye Su"
__email__ = "suqianye2000@gmail.com"
__url__ = "https://github.com/QianyeSu/PyFishPack"
