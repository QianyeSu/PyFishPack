"""Equation-level inversion helpers built on the compiled Fishpack backend."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from . import fishpack


_SUPPORTED_BCS = {"fixed", "periodic"}
_SUPPORTED_SOR_BCS = {"fixed", "extend", "periodic"}
_UNDEF = -9.99e8
_DEFAULT_MPARAMS = {
    "f0": 1e-5,
    "beta": 2e-11,
    "N2": 2e-4,
    "D": 100.0,
    "depth": 100.0,
    "lambda": 1e-8,
    "c0": 8e-9,
    "c1": 8e-5,
    "A": 1.0,
    "B": 0.0,
    "C": 1.0,
    "R": 1.0,
    "A4": 0.0,
    "rho0": 1025.0,
    "epsilon": 1e-5,
    "Phi": 1.0,
    "k": 1.0,
    "ang0": 2e5,
    "Gamma": 1.0,
    "M0": 0.0,
    "C0": 1.0,
    "Rearth": 6371200.0,
    "Omega": 7.292e-5,
    "g": 9.80665,
}
_DEFAULT_IPARAMS = {
    "BCs": ("fixed", "fixed"),
    "undef": np.nan,
    "mxLoop": 5000,
    "tolerance": 1e-8,
    "optArg": None,
    "printInfo": False,
    "debug": False,
}


def invert_Poisson(
    F: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian Poisson equation on a uniform grid.

    This is the equation-level xinvert-style wrapper around the modern Fortran
    Fishpack ``genbun`` backend.  It supports only Cartesian, uniform-grid,
    constant-coefficient Poisson solves.  Lat-lon, variable-coefficient, and
    other non-Cartesian formulations are not supported.

    Parameters
    ----------
    F : array-like or xarray.DataArray
        Forcing field to invert.  NumPy inputs are solved along ``dims`` as
        axis indices; xarray inputs use named dimensions.
    dims : sequence of str or int, optional
        Two inversion dimensions.  Defaults to the last two dimensions.
    coords : str, default "cartesian"
        Coordinate system.  Only ``"cartesian"`` is supported.
    icbc : any, optional
        Accepted for xinvert-style compatibility and ignored.
    mParams, iParams : dict, optional
        Compatibility mappings for equation and inversion parameters.  Boundary
        conditions may be supplied through ``iParams["BCs"]`` when ``BCs`` is
        not passed explicitly.
    BCs : sequence of {"fixed", "periodic"}, optional
        Boundary conditions for the two inversion dimensions.  Defaults to
        ``("fixed", "fixed")``.
    spacing : sequence of float, optional
        Grid spacing ``(dy, dx)`` for the inversion dimensions.  Defaults to
        unit spacing.
    raise_on_error : bool, default True
        Raise ``RuntimeError`` if Fishpack reports a nonzero solver error code.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    del icbc, mParams  # Kept for API compatibility with xinvert-style calls.

    if coords != "cartesian":
        raise NotImplementedError(
            "invert_Poisson currently supports uniform Cartesian coordinates only"
        )

    params = dict(iParams or {})
    bcs = _normalize_bcs(BCs if BCs is not None else params.get("BCs", None))

    if _is_dataarray(F):
        return _invert_poisson_labeled(
            F, dims=dims, bcs=bcs, spacing=spacing, raise_on_error=raise_on_error
        )

    return _invert_poisson_ndarray(
        F, axes=dims, bcs=bcs, spacing=spacing, raise_on_error=raise_on_error
    )


def invert_RefState(
    PV: Any,
    dims: Sequence[str] | None,
    coords: str = "z-lat",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
) -> Any:
    r"""Invert xinvert's balanced symmetric-vortex reference-state equation.

    The solve is backed by the modern Fortran SOR kernel exposed as
    ``fishpack.sor_standard2d``.  The wrapper follows xinvert's coefficient
    construction for the ``"z-lat"`` and ``"cartesian"`` coordinate forms and
    preserves xarray coordinates and metadata.

    Parameters
    ----------
    PV : xarray.DataArray
        Two-dimensional potential-vorticity distribution, optionally with
        non-core dimensions.
    dims : sequence of str
        Two inversion dimensions.  For ``"z-lat"``, the second dimension is
        interpreted as latitude in degrees.  For ``"cartesian"``, the second
        dimension is interpreted as radius.
    coords : {"z-lat", "cartesian"}, default "z-lat"
        Coordinate form used for xinvert-compatible coefficients.
    icbc : xarray.DataArray, optional
        Initial guess and fixed boundary values.
    mParams : dict, optional
        Model parameters.  Uses ``ang0``, ``Gamma``, ``g``, and ``Rearth``.
    iParams : dict, optional
        Iteration parameters including ``BCs``, ``undef``, ``mxLoop``,
        ``tolerance``, and optional ``optArg``.
    BCs : sequence of {"fixed", "extend", "periodic"}, optional
        Boundary conditions for the two inversion dimensions.

    Returns
    -------
    xarray.DataArray
        Inverted angular-momentum field named ``"inverted"``.
    """

    _require_dataarray(PV, "invert_RefState")
    if dims is None or len(dims) != 2:
        raise ValueError("invert_RefState requires two xarray dimension names")
    if coords.lower() not in {"z-lat", "cartesian"}:
        raise NotImplementedError("invert_RefState supports only 'z-lat' and 'cartesian'")

    params = _merged_mparams(mParams)
    iparams = _merged_iparams(iParams, ndim=2, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 2)
    mask_f, init_s, zero = _mask_labeled_field(PV, dims, iparams, bcs, icbc)
    ydim, xdim = dims
    xcoord = mask_f.coords[xdim]

    if coords.lower() == "z-lat":
        acoef = zero + np.sin(np.deg2rad(xcoord))
        ccoef = zero + params["Gamma"] * float(params["g"]) / mask_f / xcoord
        dy, dx = _sor_spacing2d(mask_f, (ydim, xdim), coords, float(params["Rearth"]))
    else:
        acoef = zero + 2.0 * params["ang0"] / (xcoord**3.0)
        ccoef = zero + params["Gamma"] * float(params["g"]) / mask_f / xcoord
        dy, dx = _sor_spacing2d(mask_f, (ydim, xdim), coords, float(params["Rearth"]))

    bcoef = zero
    solved = _solve_sor2d_labeled(
        init_s,
        acoef,
        bcoef,
        ccoef,
        mask_f,
        dims=dims,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def invert_RefStateSWM(
    Q: Any,
    dims: Sequence[str] | None,
    coords: str = "lat",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
) -> Any:
    r"""Invert xinvert's one-dimensional shallow-water reference state.

    The variable-coefficient 1-D solve is performed by the modern Fortran SOR
    kernel ``fishpack.sor_standard1d``.  Coefficients follow xinvert's
    ``invert_RefStateSWM`` construction for latitude coordinates.

    Parameters
    ----------
    Q : xarray.DataArray
        Potential-vorticity contour field with one inversion dimension.
    dims : sequence of str
        Single latitude dimension.
    coords : {"lat"}, default "lat"
        Coordinate form.  Only xinvert's latitude form is supported.
    icbc : xarray.DataArray, optional
        Initial guess and fixed boundary values.
    mParams : dict, optional
        Model parameters.  Uses ``M0``, ``C0``, ``g``, ``Rearth``, and
        ``Omega``.
    iParams : dict, optional
        Iteration parameters including ``BCs``, ``undef``, ``mxLoop``,
        ``tolerance``, and optional ``optArg``.
    BCs : sequence of {"fixed", "extend", "periodic"}, optional
        Boundary condition for the latitude dimension.

    Returns
    -------
    xarray.DataArray
        Inverted mass-correction field named ``"inverted"``.
    """

    _require_dataarray(Q, "invert_RefStateSWM")
    if dims is None or len(dims) != 1:
        raise ValueError("invert_RefStateSWM requires one xarray dimension name")
    if coords.lower() != "lat":
        raise NotImplementedError("invert_RefStateSWM supports only 'lat' coordinates")

    params = _merged_mparams(mParams)
    iparams = _merged_iparams(iParams, ndim=1, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 1)
    mask_f, init_s, zero = _mask_labeled_field(Q, dims, iparams, bcs, icbc)
    dim = dims[0]
    lat = np.deg2rad(mask_f.coords[dim])
    cos_g = np.cos(lat)
    sin_g = np.sin(lat)
    cos_h = np.cos((lat + lat.shift({dim: 1})) / 2.0)
    asin = float(params["Rearth"]) * sin_g
    acos = float(params["Rearth"]) * cos_g
    acos = acos.where(acos >= 0.0, other=-acos * 0.1)

    m0 = _as_labeled_like(params["M0"], mask_f)
    c0 = _as_labeled_like(params["C0"], mask_f)
    delx = _sor_spacing1d(mask_f, dim, coords, float(params["Rearth"]))
    diff = _second_diff_swm(m0, cos_h, delx, dim)

    acoef = zero + 1.0 / cos_h
    bcoef = zero - c0 * mask_f * asin / (np.pi * float(params["g"]) * acos**3.0)
    force = (
        zero
        - (asin * c0**2.0 / (2.0 * np.pi * float(params["g"]) * acos**3.0))
        + (2.0 * np.pi * float(params["Omega"]) ** 2.0 * asin * acos) / float(params["g"])
        - diff
    )
    solved = _solve_sor1d_labeled(
        init_s,
        acoef,
        bcoef,
        force,
        dims=dims,
        dx=delx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def invert_geostrophic(
    lapPhi: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian geostrophic balance equation.

    The constant-Coriolis Cartesian subset uses the direct modern Fortran
    Fishpack path.  Cartesian beta-plane inputs use the modern Fortran
    ``sor_standard2d`` backend with xinvert-style flux coefficients.  Lat-lon
    formulations remain unsupported.

    Parameters
    ----------
    lapPhi : array-like or xarray.DataArray
        Laplacian forcing field to divide by ``f0`` before the 2-D inversion.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``f0`` and
        ``beta``; the direct ``beta = 0`` path requires nonzero ``f0``.

    Returns
    -------
    array-like or xarray.DataArray
        The geostrophic streamfunction or velocity potential field.
    """

    if coords != "cartesian":
        raise NotImplementedError(
            "invert_geostrophic currently supports Cartesian coordinates only"
        )
    params = _merged_mparams(mParams)
    beta = _scalar_param(params, "beta")
    f0 = _scalar_param(params, "f0")
    if beta != 0.0:
        return _invert_standard_2d(
            lapPhi,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=_cartesian_geostrophic_coefficients(
                lapPhi,
                dims=dims,
                spacing=spacing,
                f0=f0,
                beta=beta,
            ),
        )
    if f0 == 0.0:
        raise ValueError("f0 must be non-zero")
    return _invert_constant_helmholtz(
        np.asarray(lapPhi, dtype=np.float64) / f0 if not _is_dataarray(lapPhi) else lapPhi / f0,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def invert_PV2D(
    PV: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian scalar-``N2`` QG potential-vorticity equation.

    Only the Cartesian, uniform-grid, constant-coefficient subset is
    supported.

    Parameters
    ----------
    PV : array-like or xarray.DataArray
        Potential-vorticity forcing field.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``f0`` and
        ``N2``; ``N2`` must be a positive scalar.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    del icbc
    if coords != "cartesian":
        raise NotImplementedError("invert_PV2D currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    f0 = _scalar_param(params, "f0")
    n2 = _scalar_param(params, "N2")
    if n2 <= 0.0:
        raise ValueError("N2 must be a positive scalar")
    return _invert_constant_2d(
        PV,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(f0 * f0 / n2, 1.0),
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def invert_Eliassen(
    F: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian Eliassen equation.

    Cartesian inputs are supported for both NumPy arrays and
    :class:`xarray.DataArray` objects.  Non-Cartesian formulations remain
    unsupported.  The separable ``B = 0`` case uses the direct Fishpack /
    ``genbun`` path, while constant-coefficient ``B != 0`` cases route to the
    modern Fortran ``sor_general2d`` backend for the equivalent
    cross-derivative form.

    Parameters
    ----------
    F : array-like or xarray.DataArray
        Forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``A``, ``B``, and
        ``C``.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    if coords != "cartesian":
        raise NotImplementedError("invert_Eliassen currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    cross = _scalar_param(params, "B")
    if cross != 0.0:
        return _invert_general_2d(
            F,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=(
                _scalar_param(params, "A"),
                2.0 * cross,
                _scalar_param(params, "C"),
                0.0,
                0.0,
                0.0,
            ),
        )
    return _invert_constant_2d(
        F,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(_scalar_param(params, "A"), _scalar_param(params, "C")),
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def invert_Fofonoff(
    F: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian Fofonoff Helmholtz equation.

    Cartesian uniform-grid inputs are backed by the modern Fortran Fishpack
    Helmholtz path.  ``beta = 0`` uses a constant right-hand side; beta-plane
    inputs use the Cartesian y coordinate from xarray coordinates or, for
    NumPy arrays, from ``dims`` and ``spacing``.

    Parameters
    ----------
    F : array-like or xarray.DataArray
        Forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``f0``, ``beta``,
        ``c0``, and ``c1``.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    del icbc
    if coords != "cartesian":
        raise NotImplementedError("invert_Fofonoff currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    forcing = _cartesian_coriolis_forcing(
        F,
        dims,
        spacing,
        params,
        sign=-1.0,
        offset=float(params["c1"]),
    )
    return _invert_constant_helmholtz(
        forcing,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        helmholtz=-float(params["c0"]),
        raise_on_error=raise_on_error,
    )


def invert_GillMatsuno(
    Q: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian Gill-Matsuno Helmholtz subset.

    Cartesian inputs are supported for both NumPy arrays and
    :class:`xarray.DataArray` objects.  Non-Cartesian formulations remain
    unsupported.  The ``beta = 0`` case uses the direct Fishpack / ``genbun``
    path, while the Cartesian beta-plane case routes to the modern Fortran
    ``sor_general2d`` backend.

    Parameters
    ----------
    Q : array-like or xarray.DataArray
        Forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``epsilon``,
        ``Phi``, ``f0``, and ``beta``.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    if coords != "cartesian":
        raise NotImplementedError("invert_GillMatsuno currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    beta = _scalar_param(params, "beta")
    epsilon = _scalar_param(params, "epsilon")
    phi = _scalar_param(params, "Phi")
    f0 = _scalar_param(params, "f0")
    denom = epsilon * epsilon + f0 * f0
    if denom == 0.0:
        raise ValueError("epsilon and f0 cannot both be zero")
    if beta != 0.0:
        return _invert_general_2d(
            Q,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=_cartesian_beta_general_coefficients(
                Q,
                dims=dims,
                spacing=spacing,
                epsilon=epsilon,
                f0=f0,
                beta=beta,
                scale=phi,
                helmholtz=-epsilon,
            ),
        )
    alpha = epsilon * phi / denom
    return _invert_constant_2d(
        Q,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(alpha, alpha),
        helmholtz=-epsilon,
        raise_on_error=raise_on_error,
    )


def invert_GillMatsuno_test(
    Q: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the xinvert Gill-Matsuno test-form subset for ``beta = 0``.

    This wrapper exposes the Cartesian, uniform-grid, constant-coefficient test
    subset backed by the modern Fortran Fishpack implementation.  Beta-plane,
    non-Cartesian, and other nonseparable cases are not supported here and
    raise through the delegated implementation.

    Parameters
    ----------
    Q : array-like or xarray.DataArray
        Forcing field to invert.
    dims : sequence of str or int, optional
        Inversion dimensions.  Defaults to the last two dimensions.
    coords : str, default "cartesian"
        Coordinate system.  Only ``"cartesian"`` is supported.
    icbc : any, optional
        Accepted for xinvert-style compatibility and ignored.
    mParams : dict, optional
        Equation-parameter mapping.  This subset expects the usual
        Gill-Matsuno parameters and requires ``beta = 0``.
    iParams : dict, optional
        Inversion-parameter mapping.  Boundary conditions may be supplied
        through ``iParams["BCs"]`` when ``BCs`` is not passed explicitly.
    BCs : sequence of {"fixed", "periodic"}, optional
        Boundary conditions for the two inversion dimensions.
    spacing : sequence of float, optional
        Grid spacing for the inversion dimensions.  Defaults to unit spacing.
    raise_on_error : bool, default True
        Raise ``RuntimeError`` if Fishpack reports a nonzero solver error code.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    return invert_GillMatsuno(
        Q,
        dims=dims,
        coords=coords,
        icbc=icbc,
        mParams=mParams,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        raise_on_error=raise_on_error,
    )


def invert_BrethertonHaidvogel(
    h: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian constant-depth Bretherton-Haidvogel equation.

    The constant-depth Cartesian subset is backed by the modern Fortran
    Fishpack Helmholtz path.  ``beta = 0`` uses a constant Coriolis parameter;
    beta-plane inputs use the Cartesian y coordinate from xarray coordinates
    or, for NumPy arrays, from ``dims`` and ``spacing``.

    Parameters
    ----------
    h : array-like or xarray.DataArray
        Layer-thickness forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``f0``, ``beta``,
        and ``lambda``; ``D`` or ``depth`` must be nonzero.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    del icbc
    if coords != "cartesian":
        raise NotImplementedError(
            "invert_BrethertonHaidvogel currently supports Cartesian coordinates only"
        )
    params = _merged_mparams(mParams)
    depth = float(params.get("D", params["depth"]))
    if depth == 0.0:
        raise ValueError("D/depth must be non-zero")
    coriolis = _cartesian_coriolis_field(h, dims, spacing, params)
    if _is_dataarray(h):
        forcing = -(h * coriolis) / depth
    else:
        forcing = -(np.asarray(h, dtype=np.float64) * coriolis) / depth
    return _invert_constant_helmholtz(
        forcing,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        helmholtz=-float(params["lambda"]) * depth,
        raise_on_error=raise_on_error,
    )


def invert_Stommel(
    curl: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian Stommel equation on the supported subset.

    The Cartesian ``beta = 0`` subset uses the direct Fishpack / ``genbun``
    solve path.  When ``beta != 0``, the solve is routed to the modern
    Fortran ``sor_general2d`` backend.  Non-Cartesian formulations remain
    unsupported.

    Parameters
    ----------
    curl : array-like or xarray.DataArray
        Wind-stress curl forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``D``, ``rho0``,
        ``R``, and ``beta``; ``D`` and ``rho0`` must be nonzero.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    if coords != "cartesian":
        raise NotImplementedError("invert_Stommel currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    beta = _scalar_param(params, "beta")
    depth = _scalar_param(params, "D")
    rho0 = _scalar_param(params, "rho0")
    resistance = _scalar_param(params, "R")
    if depth == 0.0 or rho0 == 0.0:
        raise ValueError("D and rho0 must be non-zero")
    alpha = -resistance / depth
    forcing = curl * (-1.0 / (depth * rho0)) if _is_dataarray(curl) else np.asarray(curl, dtype=np.float64) * (-1.0 / (depth * rho0))
    if beta != 0.0:
        return _invert_general_2d(
            forcing,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=(alpha, 0.0, alpha, 0.0, -beta, 0.0),
        )
    return _invert_constant_2d(
        forcing,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(alpha, alpha),
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def invert_Stommel_test(
    curl: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the xinvert Stommel test-form subset.

    This wrapper preserves xinvert's Cartesian, uniform-grid,
    constant-coefficient test form and delegates to :func:`invert_Stommel`.
    As a result, ``beta = 0`` uses the direct Fishpack / ``genbun`` path while
    ``beta != 0`` uses the modern Fortran ``sor_general2d`` backend.
    Non-Cartesian formulations remain unsupported.

    Parameters
    ----------
    curl : array-like or xarray.DataArray
        Wind-stress curl forcing field to invert.
    dims : sequence of str or int, optional
        Inversion dimensions.  Defaults to the last two dimensions.
    coords : str, default "cartesian"
        Coordinate system.  Only ``"cartesian"`` is supported.
    icbc : any, optional
        Accepted for xinvert-style compatibility and ignored.
    mParams : dict, optional
        Equation-parameter mapping.  This subset expects the usual Stommel
        parameters; ``D`` and ``rho0`` must be nonzero, and the ``beta``
        handling follows :func:`invert_Stommel`.
    iParams : dict, optional
        Inversion-parameter mapping.  Boundary conditions may be supplied
        through ``iParams["BCs"]`` when ``BCs`` is not passed explicitly.
    BCs : sequence of {"fixed", "periodic"}, optional
        Boundary conditions for the two inversion dimensions.
    spacing : sequence of float, optional
        Grid spacing for the inversion dimensions.  Defaults to unit spacing.
    raise_on_error : bool, default True
        Raise ``RuntimeError`` if Fishpack reports a nonzero solver error code.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    return invert_Stommel(
        curl,
        dims=dims,
        coords=coords,
        icbc=icbc,
        mParams=mParams,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        raise_on_error=raise_on_error,
    )


def invert_StommelMunk(
    curl: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian Stommel-Munk equation in xinvert form.

    The ``A4 = 0`` subset delegates to :func:`invert_Stommel`.  Nonzero
    ``A4`` uses the modern Fortran ``sor_biharmonic2d`` backend for the
    Cartesian fourth-order Stommel-Munk equation.  Lat-lon and other
    non-Cartesian formulations remain unsupported.

    Parameters
    ----------
    curl : array-like or xarray.DataArray
        Wind-stress curl forcing field to invert.
    dims : sequence of str or int, optional
        Inversion dimensions.  Defaults to the last two dimensions.
    coords : str, default "cartesian"
        Coordinate system.  Only ``"cartesian"`` is supported.
    icbc : any, optional
        Accepted for xinvert-style compatibility and ignored.
    mParams : dict, optional
        Equation-parameter mapping.  This wrapper uses ``A4``, ``beta``,
        ``D``, ``rho0``, and ``R``.
    iParams : dict, optional
        Inversion-parameter mapping.  Boundary conditions may be supplied
        through ``iParams["BCs"]`` when ``BCs`` is not passed explicitly.
    BCs : sequence of {"fixed", "periodic"}, optional
        Boundary conditions for the two inversion dimensions.
    spacing : sequence of float, optional
        Grid spacing for the inversion dimensions.  Defaults to unit spacing.
    raise_on_error : bool, default True
        Raise ``RuntimeError`` if Fishpack reports a nonzero solver error code.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    params = _merged_mparams(mParams)
    a4 = _scalar_param(params, "A4")
    if a4 != 0.0:
        if coords != "cartesian":
            raise NotImplementedError("invert_StommelMunk currently supports Cartesian coordinates only")
        depth = _scalar_param(params, "D")
        rho0 = _scalar_param(params, "rho0")
        if depth == 0.0 or rho0 == 0.0:
            raise ValueError("D and rho0 must be non-zero")
        forcing = (
            curl * (-1.0 / (depth * rho0))
            if _is_dataarray(curl)
            else np.asarray(curl, dtype=np.float64) * (-1.0 / (depth * rho0))
        )
        resistance = _scalar_param(params, "R")
        beta = _scalar_param(params, "beta")
        return _invert_biharmonic_2d(
            forcing,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=(a4, 0.0, a4, -resistance / depth, 0.0, -resistance / depth, 0.0, -beta, 0.0),
        )
    return invert_Stommel(
        curl,
        dims=dims,
        coords=coords,
        icbc=icbc,
        mParams=params,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        raise_on_error=raise_on_error,
    )


def invert_StommelArons(
    Q: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian constant-Coriolis Stommel-Arons subset.

    Cartesian inputs are supported for both NumPy arrays and
    :class:`xarray.DataArray` objects.  Non-Cartesian formulations remain
    unsupported.  The ``beta = 0`` case uses the direct Fishpack / ``genbun``
    path, while the Cartesian beta-plane case routes to the modern Fortran
    ``sor_general2d`` backend.

    Parameters
    ----------
    Q : array-like or xarray.DataArray
        Forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``epsilon``,
        ``f0``, and ``beta``.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    if coords != "cartesian":
        raise NotImplementedError("invert_StommelArons currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    beta = _scalar_param(params, "beta")
    epsilon = _scalar_param(params, "epsilon")
    f0 = _scalar_param(params, "f0")
    denom = epsilon * epsilon + f0 * f0
    if denom == 0.0:
        raise ValueError("epsilon and f0 cannot both be zero")
    if beta != 0.0:
        return _invert_general_2d(
            Q,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=_cartesian_beta_general_coefficients(
                Q,
                dims=dims,
                spacing=spacing,
                epsilon=epsilon,
                f0=f0,
                beta=beta,
                scale=1.0,
                helmholtz=0.0,
            ),
        )
    alpha = epsilon / denom
    return _invert_constant_2d(
        Q,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(alpha, alpha),
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def invert_omega(
    F: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian QG omega equation.

    The constant-Coriolis subset uses the direct Fishpack 3-D solver.  The
    Cartesian beta-plane subset uses the modern Fortran ``sor_standard3d``
    backend with xinvert-style flux coefficients.

    Parameters
    ----------
    F : array-like or xarray.DataArray
        Forcing field to invert.
    dims, coords, icbc, mParams, iParams, BCs, spacing, raise_on_error
        See :func:`invert_Poisson`.  ``mParams`` may provide ``f0``, ``beta``,
        and ``N2``; this subset requires ``beta = 0`` and positive ``N2``.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    if coords != "cartesian":
        raise NotImplementedError("invert_omega currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    beta = float(params["beta"])
    n2 = float(params["N2"])
    if n2 <= 0.0:
        raise ValueError("N2 must be a positive scalar")
    f0 = float(params["f0"])
    if beta != 0.0:
        return _invert_standard_3d(
            F,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=_cartesian_omega_coefficients(
                F, dims=dims, spacing=spacing, f0=f0, beta=beta, n2=n2
            ),
        )
    del icbc
    rhs = F / n2 if _is_dataarray(F) else np.asarray(F, dtype=np.float64) / n2
    return _invert_constant_3d(
        rhs,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(f0 * f0 / n2, 1.0, 1.0),
        raise_on_error=raise_on_error,
    )


def invert_3DOcean(
    F: Any,
    dims: Sequence[str] | Sequence[int] | None = None,
    coords: str = "cartesian",
    icbc: Any = None,
    mParams: dict[str, Any] | None = None,
    iParams: dict[str, Any] | None = None,
    *,
    BCs: Sequence[str] | None = None,
    spacing: Sequence[float] | None = None,
    raise_on_error: bool = True,
) -> Any:
    r"""Invert the Cartesian constant-coefficient 3-D ocean equation subset.

    Only the Cartesian, uniform-grid, constant-coefficient Fishpack subset is
    supported.  Beta-plane and other nonseparable formulations intentionally
    raise ``NotImplementedError``.

    Parameters
    ----------
    F : array-like or xarray.DataArray
        Forcing field to invert.
    dims : sequence of str or int, optional
        Three inversion dimensions.  Defaults to the last three dimensions.
    coords : str, default "cartesian"
        Coordinate system.  Only ``"cartesian"`` is supported.
    icbc : any, optional
        Accepted for xinvert-style compatibility and ignored.
    mParams, iParams : dict, optional
        Compatibility mappings for equation and inversion parameters.  Boundary
        conditions may be supplied through ``iParams["BCs"]`` when ``BCs`` is
        not passed explicitly.
    BCs : sequence of {"fixed", "periodic"}, optional
        Boundary conditions for the three inversion dimensions.  Defaults to
        ``("fixed", "fixed", "fixed")``.
    spacing : sequence of float, optional
        Grid spacing ``(dz, dy, dx)`` for the inversion dimensions.  Defaults
        to unit spacing.
    raise_on_error : bool, default True
        Raise ``RuntimeError`` if Fishpack reports a nonzero solver error code.

    Returns
    -------
    array-like or xarray.DataArray
        The inverted field with the same array type as the input.
    """

    if coords != "cartesian":
        raise NotImplementedError("invert_3DOcean currently supports Cartesian coordinates only")
    params = _merged_mparams(mParams)
    beta = _scalar_param(params, "beta")
    epsilon = _scalar_param(params, "epsilon")
    f0 = _scalar_param(params, "f0")
    n2 = _scalar_param(params, "N2")
    buoyancy_damping = _scalar_param(params, "k")
    if n2 <= 0.0:
        raise ValueError("N2 must be a positive scalar")
    denom = epsilon * epsilon + f0 * f0
    if denom == 0.0:
        raise ValueError("epsilon and f0 cannot both be zero")
    if beta != 0.0:
        return _invert_general_3d(
            F,
            dims=dims,
            coords=coords,
            iParams=iParams,
            BCs=BCs,
            spacing=spacing,
            icbc=icbc,
            coefficients=_cartesian_3d_ocean_coefficients(
                F,
                dims=dims,
                spacing=spacing,
                epsilon=epsilon,
                f0=f0,
                beta=beta,
                n2=n2,
                buoyancy_damping=buoyancy_damping,
            ),
        )
    del icbc
    horizontal = epsilon / denom
    vertical = buoyancy_damping / n2
    return _invert_constant_3d(
        F,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(vertical, horizontal, horizontal),
        raise_on_error=raise_on_error,
    )


def invert_MultiGrid(
    invert_func: Any,
    *args: Any,
    ratio: int = 3,
    gridNo: int = 3,
    **kwargs: Any,
) -> tuple[Any, list[list[Any]], list[Any]]:
    """Compatibility helper for xinvert's ``MultiGrid`` entry point.

    This is not a separate equation solver. PyFishPack uses direct Fishpack
    solvers, so this wrapper delegates once to ``invert_func`` and returns a
    xinvert-style ``(solution, grids, history)`` tuple.

    Parameters
    ----------
    invert_func : callable
        Solver function to call once with ``*args`` and ``**kwargs``.
    *args : Any
        Positional arguments passed through to ``invert_func``.
    ratio : int, optional
        Accepted for API compatibility with xinvert and currently ignored.
    gridNo : int, optional
        Accepted for API compatibility with xinvert and currently ignored.
    **kwargs : Any
        Keyword arguments passed through to ``invert_func``.

    Returns
    -------
    tuple[Any, list[list[Any]], list[Any]]
        ``(solution, grids, history)`` where ``solution`` is the direct
        result from ``invert_func``, ``grids`` records the input arguments, and
        ``history`` records the single returned solution.
    """

    del ratio, gridNo
    if not callable(invert_func):
        raise TypeError("invert_func must be callable")
    solution = invert_func(*args, **kwargs)
    return solution, [list(args)], [solution]


def spectral_transform(
    data: Any,
    *,
    kind: str = "rfft",
    direction: str = "forward",
    axis: int = -1,
    normalize: bool = False,
) -> np.ndarray:
    """Apply a lightweight one-dimensional FFTPACK transform along an array axis."""

    transform = kind.lower()
    direct = direction.lower()
    if direct not in {"forward", "backward", "inverse"}:
        raise ValueError("direction must be 'forward', 'backward', or 'inverse'")
    inverse = direct in {"backward", "inverse"}

    real_methods = {
        ("rfft", False): fishpack.rfftf,
        ("rfft", True): fishpack.rfftb,
        ("sint", False): fishpack.sint,
        ("sint", True): fishpack.sint,
        ("cost", False): fishpack.cost,
        ("cost", True): fishpack.cost,
        ("sinq", False): fishpack.sinqf,
        ("sinq", True): fishpack.sinqb,
        ("cosq", False): fishpack.cosqf,
        ("cosq", True): fishpack.cosqb,
    }

    if transform == "cfft":
        arr = np.asarray(data, dtype=np.complex128)
        if arr.ndim == 0:
            raise ValueError("spectral_transform expects at least one-dimensional input")
        axis = axis % arr.ndim
        moved = np.moveaxis(arr, axis, -1)
        flat = moved.reshape((-1, moved.shape[-1]))
        out = np.empty_like(flat)
        method = fishpack.cfftb if inverse else fishpack.cfftf
        for idx, row in enumerate(flat):
            interleaved = np.ascontiguousarray(row, dtype=np.complex128).view(np.float64)
            transformed = np.asarray(method(interleaved), dtype=np.float64).view(np.complex128)
            out[idx] = transformed
        result = out.reshape(moved.shape)
        if normalize and inverse:
            result = result / moved.shape[-1]
        return np.moveaxis(result, -1, axis)

    method = real_methods.get((transform, inverse))
    if method is None:
        raise ValueError("kind must be one of 'rfft', 'cfft', 'sint', 'cost', 'sinq', or 'cosq'")

    arr = np.asarray(data, dtype=np.float64)
    if arr.ndim == 0:
        raise ValueError("spectral_transform expects at least one-dimensional input")
    axis = axis % arr.ndim
    moved = np.moveaxis(arr, axis, -1)
    flat = moved.reshape((-1, moved.shape[-1]))
    out = np.empty_like(flat)
    for idx, row in enumerate(flat):
        out[idx] = method(np.ascontiguousarray(row, dtype=np.float64))
    result = out.reshape(moved.shape)
    if normalize and inverse and transform == "rfft":
        result = result / moved.shape[-1]
    return np.moveaxis(result, -1, axis)


def _require_dataarray(field: Any, func_name: str) -> None:
    if not _is_dataarray(field):
        raise NotImplementedError(f"{func_name} currently requires an xarray.DataArray input")


def _merged_iparams(
    params: dict[str, Any] | None,
    *,
    ndim: int,
    BCs: Sequence[str] | None,
) -> dict[str, Any]:
    merged = dict(_DEFAULT_IPARAMS)
    merged["BCs"] = tuple("fixed" for _ in range(ndim))
    if params is not None:
        merged.update(params)
    if BCs is not None:
        merged["BCs"] = tuple(BCs)
    return merged


def _mask_labeled_field(
    field: Any,
    dims: Sequence[str],
    iparams: dict[str, Any],
    bcs: tuple[str, ...],
    icbc: Any,
) -> tuple[Any, Any, Any]:
    missing = [dim for dim in dims if dim not in field.dims]
    if missing:
        raise ValueError(f"dims not present in input DataArray: {missing}")

    undef = iparams["undef"]
    if np.isnan(undef):
        mask_f = field.fillna(_UNDEF)
    else:
        mask_f = field.where(field != undef, other=_UNDEF)
    zero = mask_f - mask_f

    if icbc is None:
        init_s = zero.copy()
    else:
        mask = mask_f == _UNDEF
        for dim, bc in zip(dims, bcs):
            if bc != "periodic":
                coord = mask_f.coords[dim]
                cond = coord.isin([coord[0], coord[-1]])
                mask = np.logical_or(mask, cond)
        init_s = field.__class__(icbc, coords=field.coords, dims=field.dims) if not _is_dataarray(icbc) else icbc
        init_s = init_s.where(mask, other=0)
    return mask_f, init_s.load(), zero


def _restore_labeled_result(result: Any, mask_f: Any, iparams: dict[str, Any], icbc: Any) -> Any:
    result = result.rename("inverted")
    if icbc is None:
        return result.where(mask_f != _UNDEF, other=iparams["undef"])
    return result


def _solve_sor2d_labeled(
    init_s: Any,
    acoef: Any,
    bcoef: Any,
    ccoef: Any,
    force: Any,
    *,
    dims: Sequence[str],
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> Any:
    import xarray as xr

    ydim, xdim = dims
    acoef, bcoef, ccoef, force, init_s = xr.broadcast(acoef, bcoef, ccoef, force, init_s)
    outer_dims = tuple(dim for dim in force.dims if dim not in dims)
    order = (*outer_dims, ydim, xdim)
    a_t = acoef.transpose(*order)
    b_t = bcoef.transpose(*order)
    c_t = ccoef.transpose(*order)
    f_t = force.transpose(*order)
    s_t = init_s.transpose(*order)
    values = np.empty(s_t.shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, s_t.shape[-2:])

    for index in np.ndindex(s_t.shape[:-2] or ()):
        key = index if s_t.ndim > 2 else (...,)
        if s_t.ndim == 2:
            solved, relerr, overflow, loops = fishpack.sor_standard2d(
                s_t.values,
                a_t.values,
                b_t.values,
                c_t.values,
                f_t.values,
                dy,
                dx,
                bcs[0],
                bcs[1],
                optarg,
                _UNDEF,
                int(iparams["mxLoop"]),
                float(iparams["tolerance"]),
            )
            values[...] = solved
            if overflow:
                raise RuntimeError("Fortran SOR standard2d overflowed")
            break
        solved, relerr, overflow, loops = fishpack.sor_standard2d(
            s_t.values[key],
            a_t.values[key],
            b_t.values[key],
            c_t.values[key],
            f_t.values[key],
            dy,
            dx,
            bcs[0],
            bcs[1],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        values[key] = solved
        if overflow:
            raise RuntimeError("Fortran SOR standard2d overflowed")

    result = force.__class__(
        values,
        coords=s_t.coords,
        dims=s_t.dims,
        attrs=dict(force.attrs),
        name="inverted",
    )
    return result.transpose(*force.dims)


def _solve_sor3d_labeled(
    init_s: Any,
    acoef: Any,
    bcoef: Any,
    ccoef: Any,
    force: Any,
    *,
    dims: Sequence[str],
    dz: float,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> Any:
    import xarray as xr

    zdim, ydim, xdim = dims
    acoef, bcoef, ccoef, force, init_s = xr.broadcast(acoef, bcoef, ccoef, force, init_s)
    outer_dims = tuple(dim for dim in force.dims if dim not in dims)
    order = (*outer_dims, zdim, ydim, xdim)
    a_t = acoef.transpose(*order)
    b_t = bcoef.transpose(*order)
    c_t = ccoef.transpose(*order)
    f_t = force.transpose(*order)
    s_t = init_s.transpose(*order)
    values = np.empty(s_t.shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, s_t.shape[-3:])

    if s_t.ndim == 3:
        solved, relerr, overflow, loops = fishpack.sor_standard3d(
            s_t.values,
            a_t.values,
            b_t.values,
            c_t.values,
            f_t.values,
            dz,
            dy,
            dx,
            bcs[0],
            bcs[1],
            bcs[2],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        values[...] = solved
        if overflow:
            raise RuntimeError("Fortran SOR standard3d overflowed")
    else:
        for index in np.ndindex(s_t.shape[:-3]):
            solved, relerr, overflow, loops = fishpack.sor_standard3d(
                s_t.values[index],
                a_t.values[index],
                b_t.values[index],
                c_t.values[index],
                f_t.values[index],
                dz,
                dy,
                dx,
                bcs[0],
                bcs[1],
                bcs[2],
                optarg,
                _UNDEF,
                int(iparams["mxLoop"]),
                float(iparams["tolerance"]),
            )
            values[index] = solved
            if overflow:
                raise RuntimeError("Fortran SOR standard3d overflowed")

    result = force.__class__(
        values,
        coords=s_t.coords,
        dims=s_t.dims,
        attrs=dict(force.attrs),
        name="inverted",
    )
    return result.transpose(*force.dims)


def _invert_standard_2d(
    F: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    icbc: Any,
    coefficients: tuple[Any, Any, Any],
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("standard-form SOR currently supports Cartesian coordinates only")
    iparams = _merged_iparams(iParams, ndim=2, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 2)
    if _is_dataarray(F):
        return _invert_standard_2d_labeled(
            F,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            iparams=iparams,
            icbc=icbc,
            coefficients=coefficients,
        )
    return _invert_standard_2d_ndarray(
        F,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        iparams=iparams,
        icbc=icbc,
        coefficients=coefficients,
    )


def _invert_standard_2d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any],
) -> Any:
    if dims is None:
        if field.ndim < 2:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
        dims = field.dims[-2:]
    if len(dims) != 2 or not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray standard-form inversion requires two dimension names")
    mask_f, init_s, _zero = _mask_labeled_field(field, dims, iparams, bcs, icbc)
    dy, dx = _spacing_for_labeled(mask_f, (dims[0], dims[1]), spacing)
    coefs = tuple(_as_labeled_like(coef, mask_f) for coef in coefficients)
    solved = _solve_sor2d_labeled(
        init_s,
        coefs[0],
        coefs[1],
        coefs[2],
        mask_f,
        dims=dims,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def _invert_standard_2d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any],
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("standard-form inversion requires at least a two-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 2 or not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy standard-form inversion requires two integer axes")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if axes_tuple[0] == axes_tuple[1]:
        raise ValueError("inversion axes must be distinct")
    dy, dx = _normalize_spacing(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-2, -1))
    init = (
        np.zeros_like(moved)
        if icbc is None
        else np.moveaxis(np.asarray(icbc, dtype=np.float64), axes_tuple, (-2, -1))
    )
    moved_coefficients = []
    for coef in coefficients:
        coef_arr = np.asarray(coef, dtype=np.float64)
        if coef_arr.ndim == arr.ndim and coef_arr.shape == arr.shape:
            coef_arr = np.moveaxis(coef_arr, axes_tuple, (-2, -1))
        moved_coefficients.append(coef_arr)
    coef_arrays = tuple(
        _broadcast_ndarray_coefficient(coef, moved.shape) for coef in moved_coefficients
    )
    solved = _solve_sor_standard2d_batched(
        init,
        coef_arrays[0],
        coef_arrays[1],
        coef_arrays[2],
        moved,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return np.moveaxis(solved, (-2, -1), axes_tuple)


def _invert_standard_3d(
    F: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    icbc: Any,
    coefficients: tuple[Any, Any, Any],
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("standard-form 3D SOR currently supports Cartesian coordinates only")
    iparams = _merged_iparams(iParams, ndim=3, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 3)
    if bcs[0] == "periodic" or bcs[1] == "periodic":
        raise NotImplementedError("standard-form 3D SOR currently supports periodic boundaries only in x")
    if _is_dataarray(F):
        return _invert_standard_3d_labeled(
            F,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            iparams=iparams,
            icbc=icbc,
            coefficients=coefficients,
        )
    return _invert_standard_3d_ndarray(
        F,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        iparams=iparams,
        icbc=icbc,
        coefficients=coefficients,
    )


def _invert_standard_3d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any],
) -> Any:
    if dims is None:
        if field.ndim < 3:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 3 dimensions")
        dims = field.dims[-3:]
    if len(dims) != 3 or not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray standard-form 3D inversion requires three dimension names")
    mask_f, init_s, _zero = _mask_labeled_field(field, dims, iparams, bcs, icbc)
    dz, dy, dx = _spacing_for_labeled3d(mask_f, (dims[0], dims[1], dims[2]), spacing)
    coefs = tuple(_as_labeled_like(coef, mask_f) for coef in coefficients)
    solved = _solve_sor3d_labeled(
        init_s,
        coefs[0],
        coefs[1],
        coefs[2],
        mask_f,
        dims=dims,
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def _invert_standard_3d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any],
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 3:
        raise ValueError("standard-form 3D inversion requires at least a three-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 3, arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 3 or not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy standard-form 3D inversion requires three integer axes")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if len(set(axes_tuple)) != 3:
        raise ValueError("inversion axes must be distinct")

    dz, dy, dx = _normalize_spacing3d(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-3, -2, -1))
    init = (
        np.zeros_like(moved)
        if icbc is None
        else np.moveaxis(np.asarray(icbc, dtype=np.float64), axes_tuple, (-3, -2, -1))
    )
    moved_coefficients = []
    for coef in coefficients:
        coef_arr = np.asarray(coef, dtype=np.float64)
        if coef_arr.ndim == arr.ndim and coef_arr.shape == arr.shape:
            coef_arr = np.moveaxis(coef_arr, axes_tuple, (-3, -2, -1))
        moved_coefficients.append(coef_arr)
    coef_arrays = tuple(
        _broadcast_ndarray_coefficient(coef, moved.shape) for coef in moved_coefficients
    )
    solved = _solve_sor_standard3d_batched(
        init,
        coef_arrays[0],
        coef_arrays[1],
        coef_arrays[2],
        moved,
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return np.moveaxis(solved, (-3, -2, -1), axes_tuple)


def _solve_sor_standard2d_batched(
    init_s: np.ndarray,
    acoef: np.ndarray,
    bcoef: np.ndarray,
    ccoef: np.ndarray,
    force: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> np.ndarray:
    arrays = [np.asarray(item, dtype=np.float64) for item in (init_s, acoef, bcoef, ccoef, force)]
    shape = np.broadcast_shapes(*(item.shape for item in arrays))
    arrays = [np.broadcast_to(item, shape) for item in arrays]
    if len(shape) < 2:
        raise ValueError("standard-form SOR expects the last two dimensions to be spatial")
    values = np.empty(shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, shape[-2:])
    if len(shape) == 2:
        solved, relerr, overflow, loops = fishpack.sor_standard2d(
            arrays[0],
            arrays[1],
            arrays[2],
            arrays[3],
            arrays[4],
            dy,
            dx,
            bcs[0],
            bcs[1],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        if overflow:
            raise RuntimeError("Fortran SOR standard2d overflowed")
        values[...] = solved
        return values
    for index in np.ndindex(shape[:-2]):
        solved, relerr, overflow, loops = fishpack.sor_standard2d(
            arrays[0][index],
            arrays[1][index],
            arrays[2][index],
            arrays[3][index],
            arrays[4][index],
            dy,
            dx,
            bcs[0],
            bcs[1],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        if overflow:
            raise RuntimeError("Fortran SOR standard2d overflowed")
        values[index] = solved
    return values


def _solve_sor_standard3d_batched(
    init_s: np.ndarray,
    acoef: np.ndarray,
    bcoef: np.ndarray,
    ccoef: np.ndarray,
    force: np.ndarray,
    *,
    dz: float,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> np.ndarray:
    arrays = [np.asarray(item, dtype=np.float64) for item in (init_s, acoef, bcoef, ccoef, force)]
    shape = np.broadcast_shapes(*(item.shape for item in arrays))
    arrays = [np.broadcast_to(item, shape) for item in arrays]
    if len(shape) < 3:
        raise ValueError("standard-form 3D SOR expects the last three dimensions to be spatial")
    values = np.empty(shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, shape[-3:])
    if len(shape) == 3:
        solved, relerr, overflow, loops = fishpack.sor_standard3d(
            arrays[0],
            arrays[1],
            arrays[2],
            arrays[3],
            arrays[4],
            dz,
            dy,
            dx,
            bcs[0],
            bcs[1],
            bcs[2],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        if overflow:
            raise RuntimeError("Fortran SOR standard3d overflowed")
        values[...] = solved
        return values
    for index in np.ndindex(shape[:-3]):
        solved, relerr, overflow, loops = fishpack.sor_standard3d(
            arrays[0][index],
            arrays[1][index],
            arrays[2][index],
            arrays[3][index],
            arrays[4][index],
            dz,
            dy,
            dx,
            bcs[0],
            bcs[1],
            bcs[2],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        if overflow:
            raise RuntimeError("Fortran SOR standard3d overflowed")
        values[index] = solved
    return values


def _invert_general_2d(
    G: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any],
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("general-form SOR currently supports Cartesian coordinates only")
    iparams = _merged_iparams(iParams, ndim=2, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 2)
    if _is_dataarray(G):
        return _invert_general_2d_labeled(
            G,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            iparams=iparams,
            icbc=icbc,
            coefficients=coefficients,
        )
    return _invert_general_2d_ndarray(
        G,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        iparams=iparams,
        icbc=icbc,
        coefficients=coefficients,
    )


def _invert_general_2d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any],
) -> Any:
    if dims is None:
        if field.ndim < 2:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
        dims = field.dims[-2:]
    if len(dims) != 2 or not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray general-form inversion requires two dimension names")
    mask_f, init_s, zero = _mask_labeled_field(field, dims, iparams, bcs, icbc)
    dy, dx = _spacing_for_labeled(mask_f, (dims[0], dims[1]), spacing)
    coefs = tuple(_as_labeled_like(coef, mask_f) for coef in coefficients)
    solved = _solve_sor_general2d_labeled(
        init_s,
        *coefs,
        mask_f,
        dims=dims,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def _invert_general_2d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any],
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("general-form inversion requires at least a two-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 2 or not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy general-form inversion requires two integer axes")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if axes_tuple[0] == axes_tuple[1]:
        raise ValueError("inversion axes must be distinct")
    dy, dx = _normalize_spacing(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-2, -1))
    init = np.zeros_like(moved) if icbc is None else np.moveaxis(np.asarray(icbc, dtype=np.float64), axes_tuple, (-2, -1))
    moved_coefficients = []
    for coef in coefficients:
        coef_arr = np.asarray(coef, dtype=np.float64)
        if coef_arr.ndim == arr.ndim and coef_arr.shape == arr.shape:
            coef_arr = np.moveaxis(coef_arr, axes_tuple, (-2, -1))
        moved_coefficients.append(coef_arr)
    coef_arrays = tuple(_broadcast_ndarray_coefficient(coef, moved.shape) for coef in moved_coefficients)
    solved = _solve_sor_general2d_batched(
        init,
        *coef_arrays,
        moved,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return np.moveaxis(solved, (-2, -1), axes_tuple)


def _invert_general_3d(
    H: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any, Any],
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("general-form 3D SOR currently supports Cartesian coordinates only")
    iparams = _merged_iparams(iParams, ndim=3, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 3)
    if bcs[0] == "periodic" or bcs[1] == "periodic":
        raise NotImplementedError("general-form 3D SOR currently supports periodic boundaries only in x")
    if _is_dataarray(H):
        return _invert_general_3d_labeled(
            H,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            iparams=iparams,
            icbc=icbc,
            coefficients=coefficients,
        )
    return _invert_general_3d_ndarray(
        H,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        iparams=iparams,
        icbc=icbc,
        coefficients=coefficients,
    )


def _invert_general_3d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any, Any],
) -> Any:
    if dims is None:
        if field.ndim < 3:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 3 dimensions")
        dims = field.dims[-3:]
    if len(dims) != 3 or not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray general-form 3D inversion requires three dimension names")
    mask_f, init_s, _zero = _mask_labeled_field(field, dims, iparams, bcs, icbc)
    dz, dy, dx = _spacing_for_labeled3d(mask_f, (dims[0], dims[1], dims[2]), spacing)
    coefs = tuple(_as_labeled_like(coef, mask_f) for coef in coefficients)
    solved = _solve_sor_general3d_labeled(
        init_s,
        *coefs,
        mask_f,
        dims=dims,
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def _invert_general_3d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any, Any],
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 3:
        raise ValueError("general-form 3D inversion requires at least a three-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 3, arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 3 or not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy general-form 3D inversion requires three integer axes")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if len(set(axes_tuple)) != 3:
        raise ValueError("inversion axes must be distinct")

    dz, dy, dx = _normalize_spacing3d(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-3, -2, -1))
    init = (
        np.zeros_like(moved)
        if icbc is None
        else np.moveaxis(np.asarray(icbc, dtype=np.float64), axes_tuple, (-3, -2, -1))
    )
    moved_coefficients = []
    for coef in coefficients:
        coef_arr = np.asarray(coef, dtype=np.float64)
        if coef_arr.ndim == arr.ndim and coef_arr.shape == arr.shape:
            coef_arr = np.moveaxis(coef_arr, axes_tuple, (-3, -2, -1))
        moved_coefficients.append(coef_arr)
    coef_arrays = tuple(
        _broadcast_ndarray_coefficient(coef, moved.shape) for coef in moved_coefficients
    )
    solved = _solve_sor_general3d_batched(
        init,
        *coef_arrays,
        moved,
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return np.moveaxis(solved, (-3, -2, -1), axes_tuple)


def _solve_sor_general2d_labeled(
    init_s: Any,
    acoef: Any,
    bcoef: Any,
    ccoef: Any,
    dcoef: Any,
    ecoef: Any,
    fcoef: Any,
    force: Any,
    *,
    dims: Sequence[str],
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> Any:
    import xarray as xr

    ydim, xdim = dims
    acoef, bcoef, ccoef, dcoef, ecoef, fcoef, force, init_s = xr.broadcast(
        acoef, bcoef, ccoef, dcoef, ecoef, fcoef, force, init_s
    )
    outer_dims = tuple(dim for dim in force.dims if dim not in dims)
    order = (*outer_dims, ydim, xdim)
    arrays = [item.transpose(*order).values for item in (init_s, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, force)]
    values = _solve_sor_general2d_batched(
        arrays[0],
        arrays[1],
        arrays[2],
        arrays[3],
        arrays[4],
        arrays[5],
        arrays[6],
        arrays[7],
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    template = force.transpose(*order)
    result = force.__class__(
        values,
        coords=template.coords,
        dims=template.dims,
        attrs=dict(force.attrs),
        name="inverted",
    )
    return result.transpose(*force.dims)


def _solve_sor_general2d_batched(
    init_s: np.ndarray,
    acoef: np.ndarray,
    bcoef: np.ndarray,
    ccoef: np.ndarray,
    dcoef: np.ndarray,
    ecoef: np.ndarray,
    fcoef: np.ndarray,
    force: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> np.ndarray:
    arrays = [np.asarray(item, dtype=np.float64) for item in (init_s, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, force)]
    shape = np.broadcast_shapes(*(item.shape for item in arrays))
    arrays = [np.broadcast_to(item, shape) for item in arrays]
    if len(shape) < 2:
        raise ValueError("general-form SOR expects the last two dimensions to be spatial")
    values = np.empty(shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, shape[-2:])
    if len(shape) == 2:
        solved, relerr, overflow, loops = fishpack.sor_general2d(
            arrays[0], arrays[1], arrays[2], arrays[3], arrays[4], arrays[5], arrays[6], arrays[7],
            dy, dx, bcs[0], bcs[1], optarg, _UNDEF, int(iparams["mxLoop"]), float(iparams["tolerance"])
        )
        if overflow:
            raise RuntimeError("Fortran SOR general2d overflowed")
        values[...] = solved
        return values
    for index in np.ndindex(shape[:-2]):
        solved, relerr, overflow, loops = fishpack.sor_general2d(
            arrays[0][index], arrays[1][index], arrays[2][index], arrays[3][index],
            arrays[4][index], arrays[5][index], arrays[6][index], arrays[7][index],
            dy, dx, bcs[0], bcs[1], optarg, _UNDEF, int(iparams["mxLoop"]), float(iparams["tolerance"])
        )
        if overflow:
            raise RuntimeError("Fortran SOR general2d overflowed")
        values[index] = solved
    return values


def _solve_sor_general3d_labeled(
    init_s: Any,
    acoef: Any,
    bcoef: Any,
    ccoef: Any,
    dcoef: Any,
    ecoef: Any,
    fcoef: Any,
    gcoef: Any,
    force: Any,
    *,
    dims: Sequence[str],
    dz: float,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> Any:
    import xarray as xr

    zdim, ydim, xdim = dims
    arrays = xr.broadcast(acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, force, init_s)
    acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, force, init_s = arrays
    outer_dims = tuple(dim for dim in force.dims if dim not in dims)
    order = (*outer_dims, zdim, ydim, xdim)
    transposed_force = force.transpose(*order)
    arrays_np = [
        item.transpose(*order).values
        for item in (init_s, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, force)
    ]
    values = _solve_sor_general3d_batched(
        arrays_np[0],
        arrays_np[1],
        arrays_np[2],
        arrays_np[3],
        arrays_np[4],
        arrays_np[5],
        arrays_np[6],
        arrays_np[7],
        arrays_np[8],
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    result = force.__class__(
        values,
        coords=transposed_force.coords,
        dims=transposed_force.dims,
        attrs=dict(force.attrs),
        name="inverted",
    )
    return result.transpose(*force.dims)


def _solve_sor_general3d_batched(
    init_s: np.ndarray,
    acoef: np.ndarray,
    bcoef: np.ndarray,
    ccoef: np.ndarray,
    dcoef: np.ndarray,
    ecoef: np.ndarray,
    fcoef: np.ndarray,
    gcoef: np.ndarray,
    force: np.ndarray,
    *,
    dz: float,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> np.ndarray:
    arrays = [
        np.asarray(item, dtype=np.float64)
        for item in (init_s, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, force)
    ]
    shape = np.broadcast_shapes(*(item.shape for item in arrays))
    arrays = [np.broadcast_to(item, shape) for item in arrays]
    if len(shape) < 3:
        raise ValueError("general-form 3D SOR expects the last three dimensions to be spatial")
    values = np.empty(shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, shape[-3:])
    if len(shape) == 3:
        solved, relerr, overflow, loops = fishpack.sor_general3d(
            arrays[0],
            arrays[1],
            arrays[2],
            arrays[3],
            arrays[4],
            arrays[5],
            arrays[6],
            arrays[7],
            arrays[8],
            dz,
            dy,
            dx,
            bcs[0],
            bcs[1],
            bcs[2],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        if overflow:
            raise RuntimeError("Fortran SOR general3d overflowed")
        values[...] = solved
        return values
    for index in np.ndindex(shape[:-3]):
        solved, relerr, overflow, loops = fishpack.sor_general3d(
            arrays[0][index],
            arrays[1][index],
            arrays[2][index],
            arrays[3][index],
            arrays[4][index],
            arrays[5][index],
            arrays[6][index],
            arrays[7][index],
            arrays[8][index],
            dz,
            dy,
            dx,
            bcs[0],
            bcs[1],
            bcs[2],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        if overflow:
            raise RuntimeError("Fortran SOR general3d overflowed")
        values[index] = solved
    return values


def _invert_biharmonic_2d(
    J: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any, Any, Any, Any],
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("biharmonic SOR currently supports Cartesian coordinates only")
    iparams = _merged_iparams(iParams, ndim=2, BCs=BCs)
    bcs = _normalize_sor_bcs(iparams["BCs"], 2)
    if _is_dataarray(J):
        return _invert_biharmonic_2d_labeled(
            J,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            iparams=iparams,
            icbc=icbc,
            coefficients=coefficients,
        )
    return _invert_biharmonic_2d_ndarray(
        J,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        iparams=iparams,
        icbc=icbc,
        coefficients=coefficients,
    )


def _invert_biharmonic_2d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any, Any, Any, Any],
) -> Any:
    if dims is None:
        if field.ndim < 2:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
        dims = field.dims[-2:]
    if len(dims) != 2 or not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray biharmonic inversion requires two dimension names")
    mask_f, init_s, _zero = _mask_labeled_field(field, dims, iparams, bcs, icbc)
    dy, dx = _spacing_for_labeled(mask_f, (dims[0], dims[1]), spacing)
    coefs = tuple(_as_labeled_like(coef, mask_f) for coef in coefficients)
    solved = _solve_sor_biharmonic2d_labeled(
        init_s,
        *coefs,
        mask_f,
        dims=dims,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return _restore_labeled_result(solved, mask_f, iparams, icbc)


def _invert_biharmonic_2d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, ...],
    spacing: Sequence[float] | None,
    iparams: dict[str, Any],
    icbc: Any,
    coefficients: tuple[Any, Any, Any, Any, Any, Any, Any, Any, Any],
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("biharmonic inversion requires at least a two-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 2 or not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy biharmonic inversion requires two integer axes")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if axes_tuple[0] == axes_tuple[1]:
        raise ValueError("inversion axes must be distinct")
    dy, dx = _normalize_spacing(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-2, -1))
    init = np.zeros_like(moved) if icbc is None else np.moveaxis(np.asarray(icbc, dtype=np.float64), axes_tuple, (-2, -1))
    moved_coefficients = []
    for coef in coefficients:
        coef_arr = np.asarray(coef, dtype=np.float64)
        if coef_arr.ndim == arr.ndim and coef_arr.shape == arr.shape:
            coef_arr = np.moveaxis(coef_arr, axes_tuple, (-2, -1))
        moved_coefficients.append(coef_arr)
    coef_arrays = tuple(_broadcast_ndarray_coefficient(coef, moved.shape) for coef in moved_coefficients)
    solved = _solve_sor_biharmonic2d_batched(
        init,
        *coef_arrays,
        moved,
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    return np.moveaxis(solved, (-2, -1), axes_tuple)


def _solve_sor_biharmonic2d_labeled(
    init_s: Any,
    acoef: Any,
    bcoef: Any,
    ccoef: Any,
    dcoef: Any,
    ecoef: Any,
    fcoef: Any,
    gcoef: Any,
    hcoef: Any,
    icoef: Any,
    force: Any,
    *,
    dims: Sequence[str],
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> Any:
    import xarray as xr

    ydim, xdim = dims
    arrays = xr.broadcast(acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, hcoef, icoef, force, init_s)
    acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, hcoef, icoef, force, init_s = arrays
    outer_dims = tuple(dim for dim in force.dims if dim not in dims)
    order = (*outer_dims, ydim, xdim)
    arrays_np = [
        item.transpose(*order).values
        for item in (init_s, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, hcoef, icoef, force)
    ]
    values = _solve_sor_biharmonic2d_batched(
        arrays_np[0],
        arrays_np[1],
        arrays_np[2],
        arrays_np[3],
        arrays_np[4],
        arrays_np[5],
        arrays_np[6],
        arrays_np[7],
        arrays_np[8],
        arrays_np[9],
        arrays_np[10],
        dy=dy,
        dx=dx,
        bcs=bcs,
        iparams=iparams,
    )
    template = force.transpose(*order)
    result = force.__class__(
        values,
        coords=template.coords,
        dims=template.dims,
        attrs=dict(force.attrs),
        name="inverted",
    )
    return result.transpose(*force.dims)


def _solve_sor_biharmonic2d_batched(
    init_s: np.ndarray,
    acoef: np.ndarray,
    bcoef: np.ndarray,
    ccoef: np.ndarray,
    dcoef: np.ndarray,
    ecoef: np.ndarray,
    fcoef: np.ndarray,
    gcoef: np.ndarray,
    hcoef: np.ndarray,
    icoef: np.ndarray,
    force: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> np.ndarray:
    arrays = [
        np.asarray(item, dtype=np.float64)
        for item in (init_s, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, hcoef, icoef, force)
    ]
    shape = np.broadcast_shapes(*(item.shape for item in arrays))
    arrays = [np.broadcast_to(item, shape) for item in arrays]
    if len(shape) < 2:
        raise ValueError("biharmonic SOR expects the last two dimensions to be spatial")
    values = np.empty(shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, shape[-2:])
    if len(shape) == 2:
        solved, relerr, overflow, loops = fishpack.sor_biharmonic2d(
            arrays[0], arrays[1], arrays[2], arrays[3], arrays[4], arrays[5],
            arrays[6], arrays[7], arrays[8], arrays[9], arrays[10],
            dy, dx, bcs[0], bcs[1], optarg, _UNDEF, int(iparams["mxLoop"]), float(iparams["tolerance"])
        )
        if overflow:
            raise RuntimeError("Fortran SOR biharmonic2d overflowed")
        values[...] = solved
        return values
    for index in np.ndindex(shape[:-2]):
        solved, relerr, overflow, loops = fishpack.sor_biharmonic2d(
            arrays[0][index], arrays[1][index], arrays[2][index], arrays[3][index],
            arrays[4][index], arrays[5][index], arrays[6][index], arrays[7][index],
            arrays[8][index], arrays[9][index], arrays[10][index],
            dy, dx, bcs[0], bcs[1], optarg, _UNDEF, int(iparams["mxLoop"]), float(iparams["tolerance"])
        )
        if overflow:
            raise RuntimeError("Fortran SOR biharmonic2d overflowed")
        values[index] = solved
    return values


def _broadcast_ndarray_coefficient(coef: Any, shape: tuple[int, ...]) -> np.ndarray:
    arr = np.asarray(coef, dtype=np.float64)
    if arr.ndim == 0:
        return np.full(shape, float(arr), dtype=np.float64)
    return np.broadcast_to(arr, shape).astype(np.float64, copy=False)


def _solve_sor1d_labeled(
    init_s: Any,
    acoef: Any,
    bcoef: Any,
    force: Any,
    *,
    dims: Sequence[str],
    dx: float,
    bcs: tuple[str, ...],
    iparams: dict[str, Any],
) -> Any:
    import xarray as xr

    (dim,) = dims
    acoef, bcoef, force, init_s = xr.broadcast(acoef, bcoef, force, init_s)
    outer_dims = tuple(item for item in force.dims if item != dim)
    order = (*outer_dims, dim)
    a_t = acoef.transpose(*order)
    b_t = bcoef.transpose(*order)
    f_t = force.transpose(*order)
    s_t = init_s.transpose(*order)
    values = np.empty(s_t.shape, dtype=np.float64)
    optarg = _sor_optarg(iparams, (s_t.shape[-1],))

    for index in np.ndindex(s_t.shape[:-1] or ()):
        key = index if s_t.ndim > 1 else (...,)
        if s_t.ndim == 1:
            solved, relerr, overflow, loops = fishpack.sor_standard1d(
                s_t.values,
                a_t.values,
                b_t.values,
                f_t.values,
                dx,
                bcs[0],
                optarg,
                _UNDEF,
                int(iparams["mxLoop"]),
                float(iparams["tolerance"]),
            )
            values[...] = solved
            if overflow:
                raise RuntimeError("Fortran SOR standard1d overflowed")
            break
        solved, relerr, overflow, loops = fishpack.sor_standard1d(
            s_t.values[key],
            a_t.values[key],
            b_t.values[key],
            f_t.values[key],
            dx,
            bcs[0],
            optarg,
            _UNDEF,
            int(iparams["mxLoop"]),
            float(iparams["tolerance"]),
        )
        values[key] = solved
        if overflow:
            raise RuntimeError("Fortran SOR standard1d overflowed")

    result = force.__class__(
        values,
        coords=s_t.coords,
        dims=s_t.dims,
        attrs=dict(force.attrs),
        name="inverted",
    )
    return result.transpose(*force.dims)


def _sor_optarg(iparams: dict[str, Any], shape: tuple[int, ...]) -> float:
    if iparams.get("optArg") is not None:
        return float(iparams["optArg"])
    if len(shape) == 1:
        epsilon = np.sin(np.pi / (2.0 * shape[0] + 2.0)) ** 2
    elif len(shape) == 2:
        epsilon = (
            np.sin(np.pi / (2.0 * shape[1] + 2.0)) ** 2
            + np.sin(np.pi / (2.0 * shape[0] + 2.0)) ** 2
        )
    elif len(shape) == 3:
        epsilon = (
            np.sin(np.pi / (2.0 * shape[2] + 2.0)) ** 2
            + np.sin(np.pi / (2.0 * shape[1] + 2.0)) ** 2
            + np.sin(np.pi / (2.0 * shape[0] + 2.0)) ** 2
        )
    else:
        raise ValueError("SOR optArg supports one, two, or three dimensions")
    return float(2.0 / (1.0 + np.sqrt((2.0 - epsilon) * epsilon)))


def _sor_spacing2d(field: Any, dims: tuple[str, str], coords: str, rearth: float) -> tuple[float, float]:
    ydim, xdim = dims
    dy = _uniform_coord_spacing(field.coords[ydim].values, ydim)
    dx = _uniform_coord_spacing(field.coords[xdim].values, xdim)
    if coords.lower() == "z-lat":
        dx = np.deg2rad(dx) * rearth
    return dy, dx


def _sor_spacing1d(field: Any, dim: str, coords: str, rearth: float) -> float:
    dx = _uniform_coord_spacing(field.coords[dim].values, dim)
    if coords.lower() == "lat":
        dx = np.deg2rad(dx) * rearth
    return dx


def _as_labeled_like(value: Any, template: Any) -> Any:
    if _is_dataarray(value):
        return value
    return (template * 0.0) + value


def _second_diff_swm(m0: Any, cos_h: Any, delx: float, dim: str) -> Any:
    shifted_forward = m0.shift({dim: -1})
    shifted_backward = m0.shift({dim: 1})
    flux_forward = (shifted_forward - m0) / cos_h.shift({dim: -1})
    flux_backward = (m0 - shifted_backward) / cos_h
    diff = (flux_forward - flux_backward) / (delx * delx)
    return diff.fillna(0.0)


def _invert_constant_helmholtz(
    F: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    helmholtz: float,
    raise_on_error: bool,
) -> Any:
    return _invert_constant_2d(
        F,
        dims=dims,
        coords=coords,
        iParams=iParams,
        BCs=BCs,
        spacing=spacing,
        coefficients=(1.0, 1.0),
        helmholtz=helmholtz,
        raise_on_error=raise_on_error,
    )


def _invert_constant_2d(
    F: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    coefficients: tuple[float, float],
    helmholtz: float,
    raise_on_error: bool,
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("only Cartesian constant-coefficient equations are supported")
    params = dict(iParams or {})
    bcs = _normalize_bcs(BCs if BCs is not None else params.get("BCs", None))
    if _is_dataarray(F):
        return _invert_constant_2d_labeled(
            F,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            coefficients=coefficients,
            helmholtz=helmholtz,
            raise_on_error=raise_on_error,
        )
    return _invert_constant_2d_ndarray(
        F,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        coefficients=coefficients,
        helmholtz=helmholtz,
        raise_on_error=raise_on_error,
    )


def _invert_constant_3d(
    F: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    coords: str,
    iParams: dict[str, Any] | None,
    BCs: Sequence[str] | None,
    spacing: Sequence[float] | None,
    coefficients: tuple[float, float, float],
    raise_on_error: bool,
) -> Any:
    if coords != "cartesian":
        raise NotImplementedError("only Cartesian constant-coefficient 3D equations are supported")
    params = dict(iParams or {})
    bcs = _normalize_bcs_n(BCs if BCs is not None else params.get("BCs", None), 3)
    if _is_dataarray(F):
        return _invert_constant_3d_labeled(
            F,
            dims=dims,
            bcs=bcs,
            spacing=spacing,
            coefficients=coefficients,
            raise_on_error=raise_on_error,
        )
    return _invert_constant_3d_ndarray(
        F,
        axes=dims,
        bcs=bcs,
        spacing=spacing,
        coefficients=coefficients,
        raise_on_error=raise_on_error,
    )


def _invert_poisson_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, str],
    spacing: Sequence[float] | None,
    raise_on_error: bool,
) -> Any:
    if dims is None:
        if field.ndim < 2:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
        dims = field.dims[-2:]
    if len(dims) != 2:
        raise ValueError("invert_Poisson requires exactly two inversion dimensions")
    if not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray inputs require dimension names in dims")
    missing = [dim for dim in dims if dim not in field.dims]
    if missing:
        raise ValueError(f"dims not present in input DataArray: {missing}")

    ydim, xdim = dims
    dy, dx = _spacing_for_labeled(field, (ydim, xdim), spacing)
    outer_dims = tuple(dim for dim in field.dims if dim not in dims)
    transposed = field.transpose(*outer_dims, ydim, xdim)

    values = _solve_poisson_batched(
        transposed.values,
        dy=dy,
        dx=dx,
        bcs=bcs,
        raise_on_error=raise_on_error,
    )

    result = field.__class__(
        values,
        coords=transposed.coords,
        dims=transposed.dims,
        attrs=dict(field.attrs),
        name="inverted" if field.name is None else f"{field.name}_inverted",
    )
    return result.transpose(*field.dims)


def _invert_poisson_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, str],
    spacing: Sequence[float] | None,
    raise_on_error: bool,
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("invert_Poisson requires at least a two-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 2:
            raise ValueError("invert_Poisson requires exactly two inversion axes")
        if not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy inputs require integer axes or axes=None")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if axes_tuple[0] == axes_tuple[1]:
        raise ValueError("inversion axes must be distinct")

    dy, dx = _normalize_spacing(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-2, -1))
    solved = _solve_poisson_batched(
        moved, dy=dy, dx=dx, bcs=bcs, raise_on_error=raise_on_error
    )
    return np.moveaxis(solved, (-2, -1), axes_tuple)


def _invert_constant_2d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, str],
    spacing: Sequence[float] | None,
    coefficients: tuple[float, float],
    helmholtz: float,
    raise_on_error: bool,
) -> Any:
    if dims is None:
        if field.ndim < 2:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
        dims = field.dims[-2:]
    if len(dims) != 2:
        raise ValueError("constant-coefficient inversion requires exactly two dimensions")
    if not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray inputs require dimension names in dims")
    missing = [dim for dim in dims if dim not in field.dims]
    if missing:
        raise ValueError(f"dims not present in input DataArray: {missing}")

    ydim, xdim = dims
    dy, dx = _spacing_for_labeled(field, (ydim, xdim), spacing)
    outer_dims = tuple(dim for dim in field.dims if dim not in dims)
    transposed = field.transpose(*outer_dims, ydim, xdim)
    values = _solve_constant_2d_batched(
        transposed.values,
        dy=dy,
        dx=dx,
        bcs=bcs,
        coefficients=coefficients,
        helmholtz=helmholtz,
        raise_on_error=raise_on_error,
    )
    result = field.__class__(
        values,
        coords=transposed.coords,
        dims=transposed.dims,
        attrs=dict(field.attrs),
        name="inverted" if field.name is None else f"{field.name}_inverted",
    )
    return result.transpose(*field.dims)


def _invert_constant_2d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, str],
    spacing: Sequence[float] | None,
    coefficients: tuple[float, float],
    helmholtz: float,
    raise_on_error: bool,
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("constant-coefficient inversion requires at least a two-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 2:
            raise ValueError("constant-coefficient inversion requires exactly two inversion axes")
        if not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy inputs require integer axes or axes=None")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if axes_tuple[0] == axes_tuple[1]:
        raise ValueError("inversion axes must be distinct")

    dy, dx = _normalize_spacing(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-2, -1))
    solved = _solve_constant_2d_batched(
        moved,
        dy=dy,
        dx=dx,
        bcs=bcs,
        coefficients=coefficients,
        helmholtz=helmholtz,
        raise_on_error=raise_on_error,
    )
    return np.moveaxis(solved, (-2, -1), axes_tuple)


def _invert_constant_3d_labeled(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, str, str],
    spacing: Sequence[float] | None,
    coefficients: tuple[float, float, float],
    raise_on_error: bool,
) -> Any:
    if dims is None:
        if field.ndim < 3:
            raise ValueError("dims must be supplied for xarray inputs with fewer than 3 dimensions")
        dims = field.dims[-3:]
    if len(dims) != 3:
        raise ValueError("constant-coefficient 3D inversion requires exactly three dimensions")
    if not all(isinstance(dim, str) for dim in dims):
        raise TypeError("xarray inputs require dimension names in dims")
    missing = [dim for dim in dims if dim not in field.dims]
    if missing:
        raise ValueError(f"dims not present in input DataArray: {missing}")

    zdim, ydim, xdim = dims
    dz, dy, dx = _spacing_for_labeled3d(field, (zdim, ydim, xdim), spacing)
    outer_dims = tuple(dim for dim in field.dims if dim not in dims)
    transposed = field.transpose(*outer_dims, zdim, ydim, xdim)
    values = _solve_constant_3d_batched(
        transposed.values,
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        coefficients=coefficients,
        raise_on_error=raise_on_error,
    )
    result = field.__class__(
        values,
        coords=transposed.coords,
        dims=transposed.dims,
        attrs=dict(field.attrs),
        name="inverted" if field.name is None else f"{field.name}_inverted",
    )
    return result.transpose(*field.dims)


def _invert_constant_3d_ndarray(
    field: Any,
    *,
    axes: Sequence[str] | Sequence[int] | None,
    bcs: tuple[str, str, str],
    spacing: Sequence[float] | None,
    coefficients: tuple[float, float, float],
    raise_on_error: bool,
) -> np.ndarray:
    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 3:
        raise ValueError("constant-coefficient 3D inversion requires at least a three-dimensional array")
    if axes is None:
        axes_tuple = (arr.ndim - 3, arr.ndim - 2, arr.ndim - 1)
    else:
        if len(axes) != 3:
            raise ValueError("constant-coefficient 3D inversion requires exactly three inversion axes")
        if not all(isinstance(axis, int) for axis in axes):
            raise TypeError("NumPy inputs require integer axes or axes=None")
        axes_tuple = tuple(axis % arr.ndim for axis in axes)
    if len(set(axes_tuple)) != 3:
        raise ValueError("inversion axes must be distinct")

    dz, dy, dx = _normalize_spacing3d(spacing)
    moved = np.moveaxis(arr, axes_tuple, (-3, -2, -1))
    solved = _solve_constant_3d_batched(
        moved,
        dz=dz,
        dy=dy,
        dx=dx,
        bcs=bcs,
        coefficients=coefficients,
        raise_on_error=raise_on_error,
    )
    return np.moveaxis(solved, (-3, -2, -1), axes_tuple)


def _solve_poisson_batched(
    values: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, str],
    raise_on_error: bool,
) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("Poisson solve expects the last two dimensions to be spatial")
    return _solve_helmholtz_batched(
        arr,
        dy=dy,
        dx=dx,
        bcs=bcs,
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def _solve_helmholtz_batched(
    values: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, str],
    helmholtz: float,
    raise_on_error: bool,
) -> np.ndarray:
    return _solve_constant_2d_batched(
        values,
        dy=dy,
        dx=dx,
        bcs=bcs,
        coefficients=(1.0, 1.0),
        helmholtz=helmholtz,
        raise_on_error=raise_on_error,
    )


def _solve_constant_2d_batched(
    values: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, str],
    coefficients: tuple[float, float],
    helmholtz: float,
    raise_on_error: bool,
) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("2D constant-coefficient solve expects the last two dimensions to be spatial")
    m, n = arr.shape[-2:]
    if m <= 2 or n <= 2:
        raise ValueError("Fishpack genbun requires both spatial dimensions to be > 2")

    result = np.empty(arr.shape, dtype=np.float64, order="C")
    if arr.ndim == 2:
        result[...] = _solve_constant_2d(
            arr,
            dy=dy,
            dx=dx,
            bcs=bcs,
            coefficients=coefficients,
            helmholtz=helmholtz,
            raise_on_error=raise_on_error,
        )
    else:
        for index in np.ndindex(arr.shape[:-2]):
            result[index] = _solve_constant_2d(
                arr[index],
                dy=dy,
                dx=dx,
                bcs=bcs,
                coefficients=coefficients,
                helmholtz=helmholtz,
                raise_on_error=raise_on_error,
            )
    return result


def _solve_poisson_2d(
    force: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, str],
    raise_on_error: bool,
) -> np.ndarray:
    return _solve_helmholtz_2d(
        force,
        dy=dy,
        dx=dx,
        bcs=bcs,
        helmholtz=0.0,
        raise_on_error=raise_on_error,
    )


def _solve_helmholtz_2d(
    force: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, str],
    helmholtz: float,
    raise_on_error: bool,
) -> np.ndarray:
    return _solve_constant_2d(
        force,
        dy=dy,
        dx=dx,
        bcs=bcs,
        coefficients=(1.0, 1.0),
        helmholtz=helmholtz,
        raise_on_error=raise_on_error,
    )


def _solve_constant_2d(
    force: np.ndarray,
    *,
    dy: float,
    dx: float,
    bcs: tuple[str, str],
    coefficients: tuple[float, float],
    helmholtz: float,
    raise_on_error: bool,
) -> np.ndarray:
    alpha_y, alpha_x = (float(item) for item in coefficients)
    if not np.isfinite(alpha_y) or not np.isfinite(alpha_x):
        raise ValueError("2D coefficients must be finite")
    if alpha_x == 0.0 or alpha_y == 0.0:
        raise ValueError("2D coefficients must be non-zero")
    if alpha_y / alpha_x <= 0.0:
        raise NotImplementedError(
            "Fishpack genbun path requires elliptic 2D coefficients with the same sign"
        )

    m, n = force.shape
    ratio = (alpha_y / alpha_x) * (dx / dy) ** 2
    a = np.full(m, ratio, dtype=np.float64)
    b = np.full(m, -2.0 * ratio + float(helmholtz) * dx * dx / alpha_x, dtype=np.float64)
    c = np.full(m, ratio, dtype=np.float64)

    if bcs[0] == "periodic":
        mperod = 0
    else:
        mperod = 1
        a[0] = 0.0
        c[-1] = 0.0

    nperod = 0 if bcs[1] == "periodic" else 1
    rhs = np.array(force, dtype=np.float64, order="F", copy=True) * (dx * dx / alpha_x)
    solution, ierror = fishpack.genbun(nperod, n, mperod, m, a, b, c, rhs)
    if raise_on_error and ierror != 0:
        raise RuntimeError(f"Fishpack genbun failed with ierror={ierror}")
    return np.asarray(solution)


def _solve_constant_3d_batched(
    values: np.ndarray,
    *,
    dz: float,
    dy: float,
    dx: float,
    bcs: tuple[str, str, str],
    coefficients: tuple[float, float, float],
    raise_on_error: bool,
) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim < 3:
        raise ValueError("3D solve expects the last three dimensions to be spatial")
    if min(arr.shape[-3:]) <= 2:
        raise ValueError("Fishpack pois3d requires all spatial dimensions to be > 2")

    result = np.empty(arr.shape, dtype=np.float64, order="C")
    if arr.ndim == 3:
        result[...] = _solve_constant_3d(
            arr,
            dz=dz,
            dy=dy,
            dx=dx,
            bcs=bcs,
            coefficients=coefficients,
            raise_on_error=raise_on_error,
        )
    else:
        for index in np.ndindex(arr.shape[:-3]):
            result[index] = _solve_constant_3d(
                arr[index],
                dz=dz,
                dy=dy,
                dx=dx,
                bcs=bcs,
                coefficients=coefficients,
                raise_on_error=raise_on_error,
            )
    return result


def _solve_constant_3d(
    force_zyx: np.ndarray,
    *,
    dz: float,
    dy: float,
    dx: float,
    bcs: tuple[str, str, str],
    coefficients: tuple[float, float, float],
    raise_on_error: bool,
) -> np.ndarray:
    alpha_z, alpha_y, alpha_x = (float(item) for item in coefficients)
    nz, ny, nx = force_zyx.shape
    f_xyz = np.asfortranarray(np.transpose(force_zyx, (2, 1, 0)))

    xperod = 0 if bcs[2] == "periodic" else 1
    yperod = 0 if bcs[1] == "periodic" else 1
    zperod = 0 if bcs[0] == "periodic" else 1
    c1 = alpha_x / (dx * dx)
    c2 = alpha_y / (dy * dy)
    zcoef = alpha_z / (dz * dz)
    a = np.full(nz, zcoef, dtype=np.float64)
    b = np.full(nz, -2.0 * zcoef, dtype=np.float64)
    c = np.full(nz, zcoef, dtype=np.float64)
    if bcs[0] != "periodic":
        a[0] = 0.0
        c[-1] = 0.0

    solution_xyz, ierror = fishpack.pois3d(
        xperod, nx, c1, yperod, ny, c2, zperod, nz, a, b, c, f_xyz
    )
    if raise_on_error and ierror != 0:
        raise RuntimeError(f"Fishpack pois3d failed with ierror={ierror}")
    return np.asarray(solution_xyz).transpose(2, 1, 0)


def _normalize_bcs(bcs: Sequence[str] | None) -> tuple[str, str]:
    return _normalize_bcs_n(bcs, 2)  # type: ignore[return-value]


def _normalize_bcs_n(bcs: Sequence[str] | None, ndim: int) -> tuple[str, ...]:
    if bcs is None:
        bcs = tuple("fixed" for _ in range(ndim))
    if len(bcs) != ndim:
        raise ValueError(f"BCs must contain {ndim} entries for the inversion dimensions")
    normalized = tuple(str(bc).lower() for bc in bcs)
    unsupported = [bc for bc in normalized if bc not in _SUPPORTED_BCS]
    if unsupported:
        raise NotImplementedError(
            "Only 'fixed' and 'periodic' boundary conditions are supported"
        )
    return normalized


def _normalize_sor_bcs(bcs: Sequence[str] | None, ndim: int) -> tuple[str, ...]:
    if bcs is None:
        bcs = tuple("fixed" for _ in range(ndim))
    if len(bcs) != ndim:
        raise ValueError(f"BCs must contain {ndim} entries for the inversion dimensions")
    normalized = tuple(str(bc).lower() for bc in bcs)
    unsupported = [bc for bc in normalized if bc not in _SUPPORTED_SOR_BCS]
    if unsupported:
        raise NotImplementedError(
            "SOR-backed inversions support 'fixed', 'extend', and 'periodic' boundary conditions"
        )
    return normalized


def _merged_mparams(params: dict[str, Any] | None) -> dict[str, Any]:
    merged = dict(_DEFAULT_MPARAMS)
    if params is not None:
        merged.update(params)
    return merged


def _scalar_param(params: dict[str, Any], name: str) -> float:
    value = params[name]
    array = np.asarray(value)
    if array.size != 1:
        raise NotImplementedError(f"{name} must be a scalar for the current Fishpack-backed subset")
    scalar = float(array.reshape(()))
    if not np.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    return scalar


def _cartesian_coriolis_forcing(
    field: Any,
    dims: Sequence[str] | Sequence[int] | None,
    spacing: Sequence[float] | None,
    params: dict[str, Any],
    *,
    sign: float,
    offset: float,
) -> Any:
    coriolis = _cartesian_coriolis_field(field, dims, spacing, params)
    if _is_dataarray(field):
        return (field * 0.0) + offset + sign * coriolis
    return offset + sign * np.asarray(coriolis, dtype=np.float64)


def _cartesian_coriolis_field(
    field: Any,
    dims: Sequence[str] | Sequence[int] | None,
    spacing: Sequence[float] | None,
    params: dict[str, Any],
) -> Any:
    f0 = float(params["f0"])
    beta = float(params["beta"])
    if beta == 0.0:
        return _like_field(field, f0)
    if _is_dataarray(field):
        if dims is None:
            if field.ndim < 2:
                raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
            dims = field.dims[-2:]
        if len(dims) != 2 or not all(isinstance(dim, str) for dim in dims):
            raise TypeError("xarray beta-plane forcing requires two dimension names")
        ydim = dims[0]
        y = field.coords[ydim]
        return f0 + beta * y

    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("beta-plane forcing requires at least a two-dimensional array")
    if dims is None:
        yaxis = arr.ndim - 2
    else:
        if len(dims) != 2 or not all(isinstance(axis, int) for axis in dims):
            raise TypeError("NumPy beta-plane forcing requires two integer axes")
        yaxis = int(dims[0]) % arr.ndim
    dy, _ = _normalize_spacing(spacing)
    y = np.arange(arr.shape[yaxis], dtype=np.float64) * dy
    shape = [1] * arr.ndim
    shape[yaxis] = arr.shape[yaxis]
    coriolis = f0 + beta * y.reshape(shape)
    return np.broadcast_to(coriolis, arr.shape)


def _cartesian_geostrophic_coefficients(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    spacing: Sequence[float] | None,
    f0: float,
    beta: float,
) -> tuple[Any, Any, Any]:
    if _is_dataarray(field):
        if dims is None:
            if field.ndim < 2:
                raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
            dims = field.dims[-2:]
        if len(dims) != 2 or not all(isinstance(dim, str) for dim in dims):
            raise TypeError("xarray geostrophic beta-plane inversion requires two dimension names")
        y = field.coords[dims[0]]
        f_center = f0 + beta * y
        f_half = f0 + beta * ((y + y.shift({dims[0]: 1})) / 2.0).fillna(y)
        return f_half, 0.0, f_center

    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("geostrophic beta-plane inversion requires at least a two-dimensional array")
    if dims is None:
        yaxis = arr.ndim - 2
    else:
        if len(dims) != 2 or not all(isinstance(axis, int) for axis in dims):
            raise TypeError("NumPy geostrophic beta-plane inversion requires two integer axes")
        yaxis = int(dims[0]) % arr.ndim
    dy, _ = _normalize_spacing(spacing)
    y = np.arange(arr.shape[yaxis], dtype=np.float64) * dy
    y_half = y.copy()
    if y_half.size > 1:
        y_half[1:] = 0.5 * (y[1:] + y[:-1])
    shape = [1] * arr.ndim
    shape[yaxis] = arr.shape[yaxis]
    f_center = f0 + beta * y.reshape(shape)
    f_half = f0 + beta * y_half.reshape(shape)
    return (
        np.broadcast_to(f_half, arr.shape),
        0.0,
        np.broadcast_to(f_center, arr.shape),
    )


def _cartesian_beta_general_coefficients(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    spacing: Sequence[float] | None,
    epsilon: float,
    f0: float,
    beta: float,
    scale: float,
    helmholtz: float,
) -> tuple[Any, Any, Any, Any, Any, Any]:
    if _is_dataarray(field):
        if dims is None:
            if field.ndim < 2:
                raise ValueError("dims must be supplied for xarray inputs with fewer than 2 dimensions")
            dims = field.dims[-2:]
        if len(dims) != 2 or not all(isinstance(dim, str) for dim in dims):
            raise TypeError("xarray beta-plane inversion requires two dimension names")
        y = field.coords[dims[0]]
        f = f0 + beta * y
        denom = epsilon * epsilon + f * f
        c1 = epsilon / denom
        c2_dy = beta * (epsilon * epsilon - f * f) / (denom * denom)
        c1_dy = -2.0 * epsilon * f * beta / (denom * denom)
        return (
            scale * c1,
            0.0,
            scale * c1,
            scale * c1_dy,
            -scale * c2_dy,
            helmholtz,
        )

    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 2:
        raise ValueError("beta-plane inversion requires at least a two-dimensional array")
    if dims is None:
        yaxis = arr.ndim - 2
    else:
        if len(dims) != 2 or not all(isinstance(axis, int) for axis in dims):
            raise TypeError("NumPy beta-plane inversion requires two integer axes")
        yaxis = int(dims[0]) % arr.ndim
    dy, _ = _normalize_spacing(spacing)
    y = np.arange(arr.shape[yaxis], dtype=np.float64) * dy
    shape = [1] * arr.ndim
    shape[yaxis] = arr.shape[yaxis]
    f = f0 + beta * y.reshape(shape)
    denom = epsilon * epsilon + f * f
    c1 = epsilon / denom
    c2_dy = beta * (epsilon * epsilon - f * f) / (denom * denom)
    c1_dy = -2.0 * epsilon * f * beta / (denom * denom)
    return (
        np.broadcast_to(scale * c1, arr.shape),
        0.0,
        np.broadcast_to(scale * c1, arr.shape),
        np.broadcast_to(scale * c1_dy, arr.shape),
        np.broadcast_to(-scale * c2_dy, arr.shape),
        helmholtz,
    )


def _cartesian_omega_coefficients(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    spacing: Sequence[float] | None,
    f0: float,
    beta: float,
    n2: float,
) -> tuple[Any, Any, Any]:
    if _is_dataarray(field):
        if dims is None:
            if field.ndim < 3:
                raise ValueError("dims must be supplied for xarray inputs with fewer than 3 dimensions")
            dims = field.dims[-3:]
        if len(dims) != 3 or not all(isinstance(dim, str) for dim in dims):
            raise TypeError("xarray omega beta-plane inversion requires three dimension names")
        y = field.coords[dims[1]]
        f = f0 + beta * y
        return f * f, n2, n2

    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 3:
        raise ValueError("omega beta-plane inversion requires at least a three-dimensional array")
    if dims is None:
        yaxis = arr.ndim - 2
    else:
        if len(dims) != 3 or not all(isinstance(axis, int) for axis in dims):
            raise TypeError("NumPy omega beta-plane inversion requires three integer axes")
        yaxis = int(dims[1]) % arr.ndim
    _, dy, _ = _normalize_spacing3d(spacing)
    y = np.arange(arr.shape[yaxis], dtype=np.float64) * dy
    shape = [1] * arr.ndim
    shape[yaxis] = arr.shape[yaxis]
    acoef = (f0 + beta * y.reshape(shape)) ** 2
    return (
        np.broadcast_to(acoef, arr.shape),
        n2,
        n2,
    )


def _cartesian_3d_ocean_coefficients(
    field: Any,
    *,
    dims: Sequence[str] | Sequence[int] | None,
    spacing: Sequence[float] | None,
    epsilon: float,
    f0: float,
    beta: float,
    n2: float,
    buoyancy_damping: float,
) -> tuple[Any, Any, Any, Any, Any, Any, Any]:
    vertical = buoyancy_damping / n2
    if _is_dataarray(field):
        if dims is None:
            if field.ndim < 3:
                raise ValueError("dims must be supplied for xarray inputs with fewer than 3 dimensions")
            dims = field.dims[-3:]
        if len(dims) != 3 or not all(isinstance(dim, str) for dim in dims):
            raise TypeError("xarray 3DOcean beta-plane inversion requires three dimension names")
        y = field.coords[dims[1]]
        f = f0 + beta * y
        denom = epsilon * epsilon + f * f
        c1 = epsilon / denom
        c1_dy = -2.0 * epsilon * f * beta / (denom * denom)
        c2_dy = beta * (epsilon * epsilon - f * f) / (denom * denom)
        return vertical, c1, c1, 0.0, c1_dy, -c2_dy, 0.0

    arr = np.asarray(field, dtype=np.float64)
    if arr.ndim < 3:
        raise ValueError("3DOcean beta-plane inversion requires at least a three-dimensional array")
    if dims is None:
        yaxis = arr.ndim - 2
    else:
        if len(dims) != 3 or not all(isinstance(axis, int) for axis in dims):
            raise TypeError("NumPy 3DOcean beta-plane inversion requires three integer axes")
        yaxis = int(dims[1]) % arr.ndim
    _, dy, _ = _normalize_spacing3d(spacing)
    y = np.arange(arr.shape[yaxis], dtype=np.float64) * dy
    shape = [1] * arr.ndim
    shape[yaxis] = arr.shape[yaxis]
    f = f0 + beta * y.reshape(shape)
    denom = epsilon * epsilon + f * f
    c1 = epsilon / denom
    c1_dy = -2.0 * epsilon * f * beta / (denom * denom)
    c2_dy = beta * (epsilon * epsilon - f * f) / (denom * denom)
    return (
        vertical,
        np.broadcast_to(c1, arr.shape),
        np.broadcast_to(c1, arr.shape),
        0.0,
        np.broadcast_to(c1_dy, arr.shape),
        np.broadcast_to(-c2_dy, arr.shape),
        0.0,
    )


def _like_field(field: Any, value: float) -> Any:
    if _is_dataarray(field):
        return (field * 0.0) + value
    return np.zeros_like(np.asarray(field, dtype=np.float64)) + value


def _normalize_spacing(spacing: Sequence[float] | None) -> tuple[float, float]:
    if spacing is None:
        return 1.0, 1.0
    if len(spacing) != 2:
        raise ValueError("spacing must contain dy and dx")
    dy, dx = (float(spacing[0]), float(spacing[1]))
    if dy <= 0.0 or dx <= 0.0:
        raise ValueError("spacing values must be positive")
    return dy, dx


def _normalize_spacing3d(spacing: Sequence[float] | None) -> tuple[float, float, float]:
    if spacing is None:
        return 1.0, 1.0, 1.0
    if len(spacing) != 3:
        raise ValueError("spacing must contain dz, dy, and dx")
    dz, dy, dx = (float(spacing[0]), float(spacing[1]), float(spacing[2]))
    if dz <= 0.0 or dy <= 0.0 or dx <= 0.0:
        raise ValueError("spacing values must be positive")
    return dz, dy, dx


def _spacing_for_labeled(
    field: Any, dims: tuple[str, str], spacing: Sequence[float] | None
) -> tuple[float, float]:
    if spacing is not None:
        return _normalize_spacing(spacing)
    return tuple(_uniform_coord_spacing(field.coords[dim].values, dim) for dim in dims)  # type: ignore[return-value]


def _spacing_for_labeled3d(
    field: Any, dims: tuple[str, str, str], spacing: Sequence[float] | None
) -> tuple[float, float, float]:
    if spacing is not None:
        return _normalize_spacing3d(spacing)
    return tuple(_uniform_coord_spacing(field.coords[dim].values, dim) for dim in dims)  # type: ignore[return-value]


def _uniform_coord_spacing(coord: Any, dim: str) -> float:
    values = np.asarray(coord, dtype=np.float64)
    if values.ndim != 1 or values.size < 2:
        raise ValueError(f"coordinate {dim!r} must be one-dimensional with at least two points")
    deltas = np.diff(values)
    delta = float(deltas[0])
    if not np.allclose(deltas, delta):
        raise ValueError(f"coordinate {dim!r} must be uniformly spaced")
    if delta == 0.0:
        raise ValueError(f"coordinate {dim!r} has zero spacing")
    return abs(delta)


def _is_dataarray(obj: Any) -> bool:
    return obj.__class__.__module__.startswith("xarray.") and hasattr(obj, "dims")
