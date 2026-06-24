"""Verify the modern Fortran PyFishPack backend.

This script is intentionally executable without pytest so it can be used as a
build acceptance check in the conda environment that owns the Fortran toolchain.
It validates the compiled backend, UCAR Fishpack example errors for the exposed
entry points, a high-level Poisson inversion, size stability, and optionally a
Poisson speed comparison against a local xinvert checkout.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_XINVERT_PATH = Path(r"I:\test_xinvert\xinvert-master")


def _import_pyfishpack() -> Any:
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    import PyFishPack as pyfishpack

    return pyfishpack


def _assert_close(name: str, actual: float, expected: float, tol: float) -> dict[str, Any]:
    delta = abs(actual - expected)
    if delta > tol:
        raise AssertionError(
            f"{name}: got {actual:.17g}, expected {expected:.17g}, delta {delta:.3g}"
        )
    return {"actual": actual, "expected": expected, "abs_delta": delta, "tol": tol}


def check_backend_info(pyfishpack: Any) -> dict[str, Any]:
    info = dict(pyfishpack.fishpack.backend_info())
    expected = {
        "backend": "modern-fortran",
        "abi": "fortran-bind-c",
        "wrapper": "numpy-c-api",
    }
    if info != expected:
        raise AssertionError(f"unexpected backend info: {info!r}")
    return info


def check_backend_methods(pyfishpack: Any) -> dict[str, Any]:
    expected = {
        "backend_info",
        "genbun",
        "poistg",
        "pois3d",
        "hwscrt",
        "hstcrt",
        "hw3crt",
        "hwsplr",
        "hwscyl",
        "hstplr",
        "hstcyl",
        "hwsssp",
        "hstssp",
        "hwscsp",
        "hstcsp",
        "sor_standard1d",
        "sor_standard2d",
        "sor_standard3d",
        "sor_general2d",
        "sor_general3d",
        "sor_biharmonic2d",
        "rfftf",
        "rfftb",
        "sint",
        "cost",
        "sinqf",
        "sinqb",
        "cosqf",
        "cosqb",
        "cfftf",
        "cfftb",
    }
    methods = {name for name in dir(pyfishpack.fishpack) if not name.startswith("_")}
    missing = sorted(expected - methods)
    if missing:
        raise AssertionError(f"missing backend methods: {missing}")
    return {"expected": sorted(expected), "available": sorted(methods)}


def check_sor_standard1d_backend(pyfishpack: Any) -> dict[str, Any]:
    n = 40
    s = np.zeros(n, dtype=np.float64)
    a = np.ones(n, dtype=np.float64)
    b = np.zeros(n, dtype=np.float64)
    f = np.ones(n, dtype=np.float64)
    solution, relerr, overflow, loops = pyfishpack.fishpack.sor_standard1d(
        s, a, b, f, 1.0, "fixed", 1.7, -9.99e8, 5000, 1e-12
    )
    if overflow:
        raise AssertionError("sor_standard1d overflowed")
    residual = (solution[2:] - 2.0 * solution[1:-1] + solution[:-2]) - f[1:-1]
    max_residual = float(np.max(np.abs(residual)))
    checked = _assert_close("sor_standard1d residual", max_residual, 0.0, 1e-9)
    checked.update({"relerr": float(relerr), "loops": int(loops)})
    return checked


def check_sor_standard2d_backend(pyfishpack: Any) -> dict[str, Any]:
    ny, nx = 30, 32
    s = np.zeros((ny, nx), dtype=np.float64)
    a = np.ones((ny, nx), dtype=np.float64)
    b = np.zeros((ny, nx), dtype=np.float64)
    c = np.ones((ny, nx), dtype=np.float64)
    f = np.ones((ny, nx), dtype=np.float64)
    solution, relerr, overflow, loops = pyfishpack.fishpack.sor_standard2d(
        s, a, b, c, f, 1.0, 1.0, "fixed", "fixed", 1.7, -9.99e8, 5000, 1e-12
    )
    if overflow:
        raise AssertionError("sor_standard2d overflowed")
    residual = (
        solution[2:, 1:-1]
        - 2.0 * solution[1:-1, 1:-1]
        + solution[:-2, 1:-1]
        + solution[1:-1, 2:]
        - 2.0 * solution[1:-1, 1:-1]
        + solution[1:-1, :-2]
        - f[1:-1, 1:-1]
    )
    max_residual = float(np.max(np.abs(residual)))
    checked = _assert_close("sor_standard2d residual", max_residual, 0.0, 1e-9)
    checked.update({"relerr": float(relerr), "loops": int(loops)})
    return checked


def check_sor_general2d_backend(pyfishpack: Any) -> dict[str, Any]:
    ny, nx = 32, 34
    dy, dx = 1.1, 0.9
    alpha = -0.25
    beta = 0.015
    y = np.arange(ny, dtype=np.float64)[:, None]
    x = np.arange(nx, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * y / (ny - 1.0)) * np.sin(math.pi * x / (nx - 1.0))
    force = np.zeros_like(exact)
    force[1:-1, 1:-1] = (
        alpha
        * (
            (exact[2:, 1:-1] - 2.0 * exact[1:-1, 1:-1] + exact[:-2, 1:-1])
            / (dy * dy)
            + (exact[1:-1, 2:] - 2.0 * exact[1:-1, 1:-1] + exact[1:-1, :-2])
            / (dx * dx)
        )
        - beta * (exact[1:-1, 2:] - exact[1:-1, :-2]) / (2.0 * dx)
    )
    zeros = np.zeros_like(exact)
    solution, relerr, overflow, loops = pyfishpack.fishpack.sor_general2d(
        zeros,
        np.full_like(exact, alpha),
        zeros,
        np.full_like(exact, alpha),
        zeros,
        np.full_like(exact, -beta),
        zeros,
        force,
        dy,
        dx,
        "fixed",
        "fixed",
        1.7,
        -9.99e8,
        20000,
        1e-11,
    )
    if overflow:
        raise AssertionError("sor_general2d overflowed")
    max_error = float(np.max(np.abs(solution - exact)))
    checked = _assert_close("sor_general2d manufactured solution", max_error, 0.0, 1e-8)
    checked.update({"relerr": float(relerr), "loops": int(loops)})
    return checked


def check_sor_standard3d_backend(pyfishpack: Any) -> dict[str, Any]:
    nz, ny, nx = 12, 13, 14
    dz, dy, dx = 1.2, 0.9, 1.1
    alpha_z, alpha_y, alpha_x = 0.7, 1.2, 0.9
    k = np.arange(nz, dtype=np.float64)[:, None, None]
    j = np.arange(ny, dtype=np.float64)[None, :, None]
    i = np.arange(nx, dtype=np.float64)[None, None, :]
    exact = (
        np.sin(math.pi * k / (nz - 1.0))
        * np.sin(math.pi * j / (ny - 1.0))
        * np.sin(math.pi * i / (nx - 1.0))
    )
    force = _discrete_constant_3d_fixed(exact, dz, dy, dx, alpha_z, alpha_y, alpha_x)
    solution, relerr, overflow, loops = pyfishpack.fishpack.sor_standard3d(
        np.zeros_like(exact),
        np.full_like(exact, alpha_z),
        np.full_like(exact, alpha_y),
        np.full_like(exact, alpha_x),
        force,
        dz,
        dy,
        dx,
        "fixed",
        "fixed",
        "fixed",
        1.6,
        -9.99e8,
        50000,
        1e-12,
    )
    if overflow:
        raise AssertionError("sor_standard3d overflowed")
    max_error = float(np.max(np.abs(solution - exact)))
    checked = _assert_close("sor_standard3d manufactured solution", max_error, 0.0, 1e-8)
    checked.update({"relerr": float(relerr), "loops": int(loops)})
    return checked


def check_sor_general3d_backend(pyfishpack: Any) -> dict[str, Any]:
    nz, ny, nx = 12, 13, 14
    dz, dy, dx = 1.2, 0.9, 1.1
    alpha_z, alpha_y, alpha_x = 0.7, 1.2, 0.9
    dzcoef, dycoef, dxcoef = 0.05, -0.03, 0.02
    helmholtz = -0.11
    k = np.arange(nz, dtype=np.float64)[:, None, None]
    j = np.arange(ny, dtype=np.float64)[None, :, None]
    i = np.arange(nx, dtype=np.float64)[None, None, :]
    exact = (
        np.sin(math.pi * k / (nz - 1.0))
        * np.sin(math.pi * j / (ny - 1.0))
        * np.sin(math.pi * i / (nx - 1.0))
    )
    force = np.zeros_like(exact)
    force[1:-1, 1:-1, 1:-1] = (
        alpha_z
        * (exact[2:, 1:-1, 1:-1] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[:-2, 1:-1, 1:-1])
        / (dz * dz)
        + alpha_y
        * (exact[1:-1, 2:, 1:-1] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[1:-1, :-2, 1:-1])
        / (dy * dy)
        + alpha_x
        * (exact[1:-1, 1:-1, 2:] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[1:-1, 1:-1, :-2])
        / (dx * dx)
        + dzcoef * (exact[2:, 1:-1, 1:-1] - exact[:-2, 1:-1, 1:-1]) / (2.0 * dz)
        + dycoef * (exact[1:-1, 2:, 1:-1] - exact[1:-1, :-2, 1:-1]) / (2.0 * dy)
        + dxcoef * (exact[1:-1, 1:-1, 2:] - exact[1:-1, 1:-1, :-2]) / (2.0 * dx)
        + helmholtz * exact[1:-1, 1:-1, 1:-1]
    )
    solution, relerr, overflow, loops = pyfishpack.fishpack.sor_general3d(
        np.zeros_like(exact),
        np.full_like(exact, alpha_z),
        np.full_like(exact, alpha_y),
        np.full_like(exact, alpha_x),
        np.full_like(exact, dzcoef),
        np.full_like(exact, dycoef),
        np.full_like(exact, dxcoef),
        np.full_like(exact, helmholtz),
        force,
        dz,
        dy,
        dx,
        "fixed",
        "fixed",
        "fixed",
        1.6,
        -9.99e8,
        50000,
        1e-12,
    )
    if overflow:
        raise AssertionError("sor_general3d overflowed")
    max_error = float(np.max(np.abs(solution - exact)))
    checked = _assert_close("sor_general3d manufactured solution", max_error, 0.0, 1e-8)
    checked.update({"relerr": float(relerr), "loops": int(loops)})
    return checked


def _discrete_biharmonic2d_fixed(
    exact: np.ndarray,
    dy: float,
    dx: float,
    acoef: float,
    bcoef: float,
    ccoef: float,
    dcoef: float,
    ecoef: float,
    fcoef: float,
    gcoef: float,
    hcoef: float,
    icoef: float,
) -> np.ndarray:
    force = np.zeros_like(exact)
    force[2:-2, 2:-2] = (
        acoef
        * (exact[4:, 2:-2] - 4.0 * exact[3:-1, 2:-2] + 6.0 * exact[2:-2, 2:-2] - 4.0 * exact[1:-3, 2:-2] + exact[:-4, 2:-2])
        / (dy**4)
        + bcoef
        * (
            exact[4:, 4:] - 2.0 * exact[4:, 2:-2] + exact[4:, :-4]
            - 2.0 * exact[2:-2, 4:] + 4.0 * exact[2:-2, 2:-2] - 2.0 * exact[2:-2, :-4]
            + exact[:-4, 4:] - 2.0 * exact[:-4, 2:-2] + exact[:-4, :-4]
        )
        / (16.0 * dy * dy * dx * dx)
        + ccoef
        * (exact[2:-2, 4:] - 4.0 * exact[2:-2, 3:-1] + 6.0 * exact[2:-2, 2:-2] - 4.0 * exact[2:-2, 1:-3] + exact[2:-2, :-4])
        / (dx**4)
        + dcoef
        * (exact[3:-1, 2:-2] - 2.0 * exact[2:-2, 2:-2] + exact[1:-3, 2:-2])
        / (dy * dy)
        + ecoef
        * (exact[3:-1, 3:-1] - exact[1:-3, 3:-1] - exact[3:-1, 1:-3] + exact[1:-3, 1:-3])
        / (4.0 * dy * dx)
        + fcoef
        * (exact[2:-2, 3:-1] - 2.0 * exact[2:-2, 2:-2] + exact[2:-2, 1:-3])
        / (dx * dx)
        + gcoef * (exact[3:-1, 2:-2] - exact[1:-3, 2:-2]) / (2.0 * dy)
        + hcoef * (exact[2:-2, 3:-1] - exact[2:-2, 1:-3]) / (2.0 * dx)
        + icoef * exact[2:-2, 2:-2]
    )
    return force


def check_sor_biharmonic2d_backend(pyfishpack: Any) -> dict[str, Any]:
    ny, nx = 30, 32
    dy, dx = 1.1, 0.9
    acoef, bcoef, ccoef = 0.08, 0.0, 0.08
    dcoef, ecoef, fcoef = -0.65, 0.0, -0.65
    gcoef, hcoef, icoef = 0.0, -0.03, 0.0
    y = np.arange(ny, dtype=np.float64)[:, None]
    x = np.arange(nx, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * y / (ny - 1.0)) * np.sin(math.pi * x / (nx - 1.0))
    force = _discrete_biharmonic2d_fixed(
        exact, dy, dx, acoef, bcoef, ccoef, dcoef, ecoef, fcoef, gcoef, hcoef, icoef
    )
    init = np.zeros_like(exact)
    init[:2, :] = exact[:2, :]
    init[-2:, :] = exact[-2:, :]
    init[:, :2] = exact[:, :2]
    init[:, -2:] = exact[:, -2:]
    solution, relerr, overflow, loops = pyfishpack.fishpack.sor_biharmonic2d(
        init,
        np.full_like(exact, acoef),
        np.full_like(exact, bcoef),
        np.full_like(exact, ccoef),
        np.full_like(exact, dcoef),
        np.full_like(exact, ecoef),
        np.full_like(exact, fcoef),
        np.full_like(exact, gcoef),
        np.full_like(exact, hcoef),
        np.full_like(exact, icoef),
        force,
        dy,
        dx,
        "fixed",
        "fixed",
        1.55,
        -9.99e8,
        50000,
        1e-12,
    )
    if overflow:
        raise AssertionError("sor_biharmonic2d overflowed")
    max_error = float(np.max(np.abs(solution - exact)))
    checked = _assert_close("sor_biharmonic2d manufactured solution", max_error, 0.0, 1e-8)
    checked.update({"relerr": float(relerr), "loops": int(loops)})
    return checked


def check_fftpack_transforms(pyfishpack: Any) -> dict[str, Any]:
    real = np.sin(np.linspace(0.0, 2.0 * math.pi, 24, endpoint=False)) + 0.25 * np.cos(
        np.linspace(0.0, 4.0 * math.pi, 24, endpoint=False)
    )
    coeff = pyfishpack.spectral_transform(real, kind="rfft", direction="forward")
    restored = pyfishpack.spectral_transform(coeff, kind="rfft", direction="inverse", normalize=True)
    real_error = float(np.max(np.abs(restored - real)))

    complex_input = np.exp(1j * np.linspace(0.0, 2.0 * math.pi, 18, endpoint=False))
    ccoeff = pyfishpack.spectral_transform(complex_input, kind="cfft", direction="forward")
    crestored = pyfishpack.spectral_transform(ccoeff, kind="cfft", direction="inverse", normalize=True)
    complex_error = float(np.max(np.abs(crestored - complex_input)))

    return {
        "rfft_roundtrip": _assert_close("FFTPACK rfft roundtrip", real_error, 0.0, 1e-12),
        "cfft_roundtrip": _assert_close("FFTPACK cfft roundtrip", complex_error, 0.0, 1e-12),
    }


def check_genbun_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 20, 40
    idimf = m + 2
    zero, one, two = 0.0, 1.0, 2.0
    dx = one / m
    dy = two * math.pi / n
    x = np.arange(m + 1, dtype=np.float64) * dx
    y = -math.pi + np.arange(n, dtype=np.float64) * dy

    a = np.zeros(m, dtype=np.float64)
    b = np.zeros(m, dtype=np.float64)
    c = np.zeros(m, dtype=np.float64)
    rhs = np.zeros((idimf, n), dtype=np.float64, order="F")

    s = (dy / dx) ** 2
    for i in range(1, m - 1):
        t = one + x[i]
        t2 = t * t
        a[i] = (t2 + t * dx) * s
        b[i] = -two * t2 * s
        c[i] = (t2 - t * dx) * s
    a[0] = zero
    b[0] = -two * s
    c[0] = -b[0]
    b[m - 1] = -two * s * (one + x[m - 1]) ** 2
    a[m - 1] = (-b[m - 1] / two) + (one + x[m - 1]) * dx * s
    c[m - 1] = zero

    dy2 = dy * dy
    for j in range(n):
        rhs[1 : m - 1, j] = 3.0 * (one + x[1 : m - 1]) ** 4 * dy2 * np.sin(y[j])
    t = one + x[m - 1]
    t2 = t * t
    t4 = t2 * t2
    rhs[0, :] = (11.0 + 8.0 / dx) * dy2 * np.sin(y)
    rhs[m - 1, :] = (
        3.0 * t4 * dy2 - 16.0 * t2 * s + 16.0 * t * s * dx
    ) * np.sin(y)

    solution, ierror = fishpack.genbun(0, n, 1, m, a, b, c, rhs)
    if ierror != 0:
        raise AssertionError(f"genbun returned ierror={ierror}")
    exact = ((one + x[:m])[:, None] ** 4) * np.sin(y)[None, :]
    error = float(np.max(np.abs(exact - solution[:m, :n])))
    return _assert_close("genbun UCAR example error", error, 0.964062912725572e-2, 5e-12)


def check_poistg_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 40, 20
    idimf = m + 2
    dx = math.pi / m
    dy = 1.0 / n
    x = -0.5 * math.pi + (np.arange(m, dtype=np.float64) + 0.5) * dx
    y = (np.arange(n, dtype=np.float64) + 0.5) * dy

    a = np.zeros(m, dtype=np.float64)
    b = np.zeros(m, dtype=np.float64)
    c = np.zeros(m, dtype=np.float64)
    rhs = np.zeros((idimf, n), dtype=np.float64, order="F")

    s = (dy / dx) ** 2
    half_dx = dx / 2.0
    a[0] = 0.0
    b[0] = -s * math.cos(-0.5 * math.pi + dx) / math.cos(x[0])
    c[0] = -b[0]
    for i in range(1, m):
        a[i] = s * math.cos(x[i] - half_dx) / math.cos(x[i])
        c[i] = s * math.cos(x[i] + half_dx) / math.cos(x[i])
        b[i] = -(a[i] + c[i])
    a[m - 1] = -b[m - 1]
    c[m - 1] = 0.0

    for j in range(n):
        rhs[:m, j] = 2.0 * (dy**2) * (y[j] ** 2) * (6.0 - y[j] ** 2) * np.sin(x)
    rhs[:m, n - 1] -= 4.0 * dy * np.sin(x)

    solution, ierror = fishpack.poistg(2, n, 1, m, a, b, c, rhs)
    if ierror != 0:
        raise AssertionError(f"poistg returned ierror={ierror}")
    exact = (y[None, :] ** 4) * np.sin(x)[:, None]
    error = float(np.max(np.abs(exact - solution[:m, :n])))
    return _assert_close("poistg UCAR example error", error, 0.564170618941665e-3, 5e-12)


def check_pois3d_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    l, m, n = 30, 30, 10
    ldimf, mdimf = l + 2, m + 3
    dx = 2.0 * math.pi / l
    dy = 2.0 * math.pi / m
    dz = 1.0 / n
    dz2 = 1.0 / dz**2
    c1 = 1.0 / dx**2
    c2 = 1.0 / dy**2
    x = -math.pi + np.arange(l, dtype=np.float64) * dx
    y = -math.pi + np.arange(m, dtype=np.float64) * dy
    z = np.zeros(n, dtype=np.float64)

    a = np.zeros(n, dtype=np.float64)
    b = np.zeros(n, dtype=np.float64)
    c = np.zeros(n, dtype=np.float64)
    f = np.zeros((ldimf, mdimf, n), dtype=np.float64, order="F")

    a[0] = 0.0
    b[0] = -2.0 * dz2
    c[0] = -b[0]
    for k in range(1, n):
        z[k] = k * dz
        t = 1.0 + z[k]
        t2 = t * t
        a[k] = t2 * dz2 + t / dz
        b[k] = -2.0 * t2 * dz2
        c[k] = t2 * dz2 - t / dz

    sinx = np.sin(x)
    siny = np.sin(y)
    for k in range(1, n):
        f[:l, :m, k] = 2.0 * sinx[:, None] * siny[None, :] * (1.0 + z[k]) ** 4
    f[:l, :m, 0] = (10.0 + 8.0 / dz) * sinx[:, None] * siny[None, :]
    f[:l, :m, n - 1] -= c[n - 1] * 16.0 * sinx[:, None] * siny[None, :]
    c[n - 1] = 0.0

    solution, ierror = fishpack.pois3d(0, l, c1, 0, m, c2, 1, n, a, b, c, f)
    if ierror != 0:
        raise AssertionError(f"pois3d returned ierror={ierror}")
    exact = (
        sinx[:, None, None]
        * siny[None, :, None]
        * (1.0 + z[None, None, :]) ** 4
    )
    error = float(np.max(np.abs(exact - solution[:l, :m, :n])))
    return _assert_close("pois3d UCAR example error", error, 0.293277049861086e-1, 5e-12)


def check_hwscrt_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 40, 80
    idimf = m + 5
    a_domain, b_domain = 0.0, 2.0
    c_domain, d_domain = -1.0, 3.0
    mbdcnd, nbdcnd = 2, 0
    elmbda = -4.0
    x = 2.0 * np.arange(m + 1, dtype=np.float64) / m
    y = -1.0 + 4.0 * np.arange(n + 1, dtype=np.float64) / n
    bda = np.zeros(1, dtype=np.float64)
    bdb = 4.0 * np.cos((y + 1.0) * (math.pi / 2.0))
    bdc = np.zeros(1, dtype=np.float64)
    bdd = np.zeros(1, dtype=np.float64)
    rhs = np.zeros((idimf, n + 1), dtype=np.float64, order="F")

    pi2 = math.pi**2
    for j in range(n + 1):
        rhs[1 : m + 1, j] = (
            2.0 - (4.0 + pi2 / 4.0) * (x[1 : m + 1] ** 2)
        ) * math.cos((y[j] + 1.0) * (math.pi / 2.0))

    solution, pertrb, ierror = fishpack.hwscrt(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hwscrt returned ierror={ierror}")
    if abs(pertrb) > 5e-14:
        raise AssertionError(f"hwscrt returned unexpected pertrb={pertrb}")
    exact = (x[:, None] ** 2) * np.cos((y[None, :] + 1.0) * (math.pi / 2.0))
    error = float(np.max(np.abs(exact - solution[: m + 1, : n + 1])))
    return _assert_close("hwscrt UCAR example error", error, 0.536508246868017e-3, 5e-12)


def check_hstcrt_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 48, 53
    idimf = m + 2
    a_domain, b_domain = 1.0, 3.0
    c_domain, d_domain = -1.0, 1.0
    mbdcnd, nbdcnd = 2, 0
    elmbda = -2.0
    dx = (b_domain - a_domain) / m
    dy = (d_domain - c_domain) / n
    x = a_domain + (np.arange(m, dtype=np.float64) + 0.5) * dx
    y = c_domain + (np.arange(n, dtype=np.float64) + 0.5) * dy

    bda = np.zeros(n, dtype=np.float64)
    bdb = -math.pi * np.cos(math.pi * y)
    bdc = np.zeros(1, dtype=np.float64)
    bdd = np.zeros(1, dtype=np.float64)
    rhs = np.zeros((idimf, n), dtype=np.float64, order="F")
    rhs[:m, :] = -2.0 * (math.pi**2 + 1.0) * np.sin(math.pi * x)[:, None] * np.cos(
        math.pi * y
    )[None, :]

    solution, pertrb, ierror = fishpack.hstcrt(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hstcrt returned ierror={ierror}")
    if abs(pertrb) > 5e-14:
        raise AssertionError(f"hstcrt returned unexpected pertrb={pertrb}")
    exact = np.sin(math.pi * x)[:, None] * np.cos(math.pi * y)[None, :]
    error = float(np.max(np.abs(exact - solution[:m, :n])))
    return _assert_close("hstcrt UCAR example error", error, 0.126000790151959e-2, 5e-12)


def check_hw3crt_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    l, m, n = 10, 40, 15
    lp1, mp1, np1 = l + 1, m + 1, n + 1
    xs, xf = 0.0, 1.0
    ys, yf = 0.0, 2.0 * math.pi
    zs, zf = 0.0, 0.5 * math.pi
    lbdcnd, mbdcnd, nbdcnd = 1, 0, 2
    elmbda = -3.0
    dx = (xf - xs) / l
    dy = (yf - ys) / m
    dz = (zf - zs) / n
    del dx, dy, dz

    x = xs + np.arange(lp1, dtype=np.float64) * ((xf - xs) / l)
    y = ys + np.arange(mp1, dtype=np.float64) * ((yf - ys) / m)
    z = zs + np.arange(np1, dtype=np.float64) * ((zf - zs) / n)

    dummy = np.zeros((1, 1), dtype=np.float64, order="F")
    bdzf = np.asfortranarray(-(x[:, None] ** 4) * np.sin(y)[None, :])
    rhs = np.zeros((lp1, mp1, np1), dtype=np.float64, order="F")

    for k in range(np1):
        rhs[0, :, k] = 0.0
        rhs[l, :, k] = np.sin(y) * np.cos(z[k])
    rhs[:, :, 0] = (x[:, None] ** 4) * np.sin(y)[None, :]
    rhs[1:l, :, 1:np1] = (
        4.0
        * (x[1:l, None, None] ** 2)
        * (3.0 - x[1:l, None, None] ** 2)
        * np.sin(y)[None, :, None]
        * np.cos(z)[None, None, 1:np1]
    )

    solution, pertrb, ierror = fishpack.hw3crt(
        xs,
        xf,
        l,
        lbdcnd,
        dummy,
        dummy,
        ys,
        yf,
        m,
        mbdcnd,
        dummy,
        dummy,
        zs,
        zf,
        n,
        nbdcnd,
        dummy,
        bdzf,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hw3crt returned ierror={ierror}")
    if abs(pertrb) > 5e-14:
        raise AssertionError(f"hw3crt returned unexpected pertrb={pertrb}")
    exact = (x[:, None, None] ** 4) * np.sin(y)[None, :, None] * np.cos(
        z
    )[None, None, :]
    error = float(np.max(np.abs(exact - solution[:lp1, :mp1, :np1])))
    return _assert_close("hw3crt UCAR example error", error, 0.964801952068239e-2, 5e-12)


def check_hwsplr_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 50, 48
    idimf = 100
    mp1, np1 = m + 1, n + 1
    a_domain, b_domain = 0.0, 1.0
    c_domain, d_domain = 0.0, 0.5 * math.pi
    mbdcnd, nbdcnd = 5, 3
    elmbda = 0.0
    r = np.arange(mp1, dtype=np.float64) / m
    theta = np.arange(np1, dtype=np.float64) * ((0.5 * math.pi) / n)
    bda = np.zeros(1, dtype=np.float64)
    bdb = np.zeros(1, dtype=np.float64)
    bdc = np.zeros(mp1, dtype=np.float64)
    bdd = np.zeros(mp1, dtype=np.float64)
    rhs = np.zeros((idimf, m), dtype=np.float64, order="F")
    rhs[mp1 - 1, :np1] = 1.0 - np.cos(4.0 * theta)
    rhs[:m, :np1] = 16.0 * r[:m, None] ** 2

    solution, pertrb, ierror = fishpack.hwsplr(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hwsplr returned ierror={ierror}")
    if abs(pertrb) > 5e-14:
        raise AssertionError(f"hwsplr returned unexpected pertrb={pertrb}")
    exact = (r[:, None] ** 4) * (1.0 - np.cos(4.0 * theta)[None, :])
    error = float(np.max(np.abs(exact - solution[:mp1, :np1])))
    return _assert_close("hwsplr UCAR example error", error, 0.619134227874629e-3, 5e-12)


def check_hwscyl_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 50, 100
    idimf = m + 25
    mp1, np1 = m + 1, n + 1
    a_domain, b_domain = 0.0, 1.0
    c_domain, d_domain = 0.0, 1.0
    mbdcnd, nbdcnd = 6, 3
    elmbda = 0.0
    r = np.arange(mp1, dtype=np.float64) / m
    z = np.arange(np1, dtype=np.float64) / n
    bda = np.zeros(np1, dtype=np.float64)
    bdb = 4.0 * z**4
    bdc = np.zeros(mp1, dtype=np.float64)
    bdd = 4.0 * r**4
    rhs = np.zeros((idimf, np1), dtype=np.float64, order="F")
    rhs[:mp1, :np1] = (
        4.0
        * (r[:, None] ** 2)
        * (z[None, :] ** 2)
        * (4.0 * z[None, :] ** 2 + 3.0 * r[:, None] ** 2)
    )

    solution, pertrb, ierror = fishpack.hwscyl(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hwscyl returned ierror={ierror}")
    exact = (r[:, None] * z[None, :]) ** 4
    adjusted = np.array(solution[:mp1, :np1], copy=True)
    adjusted -= float(np.mean(adjusted - exact))
    error = float(np.max(np.abs(exact - adjusted)))
    pertrb_check = _assert_close("hwscyl UCAR perturbation", pertrb, 0.226742668667313e-3, 5e-12)
    error_check = _assert_close("hwscyl UCAR example error", error, 0.373672238079603e-3, 5e-12)
    return {"error": error_check, "pertrb": pertrb_check}


def check_hstplr_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 50, 48
    idimf = m + 1
    np1 = n + 1
    a_domain, b_domain = 0.0, 1.0
    c_domain, d_domain = 0.0, 0.5 * math.pi
    mbdcnd, nbdcnd = 5, 3
    elmbda = 0.0
    dr = (b_domain - a_domain) / m
    dtheta = (d_domain - c_domain) / n
    r = (np.arange(m, dtype=np.float64) + 0.5) * dr
    theta = (np.arange(n, dtype=np.float64) + 0.5) * dtheta
    bda = np.zeros(1, dtype=np.float64)
    bdb = 1.0 - np.cos(4.0 * theta)
    bdc = np.zeros(m, dtype=np.float64)
    bdd = np.zeros(m, dtype=np.float64)
    rhs = np.zeros((idimf, np1), dtype=np.float64, order="F")
    rhs[:m, :n] = 16.0 * r[:, None] ** 2

    solution, pertrb, ierror = fishpack.hstplr(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hstplr returned ierror={ierror}")
    if abs(pertrb) > 5e-14:
        raise AssertionError(f"hstplr returned unexpected pertrb={pertrb}")
    exact = (r[:, None] ** 4) * (1.0 - np.cos(4.0 * theta)[None, :])
    error = float(np.max(np.abs(exact - solution[:m, :n])))
    return _assert_close("hstplr UCAR example error", error, 0.113037945648764e-2, 5e-12)


def check_hstcyl_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 50, 52
    idimf = m + 1
    a_domain, b_domain = 0.0, 1.0
    c_domain, d_domain = 0.0, 1.0
    mbdcnd, nbdcnd = 6, 3
    elmbda = 0.0
    r = (np.arange(m, dtype=np.float64) + 0.5) / m
    z = (np.arange(n, dtype=np.float64) + 0.5) / n
    bda = np.zeros(1, dtype=np.float64)
    bdb = 4.0 * z**4
    bdc = np.zeros(m, dtype=np.float64)
    bdd = 4.0 * r**4
    rhs = np.zeros((idimf, n), dtype=np.float64, order="F")
    rhs[:m, :n] = (
        4.0
        * (r[:, None] ** 2)
        * (z[None, :] ** 2)
        * (4.0 * z[None, :] ** 2 + 3.0 * r[:, None] ** 2)
    )

    solution, pertrb, ierror = fishpack.hstcyl(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hstcyl returned ierror={ierror}")
    exact = (r[:, None] * z[None, :]) ** 4
    adjusted = np.array(solution[:m, :n], copy=True)
    adjusted -= float(np.mean(adjusted - exact))
    error = float(np.max(np.abs(exact - adjusted)))
    pertrb_check = _assert_close("hstcyl UCAR perturbation", pertrb, -0.443113920336705e-3, 5e-12)
    error_check = _assert_close("hstcyl UCAR example error", error, 0.752796331450795e-4, 5e-12)
    return {"error": error_check, "pertrb": pertrb_check}


def check_hwsssp_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 18, 72
    mp1, np1 = m + 1, n + 1
    idimf = mp1
    ts, tf = 0.0, 0.5 * math.pi
    ps, pf = 0.0, 2.0 * math.pi
    mbdcnd, nbdcnd = 6, 0
    elmbda = 0.0
    dtheta = tf / m
    dphi = (2.0 * math.pi) / n
    sint = np.sin(np.arange(mp1, dtype=np.float64) * dtheta)
    sinp = np.sin(np.arange(np1, dtype=np.float64) * dphi)
    bdts = np.zeros(np1, dtype=np.float64)
    bdtf = np.zeros(np1, dtype=np.float64)
    bdps = np.zeros(1, dtype=np.float64)
    bdpf = np.zeros(1, dtype=np.float64)
    rhs = np.zeros((idimf, np1), dtype=np.float64, order="F")
    rhs[:mp1, :np1] = 2.0 - 6.0 * (sint[:, None] * sinp[None, :]) ** 2

    solution, pertrb, ierror = fishpack.hwsssp(
        ts,
        tf,
        m,
        mbdcnd,
        bdts,
        bdtf,
        ps,
        pf,
        n,
        nbdcnd,
        bdps,
        bdpf,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hwsssp returned ierror={ierror}")
    exact = (sint[:, None] * sinp[None, :]) ** 2 - solution[0, 0]
    error = float(np.max(np.abs(exact - solution[:mp1, :np1])))
    error_check = _assert_close("hwsssp UCAR example error", error, 0.338107338525839e-2, 5e-12)
    return {"error": error_check, "pertrb": pertrb}


def check_hstssp_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 18, 72
    idimf = m
    a_domain, b_domain = 0.0, 0.5 * math.pi
    c_domain, d_domain = 0.0, 2.0 * math.pi
    mbdcnd, nbdcnd = 6, 0
    elmbda = 0.0
    dtheta = (b_domain - a_domain) / m
    dphi = (d_domain - c_domain) / n
    sint = np.sin((np.arange(m, dtype=np.float64) + 0.5) * dtheta)
    sinp = np.sin((np.arange(n, dtype=np.float64) + 0.5) * dphi)
    bda = np.zeros(1, dtype=np.float64)
    bdb = np.zeros(n, dtype=np.float64)
    bdc = np.zeros(1, dtype=np.float64)
    bdd = np.zeros(1, dtype=np.float64)
    rhs = np.zeros((idimf, n), dtype=np.float64, order="F")
    rhs[:m, :n] = 2.0 - 6.0 * (sint[:, None] * sinp[None, :]) ** 2

    solution, pertrb, ierror = fishpack.hstssp(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        elmbda,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hstssp returned ierror={ierror}")
    exact = (sint[:, None] * sinp[None, :]) ** 2 + solution[0, 0]
    error = float(np.max(np.abs(exact - solution[:m, :n])))
    pertrb_check = _assert_close("hstssp UCAR perturbation", pertrb, 0.635830001454109e-3, 5e-12)
    error_check = _assert_close("hstssp UCAR example error", error, 0.337523232257420e-2, 5e-12)
    return {"error": error_check, "pertrb": pertrb_check}


def check_hwscsp_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 36, 32
    mp1, np1 = m + 1, n + 1
    idimf = m + 12
    ts, tf = 0.0, 0.5 * math.pi
    rs, rf = 0.0, 1.0
    dtheta = tf / m
    dr = 1.0 / n
    theta = np.arange(mp1, dtype=np.float64) * dtheta
    r = np.arange(np1, dtype=np.float64) * dr
    bdts = np.zeros(1, dtype=np.float64)
    bdtf = np.zeros(np1, dtype=np.float64)
    bdrs = np.zeros(1, dtype=np.float64)
    bdrf = np.zeros(1, dtype=np.float64)

    rhs1 = np.zeros((idimf, np1), dtype=np.float64, order="F")
    rhs1[:mp1, n] = np.cos(theta) ** 4
    rhs1[:mp1, :n] = 12.0 * (np.cos(theta)[:, None] ** 2) * (r[None, :n] ** 2)
    solution1, pertrb1, ierror1 = fishpack.hwscsp(
        ts,
        tf,
        m,
        6,
        bdts,
        bdtf,
        rs,
        rf,
        n,
        5,
        bdrs,
        bdrf,
        0.0,
        rhs1,
    )
    if ierror1 != 0:
        raise AssertionError(f"hwscsp example 1 returned ierror={ierror1}")
    if abs(pertrb1) > 5e-14:
        raise AssertionError(f"hwscsp example 1 returned unexpected pertrb={pertrb1}")
    exact1 = (np.cos(theta)[:, None] ** 4) * (r[None, :n] ** 4)
    error1 = float(np.max(np.abs(exact1 - solution1[:mp1, :n])))

    dphi = math.pi / 72.0
    elmbda = -2.0 * (1.0 - math.cos(dphi)) / dphi**2
    rhs2 = np.zeros((idimf, np1), dtype=np.float64, order="F")
    rhs2[:mp1, n] = np.sin(theta)
    solution2, pertrb2, ierror2 = fishpack.hwscsp(
        ts,
        tf,
        m,
        2,
        bdts,
        bdtf,
        rs,
        rf,
        n,
        1,
        bdrs,
        bdrf,
        elmbda,
        rhs2,
    )
    if ierror2 != 0:
        raise AssertionError(f"hwscsp example 2 returned ierror={ierror2}")
    if abs(pertrb2) > 5e-14:
        raise AssertionError(f"hwscsp example 2 returned unexpected pertrb={pertrb2}")
    exact2 = np.sin(theta)[:, None] * r[None, :]
    error2 = float(np.max(np.abs(exact2 - solution2[:mp1, :np1])))

    return {
        "example_1": _assert_close(
            "hwscsp UCAR example 1 error", error1, 0.799841637481730e-3, 5e-12
        ),
        "example_2": _assert_close(
            "hwscsp UCAR example 2 error", error2, 0.586824289033339e-4, 5e-12
        ),
    }


def check_hstcsp_example(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    m, n = 45, 15
    idimf = m + 2
    np1 = n + 1
    a_domain, b_domain = 0.0, math.pi
    c_domain, d_domain = 0.0, 1.0
    mbdcnd, nbdcnd = 9, 5
    dt = (b_domain - a_domain) / m
    dr = (d_domain - c_domain) / n
    theta = a_domain + (np.arange(m, dtype=np.float64) + 0.5) * dt
    cost = np.cos(theta)
    r = c_domain + (np.arange(n, dtype=np.float64) + 0.5) * dr
    bda = np.zeros(1, dtype=np.float64)
    bdb = np.zeros(1, dtype=np.float64)
    bdc = np.zeros(1, dtype=np.float64)
    bdd = cost**4
    rhs = np.zeros((idimf, np1), dtype=np.float64, order="F")
    rhs[:m, :n] = 12.0 * (r[None, :] * cost[:, None]) ** 2

    solution, pertrb, ierror = fishpack.hstcsp(
        a_domain,
        b_domain,
        m,
        mbdcnd,
        bda,
        bdb,
        c_domain,
        d_domain,
        n,
        nbdcnd,
        bdc,
        bdd,
        0.0,
        rhs,
    )
    if ierror != 0:
        raise AssertionError(f"hstcsp returned ierror={ierror}")
    if abs(pertrb) > 5e-14:
        raise AssertionError(f"hstcsp returned unexpected pertrb={pertrb}")
    exact = (r[None, :] * cost[:, None]) ** 4
    error = float(np.max(np.abs(exact - solution[:m, :n])))
    return _assert_close("hstcsp UCAR example error", error, 0.558432375660800e-2, 5e-12)


def check_high_level_poisson(pyfishpack: Any) -> dict[str, Any]:
    m, n = 96, 96
    dy, dx = 0.7, 1.3
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))

    lap = np.zeros_like(exact)
    lap += -2.0 * exact / dy**2
    lap += -2.0 * exact / dx**2
    lap[1:, :] += exact[:-1, :] / dy**2
    lap[:-1, :] += exact[1:, :] / dy**2
    lap[:, 1:] += exact[:, :-1] / dx**2
    lap[:, :-1] += exact[:, 1:] / dx**2

    solution = pyfishpack.invert_Poisson(lap, BCs=("fixed", "fixed"), spacing=(dy, dx))
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_Poisson manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_multigrid(pyfishpack: Any) -> dict[str, Any]:
    m, n = 54, 52
    dy, dx = 0.8, 1.2
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_laplacian_fixed(exact, dy, dx)

    solution, grids, history = pyfishpack.invert_MultiGrid(
        pyfishpack.invert_Poisson,
        forcing,
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    if len(grids) != 1 or len(history) != 1:
        raise AssertionError("invert_MultiGrid compatibility wrapper returned unexpected history")
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_MultiGrid direct solution", max_error, 0.0, 5e-12)


def _discrete_laplacian_fixed(exact: np.ndarray, dy: float, dx: float) -> np.ndarray:
    lap = np.zeros_like(exact)
    lap += -2.0 * exact / dy**2
    lap += -2.0 * exact / dx**2
    lap[1:, :] += exact[:-1, :] / dy**2
    lap[:-1, :] += exact[1:, :] / dy**2
    lap[:, 1:] += exact[:, :-1] / dx**2
    lap[:, :-1] += exact[:, 1:] / dx**2
    return lap


def _discrete_constant_2d_fixed(
    exact: np.ndarray,
    dy: float,
    dx: float,
    alpha_y: float,
    alpha_x: float,
    helmholtz: float = 0.0,
) -> np.ndarray:
    op = np.zeros_like(exact)
    op += -2.0 * alpha_y * exact / dy**2
    op += -2.0 * alpha_x * exact / dx**2
    op += helmholtz * exact
    op[1:, :] += alpha_y * exact[:-1, :] / dy**2
    op[:-1, :] += alpha_y * exact[1:, :] / dy**2
    op[:, 1:] += alpha_x * exact[:, :-1] / dx**2
    op[:, :-1] += alpha_x * exact[:, 1:] / dx**2
    return op


def _discrete_standard_2d_fixed(
    exact: np.ndarray,
    dy: float,
    dx: float,
    a: np.ndarray,
    c: np.ndarray,
) -> np.ndarray:
    op = np.zeros_like(exact)
    op[1:-1, 1:-1] = (
        (
            a[2:, 1:-1] * (exact[2:, 1:-1] - exact[1:-1, 1:-1])
            - a[1:-1, 1:-1] * (exact[1:-1, 1:-1] - exact[:-2, 1:-1])
        )
        / (dy * dy)
        + (
            c[1:-1, 2:] * (exact[1:-1, 2:] - exact[1:-1, 1:-1])
            - c[1:-1, 1:-1] * (exact[1:-1, 1:-1] - exact[1:-1, :-2])
        )
        / (dx * dx)
    )
    return op


def _discrete_general_2d_fixed(
    exact: np.ndarray,
    dy: float,
    dx: float,
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
    e: np.ndarray,
    fcoef: np.ndarray,
) -> np.ndarray:
    op = np.zeros_like(exact)
    op[1:-1, 1:-1] = (
        a[1:-1, 1:-1]
        * (exact[2:, 1:-1] - 2.0 * exact[1:-1, 1:-1] + exact[:-2, 1:-1])
        / (dy * dy)
        + b[1:-1, 1:-1]
        * (
            exact[2:, 2:]
            - exact[:-2, 2:]
            - exact[2:, :-2]
            + exact[:-2, :-2]
        )
        / (4.0 * dy * dx)
        + c[1:-1, 1:-1]
        * (exact[1:-1, 2:] - 2.0 * exact[1:-1, 1:-1] + exact[1:-1, :-2])
        / (dx * dx)
        + d[1:-1, 1:-1]
        * (exact[2:, 1:-1] - exact[:-2, 1:-1])
        / (2.0 * dy)
        + e[1:-1, 1:-1]
        * (exact[1:-1, 2:] - exact[1:-1, :-2])
        / (2.0 * dx)
        + fcoef[1:-1, 1:-1] * exact[1:-1, 1:-1]
    )
    return op


def _discrete_omega_operator_fixed(
    exact: np.ndarray, dz: float, dy: float, dx: float, alpha_z: float
) -> np.ndarray:
    return _discrete_constant_3d_fixed(exact, dz, dy, dx, alpha_z, 1.0, 1.0)


def _discrete_constant_3d_fixed(
    exact: np.ndarray,
    dz: float,
    dy: float,
    dx: float,
    alpha_z: float,
    alpha_y: float,
    alpha_x: float,
) -> np.ndarray:
    op = np.zeros_like(exact)
    op += -2.0 * alpha_z * exact / dz**2
    op += -2.0 * alpha_y * exact / dy**2
    op += -2.0 * alpha_x * exact / dx**2
    op[1:, :, :] += alpha_z * exact[:-1, :, :] / dz**2
    op[:-1, :, :] += alpha_z * exact[1:, :, :] / dz**2
    op[:, 1:, :] += alpha_y * exact[:, :-1, :] / dy**2
    op[:, :-1, :] += alpha_y * exact[:, 1:, :] / dy**2
    op[:, :, 1:] += alpha_x * exact[:, :, :-1] / dx**2
    op[:, :, :-1] += alpha_x * exact[:, :, 1:] / dx**2
    return op


def check_high_level_geostrophic(pyfishpack: Any) -> dict[str, Any]:
    m, n = 80, 76
    dy, dx = 0.9, 1.1
    f0 = 1.3e-4
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    lap_phi = f0 * _discrete_laplacian_fixed(exact, dy, dx)

    solution = pyfishpack.invert_geostrophic(
        lap_phi,
        coords="cartesian",
        mParams={"f0": f0, "beta": 0.0},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_geostrophic manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_geostrophic_beta(pyfishpack: Any) -> dict[str, Any]:
    m, n = 44, 42
    dy, dx = 0.8, 1.05
    f0 = 0.7
    beta = 0.02
    y = np.arange(m, dtype=np.float64)[:, None] * dy
    y_half = np.arange(m, dtype=np.float64)[:, None] * dy
    y_half[1:, :] = 0.5 * (y[1:, :] + y[:-1, :])
    a = np.broadcast_to(f0 + beta * y_half, (m, n))
    c = np.broadcast_to(f0 + beta * y, (m, n))
    jj = np.arange(m, dtype=np.float64)[:, None]
    ii = np.arange(n, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * jj / (m - 1.0)) * np.sin(math.pi * ii / (n - 1.0))
    forcing = _discrete_standard_2d_fixed(exact, dy, dx, a, c)

    solution = pyfishpack.invert_geostrophic(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": beta},
        iParams={"mxLoop": 50000, "tolerance": 1e-12},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_geostrophic beta-plane manufactured solution",
        max_error,
        0.0,
        1e-8,
    )


def check_high_level_pv2d(pyfishpack: Any) -> dict[str, Any]:
    m, n = 70, 66
    dy, dx = 0.75, 1.25
    f0 = 1.15e-4
    n2 = 2.3e-4
    alpha_y = f0 * f0 / n2
    alpha_x = 1.0
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_constant_2d_fixed(exact, dy, dx, alpha_y, alpha_x)

    solution = pyfishpack.invert_PV2D(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "N2": n2},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_PV2D manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_eliassen(pyfishpack: Any) -> dict[str, Any]:
    m, n = 68, 64
    dy, dx = 1.05, 0.85
    alpha_y = 2.4
    alpha_x = 3.1
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_constant_2d_fixed(exact, dy, dx, alpha_y, alpha_x)

    solution = pyfishpack.invert_Eliassen(
        forcing,
        coords="cartesian",
        mParams={"A": alpha_y, "B": 0.0, "C": alpha_x},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_Eliassen manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_eliassen_cross(pyfishpack: Any) -> dict[str, Any]:
    m, n = 40, 42
    dy, dx = 0.9, 1.2
    alpha_y = 2.2
    cross = 0.35
    alpha_x = 2.7
    y = np.arange(m, dtype=np.float64)[:, None]
    x = np.arange(n, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * y / (m - 1.0)) * np.sin(math.pi * x / (n - 1.0))
    a = np.full_like(exact, alpha_y)
    b = np.full_like(exact, 2.0 * cross)
    c = np.full_like(exact, alpha_x)
    zeros = np.zeros_like(exact)
    forcing = _discrete_general_2d_fixed(exact, dy, dx, a, b, c, zeros, zeros, zeros)

    solution = pyfishpack.invert_Eliassen(
        forcing,
        coords="cartesian",
        mParams={"A": alpha_y, "B": cross, "C": alpha_x},
        iParams={"mxLoop": 30000, "tolerance": 1e-11},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_Eliassen cross-derivative manufactured solution", max_error, 0.0, 1e-8
    )


def check_high_level_bretherton(pyfishpack: Any) -> dict[str, Any]:
    m, n = 72, 68
    dy, dx = 1.2, 0.8
    f0 = 1.0e-4
    depth = 250.0
    lamb = 3.0e-6
    helmholtz = -lamb * depth
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_laplacian_fixed(exact, dy, dx) + helmholtz * exact
    topography = -forcing * depth / f0

    solution = pyfishpack.invert_BrethertonHaidvogel(
        topography,
        coords="cartesian",
        mParams={"f0": f0, "beta": 0.0, "D": depth, "lambda": lamb},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_BrethertonHaidvogel manufactured solution", max_error, 0.0, 5e-12
    )


def check_high_level_bretherton_beta(pyfishpack: Any) -> dict[str, Any]:
    m, n = 42, 38
    dy, dx = 0.9, 1.1
    f0 = 0.8
    beta = 0.03
    depth = 80.0
    lamb = 2.0e-3
    helmholtz = -lamb * depth
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    y = np.arange(m, dtype=np.float64)[:, None] * dy
    coriolis = f0 + beta * y
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_laplacian_fixed(exact, dy, dx) + helmholtz * exact
    topography = -forcing * depth / coriolis

    solution = pyfishpack.invert_BrethertonHaidvogel(
        topography,
        coords="cartesian",
        mParams={"f0": f0, "beta": beta, "D": depth, "lambda": lamb},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_BrethertonHaidvogel beta-plane manufactured solution",
        max_error,
        0.0,
        5e-11,
    )


def check_high_level_fofonoff(pyfishpack: Any) -> dict[str, Any]:
    m, n = 64, 60
    field = np.zeros((m, n), dtype=np.float64)
    solution = pyfishpack.invert_Fofonoff(
        field,
        coords="cartesian",
        mParams={"f0": 8.0e-5, "beta": 0.0, "c0": 2.0e-5, "c1": 8.0e-5},
        BCs=("fixed", "fixed"),
    )
    max_abs = float(np.max(np.abs(solution)))
    return _assert_close("invert_Fofonoff zero-forcing solution", max_abs, 0.0, 5e-12)


def check_high_level_fofonoff_beta(pyfishpack: Any) -> dict[str, Any]:
    m, n = 44, 40
    dy, dx = 0.7, 1.1
    f0 = 0.5
    beta = 0.02
    c0 = 0.03
    c1 = 1.1
    field = np.zeros((m, n), dtype=np.float64)
    y = np.arange(m, dtype=np.float64)[:, None] * dy
    forcing = np.broadcast_to(c1 - (f0 + beta * y), field.shape)
    solution = pyfishpack.invert_Fofonoff(
        field,
        coords="cartesian",
        mParams={"f0": f0, "beta": beta, "c0": c0, "c1": c1},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    residual = _discrete_constant_2d_fixed(solution, dy, dx, 1.0, 1.0, -c0) - forcing
    max_residual = float(np.max(np.abs(residual)))
    return _assert_close("invert_Fofonoff beta-plane residual", max_residual, 0.0, 5e-11)


def check_high_level_gill_matsuno(pyfishpack: Any) -> dict[str, Any]:
    m, n = 62, 58
    dy, dx = 0.95, 1.15
    epsilon = 1.1e-5
    phi = 9.8
    f0 = 1.0e-4
    alpha = epsilon * phi / (epsilon * epsilon + f0 * f0)
    helmholtz = -epsilon
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_constant_2d_fixed(exact, dy, dx, alpha, alpha, helmholtz)

    solution = pyfishpack.invert_GillMatsuno(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": 0.0, "epsilon": epsilon, "Phi": phi},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_GillMatsuno manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_gill_matsuno_beta(pyfishpack: Any) -> dict[str, Any]:
    m, n = 38, 40
    dy, dx = 0.8, 1.1
    epsilon = 1.4
    phi = 2.3
    f0 = 0.7
    beta = 0.03
    y = np.arange(m, dtype=np.float64)[:, None] * dy
    x = np.arange(n, dtype=np.float64)[None, :] * dx
    exact = np.sin(math.pi * y / ((m - 1.0) * dy)) * np.sin(
        math.pi * x / ((n - 1.0) * dx)
    )
    coriolis = f0 + beta * y
    denom = epsilon * epsilon + coriolis * coriolis
    c1 = epsilon / denom
    c1_dy = -2.0 * epsilon * coriolis * beta / (denom * denom)
    c2_dy = beta * (epsilon * epsilon - coriolis * coriolis) / (denom * denom)
    a = np.broadcast_to(phi * c1, exact.shape)
    b = np.zeros_like(exact)
    c = np.broadcast_to(phi * c1, exact.shape)
    d = np.broadcast_to(phi * c1_dy, exact.shape)
    e = np.broadcast_to(-phi * c2_dy, exact.shape)
    fcoef = np.full_like(exact, -epsilon)
    forcing = _discrete_general_2d_fixed(exact, dy, dx, a, b, c, d, e, fcoef)

    solution = pyfishpack.invert_GillMatsuno(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": beta, "epsilon": epsilon, "Phi": phi},
        iParams={"mxLoop": 30000, "tolerance": 1e-11},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_GillMatsuno beta-plane manufactured solution", max_error, 0.0, 1e-7
    )


def check_high_level_gill_matsuno_test(pyfishpack: Any) -> dict[str, Any]:
    m, n = 58, 54
    dy, dx = 1.05, 0.9
    epsilon = 1.3e-5
    phi = 8.7
    f0 = 1.15e-4
    alpha = epsilon * phi / (epsilon * epsilon + f0 * f0)
    helmholtz = -epsilon
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_constant_2d_fixed(exact, dy, dx, alpha, alpha, helmholtz)

    solution = pyfishpack.invert_GillMatsuno_test(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": 0.0, "epsilon": epsilon, "Phi": phi},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_GillMatsuno_test manufactured solution", max_error, 0.0, 5e-12
    )


def check_high_level_stommel(pyfishpack: Any) -> dict[str, Any]:
    m, n = 66, 60
    dy, dx = 1.3, 0.7
    resistance = 2.5
    depth = 12.0
    rho0 = 1020.0
    alpha = -resistance / depth
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    rhs = _discrete_constant_2d_fixed(exact, dy, dx, alpha, alpha)
    curl = -rhs * depth * rho0

    solution = pyfishpack.invert_Stommel(
        curl,
        coords="cartesian",
        mParams={"beta": 0.0, "R": resistance, "D": depth, "rho0": rho0},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_Stommel manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_stommel_beta(pyfishpack: Any) -> dict[str, Any]:
    m, n = 34, 36
    dy, dx = 1.1, 0.9
    resistance = 3.0
    depth = 12.0
    rho0 = 1020.0
    beta = 0.015
    alpha = -resistance / depth
    y = np.arange(m, dtype=np.float64)[:, None]
    x = np.arange(n, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * y / (m - 1.0)) * np.sin(math.pi * x / (n - 1.0))
    rhs = np.zeros_like(exact)
    rhs[1:-1, 1:-1] = (
        alpha
        * (
            (exact[2:, 1:-1] - 2.0 * exact[1:-1, 1:-1] + exact[:-2, 1:-1])
            / (dy * dy)
            + (exact[1:-1, 2:] - 2.0 * exact[1:-1, 1:-1] + exact[1:-1, :-2])
            / (dx * dx)
        )
        - beta * (exact[1:-1, 2:] - exact[1:-1, :-2]) / (2.0 * dx)
    )
    curl = -rhs * depth * rho0

    solution = pyfishpack.invert_Stommel(
        curl,
        coords="cartesian",
        mParams={"beta": beta, "R": resistance, "D": depth, "rho0": rho0},
        iParams={"mxLoop": 20000, "tolerance": 1e-11},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_Stommel beta-plane manufactured solution", max_error, 0.0, 1e-8
    )


def check_high_level_stommel_test(pyfishpack: Any) -> dict[str, Any]:
    m, n = 64, 62
    dy, dx = 0.8, 1.2
    resistance = 1.8
    depth = 11.0
    rho0 = 1023.0
    alpha = -resistance / depth
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    rhs = _discrete_constant_2d_fixed(exact, dy, dx, alpha, alpha)
    curl = -rhs * depth * rho0

    solution = pyfishpack.invert_Stommel_test(
        curl,
        coords="cartesian",
        mParams={"f0": 1.0e-4, "beta": 0.0, "R": resistance, "D": depth, "rho0": rho0},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_Stommel_test manufactured solution", max_error, 0.0, 5e-12
    )


def check_high_level_stommel_munk(pyfishpack: Any) -> dict[str, Any]:
    m, n = 60, 58
    dy, dx = 0.9, 1.1
    resistance = 2.1
    depth = 13.0
    rho0 = 1024.0
    alpha = -resistance / depth
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    rhs = _discrete_constant_2d_fixed(exact, dy, dx, alpha, alpha)
    curl = -rhs * depth * rho0

    solution = pyfishpack.invert_StommelMunk(
        curl,
        coords="cartesian",
        mParams={
            "A4": 0.0,
            "beta": 0.0,
            "R": resistance,
            "D": depth,
            "rho0": rho0,
        },
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    accepted = _assert_close(
        "invert_StommelMunk A4=0 manufactured solution", max_error, 0.0, 5e-12
    )

    a4 = 0.06
    beta = 0.018
    rhs_bih = _discrete_biharmonic2d_fixed(
        exact,
        dy,
        dx,
        a4,
        0.0,
        a4,
        -resistance / depth,
        0.0,
        -resistance / depth,
        0.0,
        -beta,
        0.0,
    )
    curl_bih = -rhs_bih * depth * rho0
    init = np.zeros_like(exact)
    init[:2, :] = exact[:2, :]
    init[-2:, :] = exact[-2:, :]
    init[:, :2] = exact[:, :2]
    init[:, -2:] = exact[:, -2:]
    solution_bih = pyfishpack.invert_StommelMunk(
        curl_bih,
        coords="cartesian",
        icbc=init,
        mParams={"A4": a4, "beta": beta, "R": resistance, "D": depth, "rho0": rho0},
        iParams={"mxLoop": 50000, "tolerance": 1e-12, "optArg": 1.55},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error_bih = float(np.max(np.abs(solution_bih - exact)))
    accepted_bih = _assert_close(
        "invert_StommelMunk A4!=0 manufactured solution", max_error_bih, 0.0, 1e-8
    )

    return {"a4_zero": accepted, "a4_nonzero": accepted_bih}


def check_high_level_stommel_arons(pyfishpack: Any) -> dict[str, Any]:
    m, n = 60, 56
    dy, dx = 0.85, 1.05
    epsilon = 1.4e-5
    f0 = 1.1e-4
    alpha = epsilon / (epsilon * epsilon + f0 * f0)
    i = np.arange(1, m + 1, dtype=np.float64)[:, None]
    j = np.arange(1, n + 1, dtype=np.float64)[None, :]
    exact = np.sin(math.pi * i / (m + 1.0)) * np.sin(math.pi * j / (n + 1.0))
    forcing = _discrete_constant_2d_fixed(exact, dy, dx, alpha, alpha)

    solution = pyfishpack.invert_StommelArons(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": 0.0, "epsilon": epsilon},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_StommelArons manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_stommel_arons_beta(pyfishpack: Any) -> dict[str, Any]:
    m, n = 36, 38
    dy, dx = 0.9, 1.05
    epsilon = 1.2
    f0 = 0.65
    beta = 0.025
    y = np.arange(m, dtype=np.float64)[:, None] * dy
    x = np.arange(n, dtype=np.float64)[None, :] * dx
    exact = np.sin(math.pi * y / ((m - 1.0) * dy)) * np.sin(
        math.pi * x / ((n - 1.0) * dx)
    )
    coriolis = f0 + beta * y
    denom = epsilon * epsilon + coriolis * coriolis
    c1 = epsilon / denom
    c1_dy = -2.0 * epsilon * coriolis * beta / (denom * denom)
    c2_dy = beta * (epsilon * epsilon - coriolis * coriolis) / (denom * denom)
    a = np.broadcast_to(c1, exact.shape)
    b = np.zeros_like(exact)
    c = np.broadcast_to(c1, exact.shape)
    d = np.broadcast_to(c1_dy, exact.shape)
    e = np.broadcast_to(-c2_dy, exact.shape)
    fcoef = np.zeros_like(exact)
    forcing = _discrete_general_2d_fixed(exact, dy, dx, a, b, c, d, e, fcoef)

    solution = pyfishpack.invert_StommelArons(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": beta, "epsilon": epsilon},
        iParams={"mxLoop": 30000, "tolerance": 1e-11},
        BCs=("fixed", "fixed"),
        spacing=(dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close(
        "invert_StommelArons beta-plane manufactured solution", max_error, 0.0, 1e-8
    )


def check_high_level_omega(pyfishpack: Any) -> dict[str, Any]:
    nz, ny, nx = 18, 20, 22
    dz, dy, dx = 1.4, 1.1, 0.9
    f0 = 1.2e-4
    n2 = 2.5e-4
    alpha_z = f0 * f0 / n2
    k = np.arange(1, nz + 1, dtype=np.float64)[:, None, None]
    j = np.arange(1, ny + 1, dtype=np.float64)[None, :, None]
    i = np.arange(1, nx + 1, dtype=np.float64)[None, None, :]
    exact = (
        np.sin(math.pi * k / (nz + 1.0))
        * np.sin(math.pi * j / (ny + 1.0))
        * np.sin(math.pi * i / (nx + 1.0))
    )
    forcing = n2 * _discrete_omega_operator_fixed(exact, dz, dy, dx, alpha_z)

    solution = pyfishpack.invert_omega(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": 0.0, "N2": n2},
        BCs=("fixed", "fixed", "fixed"),
        spacing=(dz, dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_omega manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_omega_beta(pyfishpack: Any) -> dict[str, Any]:
    nz, ny, nx = 12, 13, 14
    dz, dy, dx = 1.2, 0.9, 1.1
    f0 = 0.7
    beta = 0.02
    n2 = 1.1
    k = np.arange(nz, dtype=np.float64)[:, None, None]
    j = np.arange(ny, dtype=np.float64)[None, :, None]
    i = np.arange(nx, dtype=np.float64)[None, None, :]
    y = np.arange(ny, dtype=np.float64)[None, :, None] * dy
    acoef = (f0 + beta * y) ** 2
    exact = (
        np.sin(math.pi * k / (nz - 1.0))
        * np.sin(math.pi * j / (ny - 1.0))
        * np.sin(math.pi * i / (nx - 1.0))
    )
    forcing = np.zeros_like(exact)
    forcing[1:-1, 1:-1, 1:-1] = (
        acoef[:, 1:-1, :]
        * (exact[2:, 1:-1, 1:-1] - exact[1:-1, 1:-1, 1:-1])
        - acoef[:, 1:-1, :]
        * (exact[1:-1, 1:-1, 1:-1] - exact[:-2, 1:-1, 1:-1])
    ) / (dz * dz)
    forcing[1:-1, 1:-1, 1:-1] += n2 * (
        (exact[1:-1, 2:, 1:-1] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[1:-1, :-2, 1:-1])
        / (dy * dy)
        + (exact[1:-1, 1:-1, 2:] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[1:-1, 1:-1, :-2])
        / (dx * dx)
    )

    solution = pyfishpack.invert_omega(
        forcing,
        coords="cartesian",
        mParams={"f0": f0, "beta": beta, "N2": n2},
        iParams={"mxLoop": 50000, "tolerance": 1e-12},
        BCs=("fixed", "fixed", "fixed"),
        spacing=(dz, dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_omega beta-plane manufactured solution", max_error, 0.0, 1e-8)


def check_high_level_3d_ocean(pyfishpack: Any) -> dict[str, Any]:
    nz, ny, nx = 16, 18, 20
    dz, dy, dx = 1.25, 0.95, 1.05
    epsilon = 1.2e-5
    f0 = 1.1e-4
    n2 = 2.8e-4
    buoyancy_damping = 3.2e-5
    alpha_z = buoyancy_damping / n2
    alpha_horizontal = epsilon / (epsilon * epsilon + f0 * f0)
    k = np.arange(1, nz + 1, dtype=np.float64)[:, None, None]
    j = np.arange(1, ny + 1, dtype=np.float64)[None, :, None]
    i = np.arange(1, nx + 1, dtype=np.float64)[None, None, :]
    exact = (
        np.sin(math.pi * k / (nz + 1.0))
        * np.sin(math.pi * j / (ny + 1.0))
        * np.sin(math.pi * i / (nx + 1.0))
    )
    forcing = _discrete_constant_3d_fixed(
        exact, dz, dy, dx, alpha_z, alpha_horizontal, alpha_horizontal
    )

    solution = pyfishpack.invert_3DOcean(
        forcing,
        coords="cartesian",
        mParams={
            "f0": f0,
            "beta": 0.0,
            "epsilon": epsilon,
            "N2": n2,
            "k": buoyancy_damping,
        },
        BCs=("fixed", "fixed", "fixed"),
        spacing=(dz, dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_3DOcean manufactured solution", max_error, 0.0, 5e-12)


def check_high_level_3d_ocean_beta(pyfishpack: Any) -> dict[str, Any]:
    nz, ny, nx = 12, 13, 14
    dz, dy, dx = 1.2, 0.9, 1.1
    epsilon = 1.3
    f0 = 0.7
    beta = 0.02
    n2 = 1.1
    buoyancy_damping = 0.8
    vertical = buoyancy_damping / n2
    k = np.arange(nz, dtype=np.float64)[:, None, None]
    j = np.arange(ny, dtype=np.float64)[None, :, None]
    i = np.arange(nx, dtype=np.float64)[None, None, :]
    y = np.arange(ny, dtype=np.float64)[None, :, None] * dy
    coriolis = f0 + beta * y
    denom = epsilon * epsilon + coriolis * coriolis
    c1 = epsilon / denom
    c1_dy = -2.0 * epsilon * coriolis * beta / (denom * denom)
    c2_dy = beta * (epsilon * epsilon - coriolis * coriolis) / (denom * denom)
    exact = (
        np.sin(math.pi * k / (nz - 1.0))
        * np.sin(math.pi * j / (ny - 1.0))
        * np.sin(math.pi * i / (nx - 1.0))
    )
    forcing = np.zeros_like(exact)
    forcing[1:-1, 1:-1, 1:-1] = (
        vertical
        * (exact[2:, 1:-1, 1:-1] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[:-2, 1:-1, 1:-1])
        / (dz * dz)
        + c1[:, 1:-1, :]
        * (exact[1:-1, 2:, 1:-1] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[1:-1, :-2, 1:-1])
        / (dy * dy)
        + c1[:, 1:-1, :]
        * (exact[1:-1, 1:-1, 2:] - 2.0 * exact[1:-1, 1:-1, 1:-1] + exact[1:-1, 1:-1, :-2])
        / (dx * dx)
        + c1_dy[:, 1:-1, :]
        * (exact[1:-1, 2:, 1:-1] - exact[1:-1, :-2, 1:-1])
        / (2.0 * dy)
        - c2_dy[:, 1:-1, :]
        * (exact[1:-1, 1:-1, 2:] - exact[1:-1, 1:-1, :-2])
        / (2.0 * dx)
    )

    solution = pyfishpack.invert_3DOcean(
        forcing,
        coords="cartesian",
        mParams={
            "f0": f0,
            "beta": beta,
            "epsilon": epsilon,
            "N2": n2,
            "k": buoyancy_damping,
        },
        iParams={"mxLoop": 50000, "tolerance": 1e-12},
        BCs=("fixed", "fixed", "fixed"),
        spacing=(dz, dy, dx),
    )
    max_error = float(np.max(np.abs(solution - exact)))
    return _assert_close("invert_3DOcean beta-plane manufactured solution", max_error, 0.0, 1e-8)


def check_high_level_refstate(pyfishpack: Any) -> dict[str, Any]:
    import xarray as xr

    z = np.linspace(1.0, 2.0, 16)
    r = np.linspace(2.0, 3.0, 18)
    pv = xr.DataArray(np.full((16, 18), 2.0), coords={"z": z, "r": r}, dims=("z", "r"))
    solution = pyfishpack.invert_RefState(
        pv,
        ["z", "r"],
        coords="cartesian",
        mParams={"Gamma": 1.0, "ang0": 2.0e5},
        iParams={"mxLoop": 500, "tolerance": 1e-10, "printInfo": False},
    )
    if solution.name != "inverted" or solution.dims != pv.dims:
        raise AssertionError("invert_RefState did not preserve expected xarray shape metadata")
    max_abs = float(np.max(np.abs(solution.values)))
    if not np.isfinite(max_abs):
        raise AssertionError("invert_RefState returned non-finite values")
    return {"max_abs": max_abs, "shape": list(solution.shape), "backend": "sor_standard2d"}


def check_high_level_refstate_swm(pyfishpack: Any) -> dict[str, Any]:
    import xarray as xr

    lat = np.linspace(-60.0, 60.0, 31)
    q = xr.DataArray(np.full(31, 1.0e-6), coords={"lat": lat}, dims=("lat",))
    m0 = xr.DataArray(np.cos(np.deg2rad(lat)), coords={"lat": lat}, dims=("lat",))
    solution = pyfishpack.invert_RefStateSWM(
        q,
        ["lat"],
        mParams={"M0": m0, "C0": 1.0},
        iParams={"mxLoop": 500, "tolerance": 1e-10, "printInfo": False},
    )
    if solution.name != "inverted" or solution.dims != q.dims:
        raise AssertionError("invert_RefStateSWM did not preserve expected xarray shape metadata")
    max_abs = float(np.max(np.abs(solution.values)))
    if not np.isfinite(max_abs):
        raise AssertionError("invert_RefStateSWM returned non-finite values")
    return {"max_abs": max_abs, "shape": list(solution.shape), "backend": "sor_standard1d"}


def check_genbun_size_sweep(pyfishpack: Any) -> dict[str, Any]:
    fishpack = pyfishpack.fishpack
    sizes = [
        88,
        89,
        90,
        91,
        92,
        93,
        94,
        95,
        96,
        97,
        98,
        99,
        100,
        101,
        102,
        103,
        104,
        105,
        106,
        107,
        108,
        109,
        110,
        111,
        112,
        127,
        128,
        129,
        255,
        256,
        257,
    ]
    for size in sizes:
        a = np.ones(size, dtype=np.float64)
        b = np.full(size, -2.0, dtype=np.float64)
        c = np.ones(size, dtype=np.float64)
        a[0] = 0.0
        c[-1] = 0.0
        y = np.zeros((size, size), dtype=np.float64, order="F")
        y[size // 2, size // 2] = 1.0
        _, ierror = fishpack.genbun(1, size, 1, size, a, b, c, y)
        if ierror != 0:
            raise AssertionError(f"genbun size {size} returned ierror={ierror}")
    return {"sizes": sizes, "count": len(sizes)}


def _time_median(func: Any, repeats: int) -> float:
    timings = []
    for _ in range(repeats):
        start = time.perf_counter()
        func()
        timings.append(time.perf_counter() - start)
    return statistics.median(timings)


def benchmark_xinvert(
    pyfishpack: Any, xinvert_path: Path, sizes: list[int], repeats: int
) -> dict[str, Any]:
    if not xinvert_path.exists():
        return {"available": False, "reason": f"not found: {xinvert_path}"}

    if str(xinvert_path) not in sys.path:
        sys.path.insert(0, str(xinvert_path))

    try:
        import xarray as xr
        from xinvert.apps import invert_Poisson as xinvert_poisson
    except Exception as exc:  # pragma: no cover - environment dependent
        return {"available": False, "reason": f"import failed: {exc}"}

    quiet_iparams = {
        "BCs": ["fixed", "fixed"],
        "undef": np.nan,
        "mxLoop": 20000,
        "tolerance": 1e-8,
        "optArg": None,
        "printInfo": False,
        "debug": False,
    }

    results = []
    for size in sizes:
        y = np.arange(size, dtype=np.float64)
        x = np.arange(size, dtype=np.float64)
        ii = np.arange(1, size + 1, dtype=np.float64)[:, None]
        jj = np.arange(1, size + 1, dtype=np.float64)[None, :]
        exact = np.sin(math.pi * ii / (size + 1.0)) * np.sin(
            math.pi * jj / (size + 1.0)
        )
        rhs = np.zeros_like(exact)
        rhs += -4.0 * exact
        rhs[1:, :] += exact[:-1, :]
        rhs[:-1, :] += exact[1:, :]
        rhs[:, 1:] += exact[:, :-1]
        rhs[:, :-1] += exact[:, 1:]
        field = xr.DataArray(rhs, coords={"y": y, "x": x}, dims=("y", "x"))

        def run_pyfishpack() -> Any:
            return pyfishpack.invert_Poisson(
                field, dims=("y", "x"), coords="cartesian", iParams=quiet_iparams
            )

        def run_xinvert() -> Any:
            return xinvert_poisson(
                field,
                dims=["y", "x"],
                coords="cartesian",
                mParams={},
                iParams=quiet_iparams,
            )

        run_pyfishpack()
        run_xinvert()
        t_pyfishpack = _time_median(run_pyfishpack, repeats)
        t_xinvert = _time_median(run_xinvert, repeats)
        results.append(
            {
                "size": size,
                "pyfishpack_seconds": t_pyfishpack,
                "xinvert_seconds": t_xinvert,
                "speedup": t_xinvert / t_pyfishpack,
            }
        )

    return {"available": True, "xinvert_path": str(xinvert_path), "results": results}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xinvert-path",
        type=Path,
        default=DEFAULT_XINVERT_PATH,
        help="Local xinvert checkout used for optional speed comparison.",
    )
    parser.add_argument(
        "--skip-xinvert",
        action="store_true",
        help="Skip optional xinvert speed comparison.",
    )
    parser.add_argument(
        "--benchmark-sizes",
        default="96,128,256",
        help="Comma-separated square grid sizes for optional xinvert benchmark.",
    )
    parser.add_argument(
        "--benchmark-repeats",
        type=int,
        default=5,
        help="Median repeats per benchmark size.",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help="Optional JSON output path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    pyfishpack = _import_pyfishpack()
    sizes = [int(item) for item in args.benchmark_sizes.split(",") if item.strip()]

    report: dict[str, Any] = {
        "backend_info": check_backend_info(pyfishpack),
        "backend_methods": check_backend_methods(pyfishpack),
        "checks": {
            "sor_standard1d_residual": check_sor_standard1d_backend(pyfishpack),
            "sor_standard2d_residual": check_sor_standard2d_backend(pyfishpack),
            "sor_standard3d_manufactured": check_sor_standard3d_backend(pyfishpack),
            "sor_general2d_manufactured": check_sor_general2d_backend(pyfishpack),
            "sor_general3d_manufactured": check_sor_general3d_backend(pyfishpack),
            "sor_biharmonic2d_manufactured": check_sor_biharmonic2d_backend(pyfishpack),
            "fftpack_transform_roundtrip": check_fftpack_transforms(pyfishpack),
            "genbun_ucar_example": check_genbun_example(pyfishpack),
            "poistg_ucar_example": check_poistg_example(pyfishpack),
            "pois3d_ucar_example": check_pois3d_example(pyfishpack),
            "hwscrt_ucar_example": check_hwscrt_example(pyfishpack),
            "hstcrt_ucar_example": check_hstcrt_example(pyfishpack),
            "hw3crt_ucar_example": check_hw3crt_example(pyfishpack),
            "hwsplr_ucar_example": check_hwsplr_example(pyfishpack),
            "hwscyl_ucar_example": check_hwscyl_example(pyfishpack),
            "hstplr_ucar_example": check_hstplr_example(pyfishpack),
            "hstcyl_ucar_example": check_hstcyl_example(pyfishpack),
            "hwsssp_ucar_example": check_hwsssp_example(pyfishpack),
            "hstssp_ucar_example": check_hstssp_example(pyfishpack),
            "hwscsp_ucar_example": check_hwscsp_example(pyfishpack),
            "hstcsp_ucar_example": check_hstcsp_example(pyfishpack),
            "invert_poisson_manufactured": check_high_level_poisson(pyfishpack),
            "invert_multigrid_poisson_direct": check_high_level_multigrid(pyfishpack),
            "invert_geostrophic_manufactured": check_high_level_geostrophic(pyfishpack),
            "invert_geostrophic_beta_fortran_sor": check_high_level_geostrophic_beta(pyfishpack),
            "invert_pv2d_manufactured": check_high_level_pv2d(pyfishpack),
            "invert_eliassen_manufactured": check_high_level_eliassen(pyfishpack),
            "invert_eliassen_cross_fortran_sor": check_high_level_eliassen_cross(pyfishpack),
            "invert_bretherton_haidvogel_manufactured": check_high_level_bretherton(pyfishpack),
            "invert_bretherton_haidvogel_beta_fortran": check_high_level_bretherton_beta(pyfishpack),
            "invert_fofonoff_zero_forcing": check_high_level_fofonoff(pyfishpack),
            "invert_fofonoff_beta_fortran": check_high_level_fofonoff_beta(pyfishpack),
            "invert_gill_matsuno_manufactured": check_high_level_gill_matsuno(pyfishpack),
            "invert_gill_matsuno_beta_fortran_sor": check_high_level_gill_matsuno_beta(pyfishpack),
            "invert_gill_matsuno_test_manufactured": check_high_level_gill_matsuno_test(pyfishpack),
            "invert_stommel_manufactured": check_high_level_stommel(pyfishpack),
            "invert_stommel_beta_fortran_sor": check_high_level_stommel_beta(pyfishpack),
            "invert_stommel_test_manufactured": check_high_level_stommel_test(pyfishpack),
            "invert_stommel_munk_a4_zero": check_high_level_stommel_munk(pyfishpack),
            "invert_stommel_arons_manufactured": check_high_level_stommel_arons(pyfishpack),
            "invert_stommel_arons_beta_fortran_sor": check_high_level_stommel_arons_beta(pyfishpack),
            "invert_omega_manufactured": check_high_level_omega(pyfishpack),
            "invert_omega_beta_fortran_sor": check_high_level_omega_beta(pyfishpack),
            "invert_3d_ocean_manufactured": check_high_level_3d_ocean(pyfishpack),
            "invert_3d_ocean_beta_fortran_sor": check_high_level_3d_ocean_beta(pyfishpack),
            "invert_refstate_fortran_sor": check_high_level_refstate(pyfishpack),
            "invert_refstate_swm_fortran_sor": check_high_level_refstate_swm(pyfishpack),
            "genbun_size_sweep": check_genbun_size_sweep(pyfishpack),
        },
    }

    if args.skip_xinvert:
        report["xinvert_benchmark"] = {"available": False, "reason": "skipped"}
    else:
        report["xinvert_benchmark"] = benchmark_xinvert(
            pyfishpack, args.xinvert_path, sizes, args.benchmark_repeats
        )

    text = json.dumps(report, indent=2, sort_keys=True)
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
