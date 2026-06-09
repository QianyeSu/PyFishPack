"""Audit PyFishPack modernization and xinvert-style API coverage."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_XINVERT_PATH = Path(r"I:\test_xinvert\xinvert-master")


F77_TO_F90 = {
    "blktri": "real_block_tridiagonal_linear_systems_solver.f90",
    "cblktri": "complex_block_tridiagonal_linear_systems_solver.f90",
    "cmgnbn": "complex_linear_systems_solver.f90",
    "comf": "type_GeneralizedCyclicReductionUtility.f90",
    "fftpack": "type_PeriodicFastFourierTransform.f90",
    "genbun": "centered_real_linear_systems_solver.f90",
    "gnbnaux": "type_CyclicReductionUtility.f90",
    "hstcrt": "staggered_cartesian_solver.f90",
    "hstcsp": "staggered_axisymmetric_spherical_solver.f90",
    "hstcyl": "staggered_cylindrical_solver.f90",
    "hstplr": "staggered_polar_solver.f90",
    "hstssp": "staggered_spherical_solver.f90",
    "hw3crt": "centered_cartesian_helmholtz_solver_3d.f90",
    "hwscrt": "centered_cartesian_solver.f90",
    "hwscsp": "centered_axisymmetric_spherical_solver.f90",
    "hwscyl": "centered_cylindrical_solver.f90",
    "hwsplr": "centered_polar_solver.f90",
    "hwsssp": "centered_spherical_solver.f90",
    "pois3d": "general_linear_systems_solver_3d.f90",
    "poistg": "staggered_real_linear_systems_solver.f90",
    "sepaux": "type_SepAux.f90",
    "sepeli": "sepeli.f90",
    "sepx4": "sepx4.f90",
}


FORTRAN_BACKED_EQUATION_APIS = {
    "invert_3DOcean": "Cartesian, uniform-grid, constant-coefficient beta=0 subset backed by pois3d.",
    "invert_BrethertonHaidvogel": "Cartesian, constant-depth Helmholtz subset backed by genbun; beta-plane forcing is supported for xarray coordinates and NumPy spacing.",
    "invert_Eliassen": "Cartesian constant-coefficient B=0 subset backed by genbun; constant cross-derivative B!=0 subset backed by the modern Fortran SOR general2d kernel.",
    "invert_Fofonoff": "Cartesian Helmholtz subset backed by genbun; beta-plane forcing is supported for xarray coordinates and NumPy spacing.",
    "invert_GillMatsuno": "Cartesian beta=0 subset backed by genbun; Cartesian beta-plane subset backed by the modern Fortran SOR general2d kernel.",
    "invert_GillMatsuno_test": "Cartesian, constant-Coriolis beta=0 test-form subset backed by genbun.",
    "invert_PV2D": "Cartesian, scalar-N2 constant-coefficient QG PV subset backed by genbun.",
    "invert_Poisson": "Cartesian, uniform-grid Poisson equation backed by genbun.",
    "invert_RefState": "Xinvert-style variable-coefficient reference-state equation backed by the modern Fortran SOR standard2d kernel.",
    "invert_RefStateSWM": "Xinvert-style one-dimensional shallow-water reference-state equation backed by the modern Fortran SOR standard1d kernel.",
    "invert_Stommel": "Cartesian beta=0 subset backed by genbun; Cartesian beta-plane subset backed by the modern Fortran SOR general2d kernel.",
    "invert_StommelArons": "Cartesian beta=0 subset backed by genbun; Cartesian beta-plane subset backed by the modern Fortran SOR general2d kernel.",
    "invert_StommelMunk": "Cartesian A4=0 subset that delegates to Stommel; A4!=0 biharmonic cases remain unsupported.",
    "invert_Stommel_test": "Cartesian test-form subset delegated to Stommel; beta=0 uses genbun and beta-plane uses the modern Fortran SOR general2d kernel.",
    "invert_geostrophic": "Cartesian constant-Coriolis beta=0 subset backed by genbun; Cartesian beta-plane subset backed by the modern Fortran SOR standard2d kernel.",
    "invert_omega": "Cartesian, scalar-N2, beta=0 3D QG omega subset backed by pois3d.",
}


COMPATIBILITY_HELPER_APIS = {
    "invert_MultiGrid": (
        "Compatibility helper for xinvert's multi-grid entry point. PyFishPack "
        "uses direct Fishpack solvers, so this wrapper delegates once to the "
        "provided Fortran-backed inversion function and returns xinvert-style "
        "(solution, grids, history) output."
    ),
}


UNSUPPORTED_XINVERT_EQUATION_APIS: dict[str, str] = {}


def _read_makefile_routines() -> list[str]:
    makefile = REPO_ROOT / "PyFishPack" / "src" / "archive" / "f77" / "Makefile"
    text = makefile.read_text(encoding="utf-8")
    match = re.search(r"SRC=(.*?)(?:\n\n|\Z)", text, flags=re.S)
    if not match:
        return sorted(F77_TO_F90)
    src_text = match.group(1).replace("\\\n", " ")
    return sorted(Path(item).stem for item in src_text.split() if item.endswith(".f"))


def _modern_source_coverage() -> dict[str, Any]:
    src_dir = REPO_ROOT / "PyFishPack" / "src"
    routines = _read_makefile_routines()
    entries = []
    for routine in routines:
        f90_name = F77_TO_F90.get(routine)
        f90_path = src_dir / f90_name if f90_name else None
        entries.append(
            {
                "routine": routine,
                "f77_source": f"PyFishPack/src/archive/f77/{routine}.f",
                "modern_f90_source": f"PyFishPack/src/{f90_name}" if f90_name else None,
                "modern_f90_exists": bool(f90_path and f90_path.exists()),
            }
        )
    return {
        "total_f77_routines": len(entries),
        "modern_f90_present": sum(1 for item in entries if item["modern_f90_exists"]),
        "entries": entries,
    }


def _source_layout() -> dict[str, Any]:
    src_dir = REPO_ROOT / "PyFishPack" / "src"
    archive_dir = src_dir / "archive" / "f77"
    f90_subdir = src_dir / "F90"
    return {
        "root_f90_sources": sorted(path.name for path in src_dir.glob("*.f90")),
        "root_f77_sources": sorted(path.name for path in src_dir.glob("*.f")),
        "archived_f77_sources": sorted(path.name for path in archive_dir.glob("*.f")),
        "stale_f90_subdir_sources": (
            sorted(path.name for path in f90_subdir.glob("*.f90"))
            if f90_subdir.exists()
            else []
        ),
    }


def _python_backend_exports() -> dict[str, Any]:
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    import PyFishPack as pyfishpack

    backend_methods = sorted(
        name for name in dir(pyfishpack.fishpack) if not name.startswith("_")
    )
    equation_apis = sorted(
        name for name in dir(pyfishpack) if name.startswith("invert_")
    )
    return {
        "backend_info": dict(pyfishpack.fishpack.backend_info()),
        "backend_methods": backend_methods,
        "equation_apis": equation_apis,
    }


def _xinvert_equation_apis(xinvert_path: Path) -> dict[str, Any]:
    apps = xinvert_path / "xinvert" / "apps.py"
    if not apps.exists():
        return {"available": False, "reason": f"not found: {apps}"}
    text = apps.read_text(encoding="utf-8")
    names = sorted(set(re.findall(r"^def (invert_[A-Za-z0-9_]+)\(", text, flags=re.M)))
    return {"available": True, "xinvert_path": str(xinvert_path), "equation_apis": names}


def _equation_gap_analysis(exports: dict[str, Any], reference: dict[str, Any]) -> dict[str, Any]:
    exported = set(exports["equation_apis"])
    reference_names = set(reference.get("equation_apis", []))
    implemented = sorted(exported & set(FORTRAN_BACKED_EQUATION_APIS))
    helpers = sorted(exported & set(COMPATIBILITY_HELPER_APIS))
    unsupported = {
        name: UNSUPPORTED_XINVERT_EQUATION_APIS[name]
        for name in sorted(reference_names - exported)
        if name in UNSUPPORTED_XINVERT_EQUATION_APIS
    }
    unexpected_missing = sorted(
        name for name in reference_names - exported if name not in unsupported
    )
    unsupported_but_exported = sorted(
        name for name in exported if name in UNSUPPORTED_XINVERT_EQUATION_APIS
    )
    return {
        "fortran_backed_implemented": {
            name: FORTRAN_BACKED_EQUATION_APIS[name] for name in implemented
        },
        "compatibility_helpers": {
            name: COMPATIBILITY_HELPER_APIS[name] for name in helpers
        },
        "unsupported_reference_apis": unsupported,
        "unexpected_missing_reference_apis": unexpected_missing,
        "unsupported_but_exported": unsupported_but_exported,
        "xinvert_reference_count": len(reference_names),
        "fortran_backed_count": len(implemented),
        "compatibility_helper_count": len(helpers),
    }


def audit(xinvert_path: Path) -> dict[str, Any]:
    exports = _python_backend_exports()
    reference = _xinvert_equation_apis(xinvert_path)
    return {
        "source_layout": _source_layout(),
        "modern_source_coverage": _modern_source_coverage(),
        "python_backend_exports": exports,
        "xinvert_reference": reference,
        "equation_gap_analysis": _equation_gap_analysis(exports, reference),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xinvert-path", type=Path, default=DEFAULT_XINVERT_PATH)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()

    report = audit(args.xinvert_path)
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
