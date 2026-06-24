# PyFishPack

PyFishPack is a modernized Python-facing Fortran solver package derived from NCAR FISHPACK.

Project links:

- Homepage: https://github.com/QianyeSu/PyFishPack
- Repository: https://github.com/QianyeSu/PyFishPack
- Bug Tracker: https://github.com/QianyeSu/PyFishPack/issues

## Source layout

- Modern Fortran `.f90` sources live under `PyFishPack/src`.
- Original F77 `.f` sources are kept as archive material under `PyFishPack/src/archive/f77`.
- The backend uses modern Fortran `bind(C)` entry points plus a thin NumPy C API wrapper.

## Build

Use the `skyborn_dev` conda environment with Meson and Ninja.

```bash
conda activate skyborn_dev
meson setup build
ninja -C build
```

## Public APIs

PyFishPack exposes xinvert-style equation helpers from `PyFishPack.apps`, including:

`invert_Poisson`, `invert_geostrophic`, `invert_PV2D`, `invert_Eliassen`, `invert_Fofonoff`, `invert_GillMatsuno`, `invert_GillMatsuno_test`, `invert_BrethertonHaidvogel`, `invert_Stommel`, `invert_StommelMunk`, `invert_Stommel_test`, `invert_StommelArons`, `invert_omega`, `invert_3DOcean`, `invert_RefState`, `invert_RefStateSWM`, `invert_MultiGrid`.

The compiled `fishpack` backend also exposes direct solver and transform entry points for the modern Fortran kernels, including the `genbun`, `poistg`, `pois3d`, Helmholtz, SOR, and FFTPACK interfaces.

## Verification

The current backend verification covers the exposed equation APIs, source-layout checks, a local xinvert benchmark, and backend import checks. The modern backend reports `modern-fortran` through `fishpack.backend_info()`.
