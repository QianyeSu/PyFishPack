# PyFishPack - Python extension backed by modern Fortran and a C wrapper.
project(
    'pyfishpack',
    'c',
    'fortran',
    version: '0.1.0',
    default_options: ['warning_level=2', 'buildtype=release', 'b_lto=true'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()
build_py = find_program('build_python', required: false)
if not build_py.found()
    build_py = py
endif

py.install_sources(
    ['__init__.py', 'apps.py'],
    subdir: 'PyFishPack',
)

incdir_numpy_cmd = run_command(
    build_py,
    ['-c', '''
import numpy
import os

paths = [
    numpy.get_include(),
    os.path.join(numpy.__path__[0], "core", "include"),
    os.path.join(numpy.__path__[0], "_core", "include"),
]
for path in paths:
    if os.path.exists(os.path.join(path, "numpy", "npy_os.h")):
        print(path)
        break
else:
    print(numpy.get_include())
'''],
    check: true,
)
incdir_numpy = incdir_numpy_cmd.stdout().strip()
message('NumPy include directory: ' + incdir_numpy)

host_system = host_machine.system()
host_cpu = host_machine.cpu_family()
message('Host system: ' + host_system)
message('Host CPU family: ' + host_cpu)

fortran_compiler = meson.get_compiler('fortran')
c_compiler = meson.get_compiler('c')
is_gcc_toolchain = fortran_compiler.get_id() == 'gcc' and c_compiler.get_id() == 'gcc'

openmp_dependencies = []
openmp_fortran_args = []
openmp_c_args = []
openmp_link_args = []
if is_gcc_toolchain
    openmp_dep = dependency('openmp', required: false)
    if openmp_dep.found()
        openmp_dependencies += openmp_dep
        openmp_fortran_args += ['-fopenmp']
        openmp_c_args += ['-fopenmp']
        openmp_link_args += ['-fopenmp']
        message('Using OpenMP dependency for parallel processing')
    endif
elif host_system != 'windows'
    openmp_dep = dependency('openmp', required: false)
    if openmp_dep.found()
        openmp_dependencies += openmp_dep
        message('Using OpenMP dependency for parallel processing')
    endif
endif

fortran_args = [
    '-O3',
    '-fPIC',
    '-funroll-loops',
    '-finline-functions',
    '-ftree-vectorize',
    '-ffinite-math-only',
    '-ffree-line-length-none',
    '-fno-common',
    '-ftree-loop-im',
    '-ftree-loop-distribution',
    '-falign-functions=32',
    '-falign-loops=32',
    '-fimplicit-none',
]

c_args = [
    '-DNPY_NO_DEPRECATED_API=NPY_1_20_API_VERSION',
    '-fPIC',
    '-fno-trapping-math',
    '-falign-functions=32',
]

if host_system == 'darwin' and host_cpu == 'aarch64'
    message('Configuring for Apple Silicon (arm64)')
    fortran_args += ['-march=armv8-a', '-mtune=apple-m1']
    c_args += ['-march=armv8-a', '-mtune=apple-m1']
else
    message('Configuring for x86-64 architecture')
    fortran_args += ['-march=x86-64', '-mtune=generic']
    c_args += ['-march=x86-64', '-mtune=generic', '-O3', '-funroll-loops', '-finline-functions', '-ftree-vectorize']
endif

fortran_args += openmp_fortran_args
c_args += openmp_c_args

message('Fortran args: ' + ' '.join(fortran_args))
message('C args: ' + ' '.join(c_args))

inc_np = include_directories(incdir_numpy)
dependencies = [py_dep] + openmp_dependencies
package_install_dir = get_option('python.purelibdir') / 'PyFishPack'

modern_fortran_sources = [
    'src/fishpack_precision.f90',
    'src/type_FishpackWorkspace.f90',
    'src/type_GeneralizedCyclicReductionUtility.f90',
    'src/type_CyclicReductionUtility.f90',
    'src/type_PeriodicFastFourierTransform.f90',
    'src/type_SepAux.f90',
    'src/type_CenteredCyclicReductionUtility.f90',
    'src/type_StaggeredCyclicReductionUtility.f90',
    'src/real_block_tridiagonal_linear_systems_solver.f90',
    'src/complex_block_tridiagonal_linear_systems_solver.f90',
    'src/complex_linear_systems_solver.f90',
    'src/three_dimensional_solvers.f90',
    'src/general_linear_systems_solver_3d.f90',
    'src/centered_cartesian_helmholtz_solver_3d.f90',
    'src/staggered_real_linear_systems_solver.f90',
    'src/centered_real_linear_systems_solver.f90',
    'src/staggered_cartesian_solver.f90',
    'src/staggered_polar_solver.f90',
    'src/staggered_cylindrical_solver.f90',
    'src/staggered_spherical_solver.f90',
    'src/staggered_axisymmetric_spherical_solver.f90',
    'src/staggered_helmholtz_solvers.f90',
    'src/centered_cartesian_solver.f90',
    'src/centered_polar_solver.f90',
    'src/centered_cylindrical_solver.f90',
    'src/centered_spherical_solver.f90',
    'src/centered_axisymmetric_spherical_solver.f90',
    'src/centered_helmholtz_solvers.f90',
    'src/sepeli.f90',
    'src/sepx4.f90',
    'src/iterative_solvers.f90',
    'src/fftpack_c_api.f90',
    'src/fishpack_c_api.f90',
]

fishpack_ext = py.extension_module(
    'fishpack',
    ['src/pyfishpack_module.c'] + modern_fortran_sources,
    include_directories: inc_np,
    dependencies: dependencies,
    fortran_args: fortran_args,
    c_args: c_args,
    link_args: openmp_link_args,
    install: true,
    install_dir: package_install_dir,
    build_by_default: true,
)

custom_target(
    'copy_to_source_for_inplace',
    input: fishpack_ext,
    output: 'inplace_copy_marker.txt',
    command: [
        py,
        '-c',
        '''
import shutil
import sys
from pathlib import Path

current_dir = Path.cwd()
module_dir = current_dir.parent
patterns = ["fishpack*.so", "fishpack*.pyd", "fishpack*.dylib"]

for stale_pattern in patterns:
    for stale_file in module_dir.glob(stale_pattern):
        try:
            stale_file.unlink()
            print(f"Removed stale source artifact {stale_file.name}")
        except Exception as exc:
            print(f"Failed to remove stale source artifact {stale_file.name}: {exc}")

copied = False
for pattern in patterns:
    for built_file in current_dir.rglob(pattern):
        target = module_dir / built_file.name
        if not target.exists() or built_file.stat().st_mtime > target.stat().st_mtime:
            shutil.copy2(str(built_file), str(target))
            print(f"Copied {built_file.name} to source")
            copied = True

if not copied:
    print("No fishpack extension artifact needed copying")

with open(sys.argv[1], "w", encoding="utf-8") as handle:
    handle.write("PyFishPack extension copied\\n")
''',
        '@OUTPUT@',
    ],
    build_by_default: true,
    console: true,
)

message('PyFishPack configured with modern Fortran backend and NumPy C API wrapper')
