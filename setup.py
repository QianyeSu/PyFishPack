"""
Setup script for PyFishPack - Python wrapper for FISHPACK Fortran library
Mixed build system with meson for Fortran modules
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from setuptools.command.install import install
import numpy as np

# Check if we're in documentation build mode
DOCS_BUILD_MODE = (
    os.environ.get("PYFISHPACK_DOCS_BUILD") == "1" or os.environ.get(
        "SKIP_FORTRAN") == "1"
)

if DOCS_BUILD_MODE:
    print("Documentation build mode detected - skipping Fortran compilation")
else:
    # Force gfortran compiler usage
    os.environ["FC"] = os.environ.get("FC", "gfortran")
    os.environ["F77"] = os.environ.get("F77", "gfortran")
    os.environ["F90"] = os.environ.get("F90", "gfortran")
    os.environ["CC"] = os.environ.get("CC", "gcc")


# Check if gfortran is available
def check_gfortran():
    """Check if gfortran is available"""
    try:
        result = subprocess.run(
            ["gfortran", "--version"], capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            print(
                f"Found gfortran: {result.stdout.split()[4] if len(result.stdout.split()) > 4 else 'unknown version'}"
            )
            return True
    except (subprocess.TimeoutExpired, subprocess.SubprocessError, FileNotFoundError):
        pass

    print("Warning: gfortran not found. Fortran extensions may not build correctly.")
    print("Please install gfortran:")
    print("  Linux: sudo apt-get install gfortran")
    print("  macOS: brew install gcc")
    print("  Windows: conda install m2w64-toolchain")
    return False


# Check gfortran availability at setup time (skip in docs mode)
if not DOCS_BUILD_MODE:
    check_gfortran()
else:
    print("Skipping gfortran check in documentation build mode")


class MesonBuildExt(build_ext):
    """Custom build extension to handle meson builds for Fortran modules"""

    def run(self):
        """Run the build process"""
        print("DEBUG: MesonBuildExt.run() called")
        # Build meson modules first
        self.build_meson_modules()
        # Then run the standard build_ext
        super().run()

    def build_meson_modules(self):
        """Build modules that use meson (like fishpack)"""
        print("DEBUG: build_meson_modules() called")

        # Skip meson builds in documentation mode
        if DOCS_BUILD_MODE:
            print("📚 Documentation build mode - skipping meson module builds")
            return

        # Determine target directory based on --inplace flag
        if self.inplace:
            print("DEBUG: --inplace detected, building to source directory")
            fishpack_target = Path("PyFishPack") / "src"
        else:
            print("DEBUG: Building to build directory")
            fishpack_target = Path(self.build_lib) / "PyFishPack"

        # Auto-discover meson modules based on directory structure
        # Each module should have a meson.build file
        meson_modules = self._discover_meson_modules()

        for module in meson_modules:
            print(f"DEBUG: Processing module {module['name']}")
            if self.should_build_meson_module(module):
                print(f"DEBUG: Building module {module['name']} with meson")
                self.build_meson_module(module)
            else:
                print(
                    f"DEBUG: Skipping module {module['name']} - no meson.build found")

    def should_build_meson_module(self, module):
        """Check if we should build this meson module"""
        meson_build_file = module["path"] / "meson.build"
        return meson_build_file.exists()

    def check_meson_available(self):
        """Check if meson and ninja are available"""
        try:
            # Check meson
            result = subprocess.run(
                ["meson", "--version"], capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                return False, "meson not found"

            meson_version = result.stdout.strip()
            print(f"Found meson version: {meson_version}")

            # Check ninja
            result = subprocess.run(
                ["ninja", "--version"], capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                return False, "ninja not found"

            ninja_version = result.stdout.strip()
            print(f"Found ninja version: {ninja_version}")

            return True, None

        except (
            subprocess.TimeoutExpired,
            subprocess.SubprocessError,
            FileNotFoundError,
        ) as e:
            return False, str(e)

    def build_meson_module(self, module):
        """
        Build a meson module using the meson build system.
        """
        print(f"Building {module['name']} with meson build system...")

        # Check if meson and ninja are available
        meson_available, error_msg = self.check_meson_available()
        if not meson_available:
            print(f"ERROR: Meson build tools not available: {error_msg}")
            print("Please install meson and ninja:")
            print("  pip install meson ninja")
            print("  or: conda install meson ninja")
            raise RuntimeError(
                f"Meson build tools required but not available: {error_msg}"
            )

        module_path = module["path"]
        # Use build subdirectory as specified in requirements
        build_dir = module_path / "build"

        try:
            # Clean build directory
            if build_dir.exists():
                print(f"Cleaning existing build directory: {build_dir}")
                shutil.rmtree(build_dir)

            # Setup build directory
            build_dir.mkdir(parents=True, exist_ok=True)

            # Configure meson build with custom install directory for wheel builds
            print(
                f"Configuring meson build in {build_dir} (cwd={module_path})")

            setup_cmd = [
                "meson",
                "setup",
                "build",  # build directory inside module_path
                ".",  # source is current directory (module_path)
                "--buildtype=release",
                "-Db_lto=true",
            ]

            # For wheel builds, configure custom install directory
            if not self.inplace and hasattr(self, "build_lib") and self.build_lib:
                # Tell meson to install to our build directory instead of system
                build_lib_path = Path(self.build_lib).resolve()
                setup_cmd.extend(
                    [
                        f"--python.purelibdir={build_lib_path}",
                        f"--python.platlibdir={build_lib_path}",
                    ]
                )
                print(
                    f"DEBUG: Configuring meson to install to: {build_lib_path}")

            print(f"Running: {' '.join(setup_cmd)} (cwd={module_path})")

            # Set up environment for conda gfortran across all platforms
            env = os.environ.copy()
            import platform

            conda_prefix = env.get("CONDA_PREFIX", "")
            if conda_prefix:
                system = platform.system()
                current_path = env.get("PATH", "")

                if system == "Windows":
                    # Windows conda environment setup
                    conda_bin = os.path.join(conda_prefix, "bin")
                    conda_library_bin = os.path.join(
                        conda_prefix, "Library", "bin")
                    env["PATH"] = f"{conda_bin};{conda_library_bin};{current_path}"
                    print(
                        f"Enhanced PATH for Windows conda environment: {conda_prefix}"
                    )

                elif system in ["Linux", "Darwin"]:
                    # Linux and macOS conda environment setup
                    conda_bin = os.path.join(conda_prefix, "bin")
                    # On Unix-like systems, use colon separator and prepend to PATH
                    env["PATH"] = f"{conda_bin}:{current_path}"

                    # Add lib directory to LD_LIBRARY_PATH (Linux) or DYLD_LIBRARY_PATH (macOS)
                    conda_lib = os.path.join(conda_prefix, "lib")
                    if system == "Linux":
                        current_lib_path = env.get("LD_LIBRARY_PATH", "")
                        env["LD_LIBRARY_PATH"] = (
                            f"{conda_lib}:{current_lib_path}"
                            if current_lib_path
                            else conda_lib
                        )
                        print(
                            f"Enhanced PATH and LD_LIBRARY_PATH for Linux conda environment: {conda_prefix}"
                        )
                    else:  # macOS
                        current_lib_path = env.get("DYLD_LIBRARY_PATH", "")
                        env["DYLD_LIBRARY_PATH"] = (
                            f"{conda_lib}:{current_lib_path}"
                            if current_lib_path
                            else conda_lib
                        )
                        print(
                            f"Enhanced PATH and DYLD_LIBRARY_PATH for macOS conda environment: {conda_prefix}"
                        )

            subprocess.run(setup_cmd, cwd=str(
                module_path), check=True, env=env)

            # Build with ninja
            print(f"Building with ninja in {build_dir} (cwd={module_path})")
            build_cmd = ["ninja", "-C", "build"]

            print(f"Running: {' '.join(build_cmd)} (cwd={module_path})")
            result = subprocess.run(
                build_cmd,
                cwd=str(module_path),
                check=True,
                capture_output=True,
                text=True,
            )

            if result.stdout:
                print("Build output:", result.stdout)
            if result.stderr:
                print("Build warnings/errors:", result.stderr)

            # Install using meson
            if not self.inplace and hasattr(self, "build_lib") and self.build_lib:
                print(f"Installing meson build outputs to {self.build_lib}")
                install_cmd = ["meson", "install",
                               "-C", "build", "--only-changed"]
                print(f"Running: {' '.join(install_cmd)} (cwd={module_path})")

                try:
                    install_result = subprocess.run(
                        install_cmd,
                        cwd=str(module_path),
                        check=True,
                        capture_output=True,
                        text=True,
                        env=env,
                    )
                    if install_result.stdout:
                        print("Install output:", install_result.stdout)
                except subprocess.CalledProcessError as e:
                    print(f"ERROR: Meson install failed: {e}")
                    print(
                        "This should not happen with the new setup. Please check meson configuration."
                    )
                    raise
            else:
                print("Inplace build - extensions handled by meson custom_target")

            print(f"Meson build for {module['name']} completed successfully!")

            self._built_modules = getattr(self, "_built_modules", set())
            self._built_modules.add(module["name"])

        except (subprocess.CalledProcessError, RuntimeError, FileNotFoundError) as e:
            print(f"ERROR: Meson build failed for {module['name']}: {e}")
            if isinstance(e, subprocess.CalledProcessError):
                print(f"Command failed with exit code: {e.returncode}")
                if hasattr(e, "stdout") and e.stdout:
                    print("Stdout:", e.stdout)
                if hasattr(e, "stderr") and e.stderr:
                    print("Stderr:", e.stderr)
            raise  # Re-raise the exception

    def _discover_meson_modules(self):
        """
        Auto-discover meson modules by looking for meson.build files
        in PyFishPack subpackages.
        """
        modules = []
        pyfishpack_root = Path("PyFishPack")

        # Look for meson.build file in PyFishPack directory
        if (pyfishpack_root / "meson.build").exists():
            print(f"DEBUG: Discovered meson module: PyFishPack")
            modules.append(
                {
                    "name": "PyFishPack",
                    "path": pyfishpack_root,
                }
            )

        return modules


class CustomDevelop(develop):
    """Custom develop command that builds meson modules"""

    def run(self):
        # Build meson modules in develop mode
        self.run_command("build_ext")
        super().run()


class CustomInstall(install):
    """Custom install command that ensures meson modules are built"""

    def run(self):
        # Ensure meson modules are built before install
        self.run_command("build_ext")
        super().run()


# Read the README file for the long description
def read_readme():
    """Read README file for long description"""
    try:
        with open("README.md", "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        return "PyFishPack - Python wrapper for FISHPACK Fortran library"


# Configuration for mixed build
setup_config = {
    "name": "PyFishPack",
    "version": "1.0.0",
    "description": "Python wrapper for FISHPACK - Fortran finite difference Poisson solvers",
    "long_description": read_readme(),
    "long_description_content_type": "text/markdown",
    "author": "PyFishPack Development Team",
    "author_email": "your.email@example.com",
    "url": "https://github.com/yourusername/PyFishPack",
    "packages": find_packages(include=['PyFishPack', 'PyFishPack.*']),
    "classifiers": [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Fortran",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    "python_requires": ">=3.8",
    "install_requires": [
        "numpy>=1.20.0",
    ],
    "extras_require": {
        "dev": [
            "pytest",
            "pytest-cov",
            "black",
            "isort",
            "flake8",
        ],
        "build": [
            "meson",
            "ninja",
            "build",
            "wheel",
        ],
    },
    "cmdclass": {
        "build_ext": MesonBuildExt,
        "develop": CustomDevelop,
        "install": CustomInstall,
    },
    # Add a dummy extension for compatibility
    "ext_modules": [
        Extension("PyFishPack._dummy", sources=[
                  "PyFishPack/_dummy.c"], optional=True)
    ],
    "zip_safe": False,  # Due to compiled extensions
}

if __name__ == "__main__":
    setup(**setup_config)
