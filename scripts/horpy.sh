#!/bin/sh
set -e  # Exit immediately on error
set -u  # Treat unset variables as errors

# -----------------------------------------------------------------------------
# Resolve key directories
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# -----------------------------------------------------------------------------
# Ask user where to place the virtual environment
# -----------------------------------------------------------------------------
DEFAULT_PYENV_DIR="$PROJECT_ROOT/pyenv"

if [ -n "${CI:-}" ]; then
    # Non-interactive mode for CI
    echo "Running in CI mode — using default virtual environment path."
    PYENV_DIR="$DEFAULT_PYENV_DIR"
else
    echo "Default virtual environment directory: $DEFAULT_PYENV_DIR"
    read -r -p "Enter custom path for virtual environment (or press Enter to use default): " USER_PYENV_DIR
    if [ -n "${USER_PYENV_DIR:-}" ]; then
        PYENV_DIR="$(realpath "$USER_PYENV_DIR")"
    else
        PYENV_DIR="$DEFAULT_PYENV_DIR"
    fi
fi

PYTHON="$PYENV_DIR/bin/python3"
F2PY_MODULE="numpy.f2py"

# Make sure venv executables (like meson, ninja, pip) are on PATH
export PATH="$PYENV_DIR/bin:$PATH"

echo "----------------------------------------"
echo " HORMONE f2py build script"
echo " Project root: $PROJECT_ROOT"
echo " Script dir:   $SCRIPT_DIR"
echo "----------------------------------------"

# -----------------------------------------------------------------------------
# Ensure Python command is available
# -----------------------------------------------------------------------------

if command -v python3 >/dev/null 2>&1; then
    PYTHON_ORG=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON_ORG=python
else
    echo "Error: Python not found"
    exit 1
fi

# -----------------------------------------------------------------------------
# Ensure Python virtual environment exists
# -----------------------------------------------------------------------------

if [ ! -d "$PYENV_DIR" ]; then
    echo "Creating virtual environment at $PYENV_DIR ..."
    $PYTHON_ORG -m venv "$PYENV_DIR"
fi

# Activate venv by using its Python directly
if [ ! -x "$PYTHON" ]; then
    echo "Error: Python not found in venv ($PYTHON)"
    exit 1
fi

# -----------------------------------------------------------------------------
# Ensure necessary python libraries are installed
# -----------------------------------------------------------------------------
echo "Checking Python environment..."
"$PYTHON" -m pip install --upgrade pip >/dev/null
"$PYTHON" -m pip install -q numpy meson ninja

# -----------------------------------------------------------------------------
# Define paths
# -----------------------------------------------------------------------------
PYF_FILE="$SCRIPT_DIR/horpy.pyf"
SO_TARGET="$SCRIPT_DIR/horpy.so"
F2PY="$PYTHON -m $F2PY_MODULE"

# -----------------------------------------------------------------------------
# Clean up old files
# -----------------------------------------------------------------------------
[ -f "$PYF_FILE" ] && rm "$PYF_FILE" && echo "Removed old $PYF_FILE"

# -----------------------------------------------------------------------------
# Generate .pyf signature interface
# -----------------------------------------------------------------------------
echo "Generating .pyf interface file..."
$F2PY -h "$PYF_FILE" -m horpy "$PROJECT_ROOT/src/main/modules.f90" "$SCRIPT_DIR/horpy.f90"

# -----------------------------------------------------------------------------
# Compile Fortran code into shared object
# -----------------------------------------------------------------------------
echo "Compiling Fortran sources..."
cd "$PROJECT_ROOT/src/main"

$F2PY -c "$SCRIPT_DIR/horpy.pyf" -m horpy \
    --fcompiler=gfortran --f90flags='-fopenmp' -lgomp \
    mpi_utils.F90 derived_types.f90 modules.f90 profiler.f90 \
    mpi_domain.F90 io.F90 utils.f90 conserve.f90 ionization.f90 \
    eos.f90 matrix_vars.F90 opacity.f90 radiation_utils.f90 matrix_coeffs.f90 \
    matrix_utils.f90 miccg.f90 fluxlimiter.f90 hlldflux.f90 dirichlet.f90 \
    particles.f90 petsc_solver.F90 matrix_solver.F90 radiation.f90 cooling.f90 \
    sinks.f90 externalforce.f90 composition.f90 star.f90 shockfind.f90 \
    output.f90 smear.f90 input.f90 setup.f90 readbin.f90 timestep.f90 \
    tests.f90 checksetup.f90 allocations.f90 \
    "$SCRIPT_DIR/horpy.f90"

# -----------------------------------------------------------------------------
# Move resulting .so file into the venv's site-packages/
# -----------------------------------------------------------------------------
SITE_PACKAGES_DIR="$($PYTHON -c 'import site; print(site.getsitepackages()[0])')"

if [ ! -d "$SITE_PACKAGES_DIR" ]; then
    echo "Error: site-packages directory not found in venv."
    exit 1
fi

echo "Installing horpy.so into $SITE_PACKAGES_DIR ..."
mv horpy*.so "$SITE_PACKAGES_DIR/horpy.so"
rm $SCRIPT_DIR/horpy.pyf

echo "✅ Installed horpy.so to virtual environment site-packages"

echo "✅ Build complete!"
echo "   Virtual env:   $PYENV_DIR"
echo "----------------------------------------"
