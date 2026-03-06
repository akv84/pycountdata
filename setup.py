"""
setup.py for pycountdata.

pip install name:  pycountdata
import name:       countdata

During install this script:
  1. Generates q15.h if missing (32768-pt Gauss-Legendre quadrature for ibb_test)
  2. Auto-detects OpenMP support
  3. Auto-detects libRmath.so in $CONDA_PREFIX/lib for exact R precision
  4. Compiles bb.c + ibb.c -> countdata/_countdata.so
"""

import os
import subprocess
import platform
import tempfile
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

HERE = os.path.dirname(os.path.abspath(__file__))
CSRC = os.path.join(HERE, "countdata", "csrc")
PKG = os.path.join(HERE, "countdata")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _try_compile(cc, flags, code="int main(){return 0;}\n"):
    """Test if compiler accepts given flags."""
    fd, tmp = tempfile.mkstemp(suffix=".c")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(code)
        return subprocess.call(
            [cc] + flags + [tmp, "-o", os.devnull],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        ) == 0
    except Exception:
        return False
    finally:
        os.unlink(tmp)


def _find_rmath():
    """Find standalone libRmath.so and Rmath.h. Returns (inc_dir, lib_dir) or (None, None)."""
    conda = os.environ.get("CONDA_PREFIX", "")
    inc_dirs = [
        os.path.join(conda, "lib", "R", "include"),
        os.path.join(conda, "include"),
        "/usr/share/R/include",
        "/usr/local/lib/R/include",
    ]
    lib_dirs = [
        os.path.join(conda, "lib"),
        os.path.join(conda, "lib64"),
        "/usr/lib",
        "/usr/lib/x86_64-linux-gnu",
        "/usr/local/lib",
    ]
    inc = next((d for d in inc_dirs if os.path.isfile(os.path.join(d, "Rmath.h"))), None)
    lib = next((d for d in lib_dirs if os.path.isfile(os.path.join(d, "libRmath.so"))), None)
    return (inc, lib) if inc and lib else (None, None)


def _generate_q15(path, n=32768):
    """Generate Gauss-Legendre quadrature nodes and log-weights."""
    print(f"    generating q15.h ({n} points, needs ~2 GB RAM) ...")
    import numpy as np
    try:
        from scipy.special import roots_legendre
        z, w = roots_legendre(n)
    except ImportError:
        z, w = np.polynomial.legendre.leggauss(n)
    lw = np.log(w)
    with open(path, "w") as f:
        f.write(f"const int _Z = {n};\n")
        f.write("const TYPE _qz[] = {")
        for i, v in enumerate(z):
            f.write(f",\n{v:.30e}" if i else f"\n{v:.30e}")
        f.write("};\nconst TYPE _qw[] = {")
        for i, v in enumerate(lw):
            f.write(f",\n{v:.30e}" if i else f"\n{v:.30e}")
        f.write("};\n")
    print(f"    done ({os.path.getsize(path) // 1024} KB)")


# ---------------------------------------------------------------------------
# Compile
# ---------------------------------------------------------------------------

def compile_countdata():
    """Compile bb.c + ibb.c into _countdata.so."""

    print("=" * 60)
    print("  pycountdata: compiling C library")
    print("=" * 60)

    # --- q15.h -----------------------------------------------------------
    q15 = os.path.join(CSRC, "q15.h")
    if not os.path.exists(q15) or os.path.getsize(q15) < 50_000:
        try:
            _generate_q15(q15)
        except Exception as e:
            print(f"  WARNING: q15.h generation failed ({e})")
            print("  ibb_test will not work. See README for manual steps.")

    # --- compiler --------------------------------------------------------
    cc = os.environ.get("CC", "gcc")
    srcs = [os.path.join(CSRC, f) for f in ("bb.c", "ibb.c")]
    ext = ".dylib" if platform.system() == "Darwin" else ".so"
    out = os.path.join(PKG, f"_countdata{ext}")

    cflags = ["-O3", "-fPIC", "-Wall", f"-I{CSRC}"]
    ldflags = ["-shared", "-lm"]

    # --- OpenMP ----------------------------------------------------------
    omp_test = "#include <omp.h>\nint main(){return omp_get_num_threads();}\n"
    if _try_compile(cc, ["-fopenmp", "-lm"], omp_test):
        cflags.append("-fopenmp")
        ldflags.append("-fopenmp")
        print("    OpenMP:  enabled")
    else:
        print("    OpenMP:  not available (single-thread fallback)")

    # --- march=native ----------------------------------------------------
    if _try_compile(cc, ["-march=native"]):
        cflags.append("-march=native")

    # --- libRmath --------------------------------------------------------
    rinc, rlib = _find_rmath()
    if rinc and rlib:
        cflags += ["-DUSE_RMATH", f"-I{rinc}"]
        ldflags += [f"-L{rlib}", "-lRmath", f"-Wl,-rpath,{rlib}"]
        print(f"    Rmath:   {rlib}/libRmath.so (exact R precision)")
    else:
        print("    Rmath:   not found (standalone C99 math)")

    # --- compile ---------------------------------------------------------
    cmd = [cc] + cflags + srcs + ["-o", out] + ldflags
    print(f"    compile: {cc} ... -o _countdata{ext}")
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        # retry without OpenMP / march=native
        print("    retrying with minimal flags ...")
        cmd = [cc, "-O2", "-fPIC", f"-I{CSRC}"] + srcs + ["-o", out, "-shared", "-lm"]
        subprocess.check_call(cmd)

    print(f"    output:  {out} ({os.path.getsize(out) // 1024} KB)")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Custom build command
# ---------------------------------------------------------------------------

class BuildWithC(build_py):
    """Compile C shared library before building Python package."""
    def run(self):
        compile_countdata()
        super().run()


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

setup(
    cmdclass={"build_py": BuildWithC},
    packages=find_packages(include=["countdata*"]),
    package_data={"countdata": [
        "_countdata.so", "_countdata.dylib",
        "csrc/*.c", "csrc/*.h",
    ]},
)