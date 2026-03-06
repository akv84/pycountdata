#!/bin/bash
# build_rmath.sh - Build standalone libRmath.so from R source.
set -e

if [ -z "$CONDA_PREFIX" ]; then
    echo "ERROR: Activate conda env first: conda activate countdata"
    exit 1
fi

rm -f "$CONDA_PREFIX/lib/libRmath.so"

R_VERSION=$(R --vanilla --slave -e "cat(paste0(R.version\$major,'.',R.version\$minor))" 2>/dev/null)
R_MAJOR=$(echo $R_VERSION | cut -d. -f1)
R_INCLUDE=$(R --vanilla --slave -e "cat(R.home('include'))" 2>/dev/null)
CC="${CC:-gcc}"

echo "=== Building standalone libRmath.so ==="
echo "  R version: $R_VERSION"
echo "  R include: $R_INCLUDE"
echo "  Compiler:  $CC"

WORKDIR=$(mktemp -d)
trap "rm -rf $WORKDIR" EXIT
cd "$WORKDIR"

echo ""
echo "Downloading R-${R_VERSION} source..."
URL="https://cran.r-project.org/src/base/R-${R_MAJOR}/R-${R_VERSION}.tar.gz"
if command -v curl &> /dev/null; then curl -fSL "$URL" -o R.tar.gz
elif command -v wget &> /dev/null; then wget -q "$URL" -O R.tar.gz
else echo "ERROR: Need curl or wget"; exit 1; fi

echo "Extracting nmath source..."
tar xzf R.tar.gz "R-${R_VERSION}/src/nmath/" "R-${R_VERSION}/src/include/" 2>/dev/null || \
tar xzf R.tar.gz "R-${R_VERSION}/src/nmath/" "R-${R_VERSION}/src/include/"

NMATH="R-${R_VERSION}/src/nmath"

# Detect available math functions
echo ""
echo "Detecting system math functions..."
cat > config.h << 'CFGEOF'
#define MATHLIB_STANDALONE 1
CFGEOF

for func in log1p expm1 hypot cospi sinpi tanpi; do
    echo "#include <math.h>" > "test_${func}.c"
    echo "double test(double x) { return ${func}(x); }" >> "test_${func}.c"
    UPPER=$(echo $func | tr a-z A-Z)
    if $CC -shared -fPIC -lm "test_${func}.c" -o /dev/null 2>/dev/null; then
        echo "#define HAVE_${UPPER} 1" >> config.h
        [ "$func" = "log1p" ] && echo "#define HAVE_WORKING_LOG1P 1" >> config.h
        echo "  $func: yes"
    else
        echo "  $func: no (R fallback)"
    fi
done

# Minimal stubs — only for symbols nmath references but doesn't define itself.
# sexp.c defines exp_rand, snorm.c defines norm_rand and N01_kind.
# We only need: unif_rand, R_unif_index, GetRNGstate, PutRNGstate.
echo ""
echo "Creating minimal RNG stubs..."
cat > rng_stubs.c << 'STUBEOF'
/* Stubs for R RNG symbols not defined in nmath.
 * Only unif_rand and R_unif_index are truly missing.
 * exp_rand, norm_rand, N01_kind are in sexp.c/snorm.c. */
double unif_rand(void) { return 0.5; }
double R_unif_index(double dn) { (void)dn; return 0.0; }
void GetRNGstate(void) {}
void PutRNGstate(void) {}
STUBEOF
$CC -O3 -fPIC -DHAVE_CONFIG_H -DMATHLIB_STANDALONE -I. -I"$R_INCLUDE" -c rng_stubs.c -o rng_stubs.o
echo "  OK"

# Compile nmath
echo ""
echo "Compiling nmath..."
CFLAGS="-O3 -fPIC -DHAVE_CONFIG_H -DMATHLIB_STANDALONE"
OBJS="rng_stubs.o"

for src in $NMATH/*.c; do
    base=$(basename "$src" .c)
    if $CC $CFLAGS -I. -I"$R_INCLUDE" -c "$src" -o "${base}.o" 2>/dev/null; then
        OBJS="$OBJS ${base}.o"
    fi
done
echo "  Compiled: $(echo $OBJS | wc -w) files"

# Link
echo ""
echo "Linking..."
$CC -shared -o libRmath.so $OBJS -lm
echo "  Built: libRmath.so ($(du -h libRmath.so | cut -f1))"

# Verify symbols
echo ""
echo "Checking key symbols..."
for sym in lgammafn digamma trigamma pnorm5 pchisq; do
    if nm -D libRmath.so 2>/dev/null | grep -q " T.*${sym}"; then
        printf "  %-12s OK\n" "$sym"
    else
        printf "  %-12s MISSING\n" "$sym"
    fi
done

# Load and functional test
echo ""
echo "Testing..."
if python3 -c "
import ctypes
lib = ctypes.CDLL('$PWD/libRmath.so')
f = lib.lgammafn
f.restype = ctypes.c_double
f.argtypes = [ctypes.c_double]
assert abs(f(1.0)) < 1e-15, 'lgammafn(1) failed'
d = lib.digamma
d.restype = ctypes.c_double
d.argtypes = [ctypes.c_double]
assert abs(d(1.0) - (-0.5772156649015329)) < 1e-12, 'digamma(1) failed'
print('  lgammafn, digamma: verified')
" 2>&1; then
    cp libRmath.so "$CONDA_PREFIX/lib/"
    echo ""
    echo "=== Installed: $CONDA_PREFIX/lib/libRmath.so ==="
    echo ""
    echo "Now run:"
    echo "  bash build.sh --rmath"
    echo "  python3 test_countdata.py"
else
    echo "  FAILED"
    exit 1
fi