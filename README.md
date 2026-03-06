# pycountdata

Python wrapper for the R [countdata](https://CRAN.R-project.org/package=countdata) package —
beta-binomial (`bb_test`) and inverted beta-binomial (`ibb_test`) tests for
differential count data analysis.

Same C code, same precision, with OpenMP multi-threading.

## Installation

### Quick install from GitHub

```bash
pip install git+https://github.com/akv84/pycountdata.git
```

### Install from cloned repository

```bash
git clone https://github.com/akv84/pycountdata.git
cd pycountdata
pip install .
```

### Install from downloaded source

```bash
cd pycountdata
pip install .

# Alternative
python setup.py install
```

### Requirements

- **Python** >= 3.8
- **numpy** >= 1.20 (installed automatically)
- **C compiler** (`gcc`) — must be available in PATH

### What happens during install

`setup.py` handles everything automatically:

1. Detects `gcc` and checks for OpenMP support (enables multi-threading if available)
2. Uses pre-generated `q15.h` from the repository (32768-point Gauss-Legendre
   quadrature needed by `ibb_test`). Regenerates from numpy if file is missing.
3. Checks for `libRmath.so` in `$CONDA_PREFIX/lib` — links against it if found
   for exact R precision (see below)
4. Compiles `bb.c` + `ibb.c` into `countdata/_countdata.so`

No manual compilation steps are needed.

### Exact R precision (optional)

By default, `bb_test` matches R exactly and `ibb_test` matches within ~0.05%.
This is sufficient for most use cases.

For bit-identical results with R's `ibb.test`, build R's standalone math library
before installing:

```bash
# Requires: conda environment with r-base
conda create -n countdata -c conda-forge python numpy r-base compilers -y
conda activate countdata

# Build libRmath.so (downloads R source from CRAN, compiles math library)
cd pycountdata
bash scripts/build_rmath.sh

# Install (setup.py auto-detects libRmath.so)
pip install .
```

### Verify installation

```python
python -c "from countdata import bb_test; print('OK')"
```

## Usage

```python
import numpy as np
from countdata import bb_test, ibb_test, normalize, fold_change

# ── Unpaired beta-binomial test (2+ groups) ──
x = np.array([[1, 5, 1, 7, 9, 11, 2, 10],
               [11, 1, 9, 1, 1, 1, 1, 1]], dtype=float)
tx = np.full_like(x, 100.0)
group = [1, 1, 1, 1, 2, 2, 2, 2]

result = bb_test(x, tx, group)
print(result['p.value'])  # [0.0837, 0.1248]

# ── Paired inverted beta-binomial test (2 equal-size groups) ──
result = ibb_test(x[:, :4], tx[:, :4], [1, 1, 2, 2])
print(result['p.value'], result['fc'])

# ── One-sided tests ──
bb_test(x, tx, group, alternative='less')
bb_test(x, tx, group, alternative='greater')

# ── Control threading ──
bb_test(x, tx, group, n_threads=8, verbose=False)
```

## API

| Function | Returns | Description |
|----------|---------|-------------|
| `bb_test(x, tx, group, ...)` | `{'p.value'}` | Unpaired beta-binomial test, >= 2 groups |
| `ibb_test(x, tx, group, ...)` | `{'p.value', 'fc'}` | Paired inverted beta-binomial test |
| `normalize(d)` | ndarray | Column-total normalization |
| `fold_change(d1, d2)` | ndarray | Per-row fold change |

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `x` | array (K x N) | — | Count matrix (K genes, N samples) |
| `tx` | array (K x N) | — | Total counts per sample, must be positive |
| `group` | array (N,) | — | Group labels |
| `alternative` | str | `'two.sided'` | `'two.sided'`, `'less'`, or `'greater'` |
| `n_threads` | int | `-1` | Number of threads (`-1` = all cores) |
| `verbose` | bool | `True` | Print progress |
| `BIG` | float | `1e4` | Fold-change cap for zero counts (`ibb_test` only) |

## Troubleshooting

**`gcc: command not found`**
Install a C compiler. On Ubuntu/Debian: `sudo apt install build-essential`.
In conda: `conda install -c conda-forge compilers`.

**`_countdata.so not found` after install**
The C compilation may have failed silently. Reinstall with verbose output:
`pip install -v .` and check for compiler errors.

**`libgomp.so: cannot open shared object file`**
OpenMP runtime missing. Fix: `conda install -c conda-forge libgomp`.
Or reinstall — `setup.py` auto-disables OpenMP and falls back to single-thread.

**`undefined symbol: lgammafn` at runtime**
You installed with `libRmath` but it is not in the library path.
Fix: `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`
Or reinstall without libRmath: `pip install --force-reinstall .`

## Citations

- **bb_test:** Pham TV, Piersma SR, Warmoes M, Jimenez CR (2010). On the beta-binomial
  model for analysis of spectral count data. *Bioinformatics*, 26(3):363-369.
- **ibb_test:** Pham TV, Jimenez CR (2012). An accurate paired sample test for count
  data. *Bioinformatics*, 28(18):i596-i602.

## License

GPL-3.0. Original C code by Thang V. Pham.