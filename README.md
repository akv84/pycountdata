# pycountdata

Python wrapper for the R [countdata](https://CRAN.R-project.org/package=countdata) package —
beta-binomial (`bb_test`) and inverted beta-binomial (`ibb_test`) tests for
differential count data analysis.

Same C code, same precision, with OpenMP multi-threading.

## Installation

```bash
# From GitHub
pip install git+https://github.com/akv84/pycountdata.git

# From local clone
git clone https://github.com/akv84/pycountdata.git
cd pycountdata
pip install .
# or
python setup.py install
```

**Requirements:** C compiler (`gcc`) and `numpy`. OpenMP is auto-detected.

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
| `bb_test(x, tx, group, ...)` | `{'p.value'}` | Unpaired beta-binomial test, ≥2 groups |
| `ibb_test(x, tx, group, ...)` | `{'p.value', 'fc'}` | Paired inverted beta-binomial test |
| `normalize(d)` | ndarray | Column-total normalization |
| `fold_change(d1, d2)` | ndarray | Per-row fold change |

**Common parameters:** `alternative` (`'two.sided'`/`'less'`/`'greater'`),
`n_threads` (int, -1=all cores), `verbose` (bool).

## Exact R precision (optional)

By default `bb_test` matches R exactly; `ibb_test` matches within ~0.05%.
For bit-identical `ibb_test`, build standalone `libRmath.so` before installing:

```bash
bash scripts/build_rmath.sh   # compiles libRmath.so into $CONDA_PREFIX/lib
pip install --force-reinstall .
```

## How it works

During `pip install`, `setup.py`:

1. Uses the pre-generated `q15.h` in the repo (32768-point Gauss-Legendre
   quadrature for `ibb_test`). Regenerates from numpy if missing.
2. Auto-detects OpenMP, `-march=native`, and `libRmath.so`.
3. Compiles `bb.c` + `ibb.c` → `_countdata.so` inside the installed package.

## Citations

- **bb_test:** Pham TV, Piersma SR, Warmoes M, Jimenez CR (2010). *Bioinformatics*, 26(3):363-369.
- **ibb_test:** Pham TV, Jimenez CR (2012). *Bioinformatics*, 28(18):i596-i602.

## License

GPL-3.0. Original C code by Thang V. Pham.