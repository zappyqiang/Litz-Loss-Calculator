# Litz Loss Calculator

Browser calculator for the AC resistance and winding loss of litz wire at a single frequency.

- **`index.html`** – the calculator (single file, no build step; open it directly or serve it with GitHub Pages).
- **`core.js`** – the same physics as plain Node-importable JavaScript (`litz()`, `strandAC()`).
- **`fem_check.py`** – 2D finite-element check in NGSolve (hex-packed bundle, every strand carrying I/N). `pip install ngsolve scipy numpy`, then `python fem_check.py 2` (19 strands) or `3` (37 strands).
- **`check_scipy.py`** – checks the Bessel-function code in `core.js` against SciPy.

## Model

Each strand (radius r) carries I/N. With γ = (1+j)/δ and δ = √(2ρ/ωμ₀):

```
Rdc    = ρ·l / (N·π·r²)
F_skin = Re[ γr·I₀(γr) / (2·I₁(γr)) ]
H_i²   = I² / (8π²·R_b²)          (mean-square internal field, rms)
P_prox = N · l · k_p(f) · H_i²    (k_p: exact loss of a round strand in a uniform transverse field)
Rac    = (F_skin·Rdc·I² + P_prox) / I²
```

Not included: fields from outside the bundle (neighbouring turns, gaps, cores), uneven current sharing, insulation losses.

## Validation against 2D FEM (37 strands of 0.1 mm, copper fill 69 %)

| f | strand d/δ | FEM Rac/Rdc | this model | diff |
|---|---|---|---|---|
| 100 kHz | 0.48 | 1.0101 | 1.0105 | +0.04 % |
| 300 kHz | 0.83 | 1.0900 | 1.0939 | +0.35 % |
| 1 MHz | 1.51 | 1.933 | 2.009 | +3.9 % |
| 2 MHz | 2.14 | 4.066 | 4.644 | +14 % |
| 4 MHz | 3.03 | 8.257 | 11.585 | +40 % |

Accurate to about 1 % up to d/δ ≈ 1; it overestimates loss (conservative) above d/δ ≈ 2.
