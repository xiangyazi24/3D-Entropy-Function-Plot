# 3D Entropy Function Plot

Shannon entropy of a three-outcome distribution,

    H(p1, p2, p3) = -Σ p_i log2 p_i,   p1 + p2 + p3 = 1,

plotted as a surface over the probability simplex.

## Interactive web version (`index.html`)

A single-page app built with [three.js](https://threejs.org/) r147. three.js is vendored in `vendor/` (MIT, see `vendor/three-LICENSE`), so the page needs no CDN.
Open `index.html` in any modern browser with WebGL, or serve the repo with GitHub Pages.

- Drag to orbit, scroll or pinch to zoom. Hover or tap the surface to see the distribution and its entropy in bits and nats.
- Height is H in bits. The maximum, log2 3 ≈ 1.585 bits, is at (1/3, 1/3, 1/3).
- Edge curves show the binary entropy h(p), which is the surface restricted to each edge of the simplex.
- Iso-entropy contours have adjustable spacing, with every fifth line drawn heavier. A flat ternary map on the floor repeats them.
- Faint grid lines mark p_i = 0.1, 0.2, and so on.

How it improves on the Mathematica version:

| Mathematica notebook | Web version |
|---|---|
| Rectangular (u, v) grid clipped by `RegionFunction`, which leaves jagged edges | Triangular grid in barycentric coordinates, so the triangle's boundary is exact |
| Surface sampled at 60 plot points, so it looks faceted near the steep edges | Color, contours and normals are computed per pixel from the exact formula and its gradient |
| `"Rainbow"` color map | Perceptually ordered sequential ramp: dark means low entropy, light means high |
| 50 × 50 mesh lines in u and v, which have no meaning for the distribution | Entropy contours and p_i grid lines, both of which you can read |
| Static image | Interactive, with live readout and adjustable height scale |

## Printable STL (`stl/`)

`stl/entropy_surface.stl` is generated directly from the formula by `stl/make_stl.py` (standard-library Python):

- units are millimetres: 100 mm triangle edge, 112 mm tall, 1.2 mm wall; rescale freely in the slicer;
- the top sheet is the exact surface on a triangular barycentric grid, refined toward the edges, with 0·log 0 = 0 on the boundary;
- the bottom sheet is offset inward along the exact normal; the normal is frozen in a thin band at the rim, where the profile s·log(1/s) curves too sharply for an exact offset, and the thickness tapers to zero at the three tips;
- the mesh is closed and manifold (Euler characteristic 2, outward normals); `python3 stl/check_stl.py stl/entropy_surface.stl` verifies this.

Options: `python3 stl/make_stl.py out.stl --edge 150 --wall 1.6 --n 300`.

## Mathematica version (`3d-entropy.nb`)

The original notebook maps the simplex x + y + z = 1 onto the (u, v) plane and uses `ParametricPlot3D`. `entropy_function.stl` is the mesh exported from it, kept for reference. Its two slanted edges come from `RegionFunction` clipping, where Mathematica interpolates boundary heights linearly instead of evaluating H; because H has infinite slope at the edges, those heights come out low (by up to 0.13 bits in the saved plot). The third edge, v = 0, is sampled directly (0·Log2[0] is `Indeterminate`, so Mathematica bisects toward it) and is exact.
