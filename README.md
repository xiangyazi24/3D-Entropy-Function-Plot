# 3D Entropy Function Plot

Shannon entropy of a three-outcome distribution,

    H(p1, p2, p3) = -Σ p_i log2 p_i,   p1 + p2 + p3 = 1,

plotted as a surface over the probability simplex.

## Interactive web version (`index.html`)

A self-contained page built with [three.js](https://threejs.org/) (loaded from the jsDelivr CDN).
Open `index.html` in any modern browser, or enable GitHub Pages for this repo to host it.

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

## Mathematica version (`3d-entropy.nb`)

The original notebook maps the simplex x + y + z = 1 onto the (u, v) plane and uses `ParametricPlot3D`. `entropy_function.stl` is a mesh exported from it.
