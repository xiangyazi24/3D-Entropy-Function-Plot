"""Generate a printable STL of the ternary entropy surface.

H(p1, p2, p3) = -sum p_i log2 p_i over the probability simplex, laid flat as an
equilateral triangle (side sqrt 2, the simplex's true size) and scaled to
millimetres. The top sheet is the exact surface; the bottom sheet is the same
surface moved inward along its exact unit normal by the wall thickness, so the
wall is equally thick on the vertical rim as on the dome. Near the three
corners the thickness tapers to zero so the legs end in points and the offset
cannot cross itself.

Usage: python3 make_stl.py [output.stl] [--n 220] [--edge 100] [--wall 1.2]
Only the standard library is needed.
"""
import argparse
import math
import struct

LN2 = math.log(2)
R2, R3, R6 = math.sqrt(2), math.sqrt(3), math.sqrt(6)

# Corners in the plane, in simplex units: B = e1 at the origin, C = e2 on the u axis.
V = [(0.0, 0.0), (R2, 0.0), (1 / R2, R6 / 2)]
ALT = R6 / 2  # altitude of the triangle

# Gradient of each barycentric coordinate in the plane (points toward its corner).
def _grads():
    (x1, y1), (x2, y2), (x3, y3) = V
    det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3)
    g1 = ((y2 - y3) / det, (x3 - x2) / det)
    g2 = ((y3 - y1) / det, (x1 - x3) / det)
    return [g1, g2, (-g1[0] - g2[0], -g1[1] - g2[1])]
G = _grads()


def entropy(p):
    return -sum(q * math.log2(q) for q in p if q > 0)


def plane(p):
    return (sum(p[i] * V[i][0] for i in range(3)), sum(p[i] * V[i][1] for i in range(3)))


def normal(p, floor):
    """Unit upward normal of z = H, evaluated at p with every coordinate raised to at least `floor`.

    Across an edge the surface profile is s log(1/s), whose radius of curvature goes to zero
    at the rim, so an exact inward offset would fold over itself there. Freezing the normal in
    a thin band along the edges turns the offset into a local translation, which cannot fold.
    """
    q = [max(x, floor) for x in p]
    s = sum(q)
    q = [x / s for x in q]
    dH = [-(math.log(x) + 1) / LN2 for x in q]
    gx = sum(dH[i] * G[i][0] for i in range(3))
    gy = sum(dH[i] * G[i][1] for i in range(3))
    n = math.sqrt(gx * gx + gy * gy + 1)
    return (-gx / n, -gy / n, 1 / n)


def warp(p, gamma):
    """Pull grid points toward the edges, where the surface is steep. Edges and corners stay put."""
    w = [q ** gamma for q in p]
    s = sum(w)
    return [q / s for q in w]


def build(n, gamma, wall, taper, floor):
    idx = {}
    top, bot = [], []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            p = warp([i / n, j / n, k / n], gamma)
            if (i, j, k).count(0) == 2:  # exact corner
                p = [float(i == n), float(j == n), float(k == n)]
            x, y = plane(p)
            z = entropy(p)
            # distance (in the plane) to the nearest corner decides the taper
            r = min(math.hypot(x - cx, y - cy) for cx, cy in V)
            t = wall * min(1.0, r / taper)
            nx, ny, nz = normal(p, floor) if t > 0 else (0.0, 0.0, 1.0)
            idx[(i, j)] = len(top)
            top.append((x, y, z))
            bot.append((x - t * nx, y - t * ny, z - t * nz) if t > 0 else None)  # None: shares the top vertex

    verts = list(top)
    bot_index = []
    for b, tv in zip(bot, range(len(top))):
        if b is None:
            bot_index.append(tv)
        else:
            bot_index.append(len(verts))
            verts.append(b)

    faces = []
    def tri_top(a, b, c):
        faces.append((a, b, c))
    def tri_bot(a, b, c):
        faces.append((bot_index[a], bot_index[c], bot_index[b]))

    for i in range(n):
        for j in range(n - i):
            a, b, c = idx[(i, j)], idx[(i + 1, j)], idx[(i, j + 1)]
            tri_top(a, b, c); tri_bot(a, b, c)
            if j < n - i - 1:
                a2, b2, c2 = idx[(i + 1, j)], idx[(i + 1, j + 1)], idx[(i, j + 1)]
                tri_top(a2, b2, c2); tri_bot(a2, b2, c2)

    # boundary loop of the top sheet, in the same direction as the top faces' boundary
    loop = [idx[(i, 0)] for i in range(n + 1)]                 # j = 0: from corner (0,0) to (n,0)
    loop += [idx[(n - s, s)] for s in range(1, n + 1)]         # k = 0: to (0,n)
    loop += [idx[(0, n - s)] for s in range(1, n)]             # i = 0: back toward (0,0)
    for a, b in zip(loop, loop[1:] + loop[:1]):
        ba, bb = bot_index[a], bot_index[b]
        # wall between top edge a->b and bottom edge; orient so each edge appears once in each direction
        if bb != b:
            faces.append((b, a, bb))
        if ba != a:
            faces.append((a, ba, bb))
    return verts, faces


def orient_top(verts, faces):
    """Make the top faces wind counter-clockwise seen from above; flip everything if needed."""
    a, b, c = faces[0]
    (x1, y1, _), (x2, y2, _), (x3, y3, _) = verts[a], verts[b], verts[c]
    if (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) < 0:
        return [(f[0], f[2], f[1]) for f in faces]
    return faces


def write_stl(path, verts, faces, scale):
    with open(path, "wb") as f:
        f.write(b"Ternary Shannon entropy surface, H(p1,p2,p3) in bits, units: mm".ljust(80, b" "))
        f.write(struct.pack("<I", len(faces)))
        for a, b, c in faces:
            p = [tuple(scale * x for x in verts[i]) for i in (a, b, c)]
            ux, uy, uz = (p[1][k] - p[0][k] for k in range(3))
            vx, vy, vz = (p[2][k] - p[0][k] for k in range(3))
            nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
            ln = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
            f.write(struct.pack("<3f", nx / ln, ny / ln, nz / ln))
            for q in p:
                f.write(struct.pack("<3f", *q))
            f.write(b"\0\0")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", nargs="?", default="entropy_surface.stl")
    ap.add_argument("--n", type=int, default=220, help="grid subdivisions per edge")
    ap.add_argument("--gamma", type=float, default=1.6, help="edge refinement strength (1 = uniform)")
    ap.add_argument("--edge", type=float, default=100.0, help="edge length of the base triangle, mm")
    ap.add_argument("--wall", type=float, default=1.2, help="wall thickness, mm")
    ap.add_argument("--floor", type=float, default=0.005, help="normal is frozen where a coordinate is below this")
    args = ap.parse_args()

    scale = args.edge / R2                      # mm per simplex unit
    wall = args.wall / scale                    # thickness in simplex units
    verts, faces = build(args.n, args.gamma, wall, 3 * wall, args.floor)
    faces = orient_top(verts, faces)
    write_stl(args.out, verts, faces, scale)
    xs, ys, zs = zip(*verts)
    print("%s: %d triangles, %.1f x %.1f x %.1f mm, wall %.2f mm"
          % (args.out, len(faces), scale * (max(xs) - min(xs)), scale * (max(ys) - min(ys)),
             scale * (max(zs) - min(zs)), args.wall))


if __name__ == "__main__":
    main()
