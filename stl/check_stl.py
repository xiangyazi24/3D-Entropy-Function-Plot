"""Sanity checks for a binary STL: watertight, consistently oriented, outward, no folded bottom sheet."""
import collections, math, struct, sys

d = open(sys.argv[1], "rb").read()
n = struct.unpack("<I", d[80:84])[0]
assert len(d) == 84 + 50 * n, "size mismatch"
key = {}
verts, faces = [], []
for i in range(n):
    o = 84 + 50 * i + 12
    f = []
    for k in range(3):
        v = struct.unpack("<3f", d[o + 12 * k:o + 12 * k + 12])
        if v not in key:
            key[v] = len(verts); verts.append(v)
        f.append(key[v])
    faces.append(f)
directed = collections.Counter()
for a, b, c in faces:
    for e in ((a, b), (b, c), (c, a)):
        directed[e] += 1
undirected = collections.Counter()
for (a, b), c in directed.items():
    undirected[tuple(sorted((a, b)))] += c
bad_multi = sum(1 for c in directed.values() if c != 1)
bad_pair = sum(1 for (a, b) in directed if (b, a) not in directed)
nonmanifold = sum(1 for c in undirected.values() if c != 2)
degenerate = sum(1 for f in faces if len(set(f)) < 3)
V, E, F = len(verts), len(undirected), len(faces)
vol = 0.0
for a, b, c in faces:
    (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = verts[a], verts[b], verts[c]
    vol += (x1 * (y2 * z3 - z2 * y3) - y1 * (x2 * z3 - z2 * x3) + z1 * (x2 * y3 - y2 * x3)) / 6
xs, ys, zs = zip(*verts)
print("triangles %d, vertices %d" % (F, V))
print("bbox %.2f x %.2f x %.2f" % (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)))
print("edges without an opposite partner: %d, repeated directed edges: %d, non-manifold edges: %d, degenerate faces: %d"
      % (bad_pair, bad_multi, nonmanifold, degenerate))
print("Euler characteristic V - E + F = %d (sphere = 2)" % (V - E + F))
print("signed volume %.1f (positive = outward normals)" % vol)
