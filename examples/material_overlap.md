# Conservative common-support integration for diagnostic comparisons

Intersects the actual initial convex Q1 material polygons of two meshes with
Sutherland-Hodgman clipping and bounding-box candidate rejection. The initial
collapsed-center quadrilaterals are triangles with repeated vertices. Each
positive-area overlap is triangulated around its interior mean point and
integrated with a positive tensor Gauss/Duffy rule, emitted in bounded chunks.
No fluid or field is remapped, evolved or altered by this analysis helper.

Cell-edge discontinuities in pressure, density, heat and magnetic gradients now
lie on integration boundaries rather than being missed by unrelated Gauss nodes.
Only zero/roundoff-area intersections below 1e-26 times the squared coordinate
scale are omitted. The comparison independently checks the integrated physical
volume against the complete smaller polygon. Integration-order controls remain
required for smooth but non-polynomial deformed-coordinate Jacobian factors.
