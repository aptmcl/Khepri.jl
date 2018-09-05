# Geometric Utilities

export show_cs

show_cs(p, scale=1) =
    let rcyl = scale/10
        rcon = scale/5
        lcyl = scale
        lcon = scale/5
        px = add_x(p, 3*lcyl)
        py = add_y(p, 2*lcyl)
        pz = add_z(p, 1*lcyl)
        union(cylinder(p, rcyl, px),
              cone(px, rcon, add_x(px, lcon)),
              cylinder(p, rcyl, py),
              cone(py, rcon, add_y(py, lcon)),
              cylinder(p, rcyl, pz))
    end

project_to_world(surf) =
    transform(surf, inverse_transformation(frame_at(surf, 0, 0)))

#project_to_world(surface_polygon(xyz(1,1,1), xyz(10,1,1), xyz(10,1,5), xyz(1,1,5)))

#=

Given a poligonal line described by its vertices, we need to compute another
polygonal line that is parallel to the first one.

=#

v_in_v(v0, v1) =
    let v = v0 + v1
        v*dot(v0, v0)/dot(v, v0)
    end

rotated_v(v, alpha) =
    vpol(pol_rho(v), pol_phi(v) + alpha)

centered_rectangle(p0, w, p1) =
    let v0 = p1 - p0
        v1 = rotated_v(v0, pi/2)
        c = loc_from_o_vx_vy(p0, v0, v1)
        rectangle(c-vy(w/2, c.cs), distance(p0, p1), w)
    end

offset(ps::Locs, d::Real) =
    let vs = map((p0, p1) -> rotated_v(unitized(p1 - p0)*d, pi/2), ps[2:end], ps[1:end-1])
        vs = [vs[1], map(v_in_v, vs[1:end-1], vs[2:end])..., vs[end]]
        map(+, ps, vs)
    end

offset(l::Line, d::Real) = line(offset(l.vertices, d))

export offset
