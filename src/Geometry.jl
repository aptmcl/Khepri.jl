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
  let v0 = p1 - p0,
      v1 = rotated_v(v0, pi/2),
      c = loc_from_o_vx_vy(p0, v0, v1)
    rectangle(c-vy(w/2, c.cs), distance(p0, p1), w)
  end
offset(ps::Locs, d::Real) = error("BUM")
offset(ps::Locs, d::Real, closed) =
  let qs = closed ? [ps[end], ps..., ps[1]] : ps,
      vs = map((p0, p1) -> rotated_v(unitized(p1 - p0)*d, pi/2), qs[2:end], qs[1:end-1]),
      ws = map(v_in_v, vs[1:end-1], vs[2:end])
    map(+, ps, closed ? ws : [vs[1], ws..., vs[end]])
  end

offset(path::OpenPolygonalPath, d::Real) = open_polygonal_path(offset(path.vertices, d, false))
offset(path::ClosedPolygonalPath, d::Real) = closed_polygonal_path(offset(path.vertices, d, true))
offset(l::Line, d::Real) = line(offset(l.vertices, d, false))

export offset

# Polygon combination

closest_vertices_indexes(pts1, pts2) =
  # This is a brute force method. There are better algorithms to do this.
  let min_dist = Inf,
      min_i = nothing,
      min_j = nothing
    for (i, pt1) in enumerate(pts1)
      for (j, pt2) in enumerate(pts2)
        let dist = distance(pt1, pt2)
          if dist < min_dist
            min_dist = dist
            min_i = i
            min_j = j
          end
        end
      end
    end
    min_i, min_j
  end

#=
using Test
pts1 = regular_polygon_vertices(4, xy(1,2), 4)
pts2 = regular_polygon_vertices(4, xy(2,2), 1)
pts3 = regular_polygon_vertices(5, xy(1,4), 1)
pts4 = regular_polygon_vertices(6, xy(0,1), 1)

polygon(pts1)
polygon(pts2)
polygon(pts3)
polygon(pts4)

@test closest_vertices_indexes(pts1, pts2) == (1,1)
@test closest_vertices_indexes(pts1, pts3) == (2,2)
@test closest_vertices_indexes(pts1, pts4) == (4,6)
@test closest_vertices_indexes(pts2, pts3) == (2,5)
=#

subtract_polygon_vertices(pts1, pts2) =
  let (i, j) = closest_vertices_indexes(pts1, pts2)
    [pts1[1:i]..., reverse([pts2[j:end]..., pts2[1:j]...])..., pts1[i:end]...]
  end

export subtract_polygon_vertices

#=
pts1 = subtract_polygon_vertices(pts1, pts2)
pts1 = subtract_polygon_vertices(pts1, pts3)
pts1 = subtract_polygon_vertices(pts1, pts4)
polygon(pts1)
for (i,p) in enumerate(pts1)
  text(string(i), p, 0.1)
  sleep(1)
end
=#
