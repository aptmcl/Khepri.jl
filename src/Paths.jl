# Paths are an important concept for BIM (and other things)
export open_path,
       closed_path,
       open_path_ops,
       closed_path_ops,
       MoveOp,
       MoveToOp,
       LineOp,
       LineToOp,
       ArcOp,
       CloseOp,
       arc_path,
       circular_path,
       rectangular_path,
       centered_rectangular_path,
       open_polygonal_path,
       closed_polygonal_path,
       polygonal_path,
       open_spline_path,
       closed_spline_path,
       spline_path,
       path_set,
       open_path_sequence,
       closed_path_sequence,
       path_sequence,
       translate,
       stroke,
       fill,
       path_length,
       location_at_length,
       path_start,
       path_end,
       subpath,
       subpaths,
       join_paths,
       subtract_paths,
       path_vertices

path_tolerance = Parameter(1e-10)
coincident_path_location(p1::Loc, p2::Loc) = distance(p1, p2) < path_tolerance()

abstract type Path end

# A path has a domain

getindex(p::Path, i::Real) = location_at_length(p, i)
firstindex(p::Path) = 0
lastindex(p::Path) = path_length(p)
getindex(p::Path, i::ClosedInterval) = subpath(p, i.left, i.right)
path_start(p::Path) = location_at_length(p, 0)
path_end(p::Path) = location_at_length(p, path_length(p))



abstract type OpenPath <: Path end
abstract type ClosedPath <: Path end

is_closed_path(path::Path) = false
is_closed_path(path::ClosedPath) = true
# To avoid duplicate convertions, we also deal with Locs
is_closed_path(pts::Locs) = coincident_path_location(pts[1], pts[end])

struct ArcPath <: OpenPath
    center::Loc
    radius::Real
    start_angle::Real
    amplitude::Real
end
arc_path(center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi*1.0) =
  # pi*1.0 just to circunvent a Julia bug https://github.com/JuliaLang/julia/issues/31949.
  false ? #amplitude < 0 ?
    ArcPath(center, radius, start_angle + amplitude, - amplitude) :
    ArcPath(center, radius, start_angle, amplitude)
struct CircularPath <: ClosedPath
    center::Loc
    radius::Real
end
circular_path(Center::Loc=u0(), Radius::Real=1; center::Loc=Center, radius::Real=Radius) = CircularPath(center, radius)
map_division(f::Function, path::CircularPath, n::Integer) =
  map_division(ϕ->f(loc_from_o_vx_vy(add_pol(path.center, path.radius, ϕ),
                                     vpol(1, ϕ+π, path.center.cs),
                                     vz(1, path.center.cs))),
               0, 2π, n, false)

struct RectangularPath <: ClosedPath
    corner::Loc
    dx::Real
    dy::Real
end
rectangular_path(corner::Loc=u0(), dx::Real=1, dy::Real=1) = RectangularPath(corner, dx, dy)
centered_rectangular_path(p, dx, dy) =
  rectangular_path(p-vxy(dx/2, dy/2), dx, dy)

# Sections also use paths, but they are centered at the origin
export rectangular_profile,
       circular_profile,
       top_aligned_rectangular_profile

rectangular_profile(Width::Real=1, Height::Real=1; width::Real=Width, height::Real=Height) =
  centered_rectangular_path(u0(), width, height)

circular_profile(Radius::Real=1; radius::Real=Radius) =
  circular_path(u0(), radius)

top_aligned_rectangular_profile(Width::Real=1, Height::Real=1; width::Real=Width, height::Real=Height) =
  rectangular_path(xy(-width/2,-height), width, height)

struct OpenPolygonalPath <: OpenPath
    vertices::Locs
end
open_polygonal_path(vertices=[u0(), x(), xy(), y()]) = OpenPolygonalPath(vertices)

struct ClosedPolygonalPath <: ClosedPath
    vertices::Locs
end
closed_polygonal_path(vertices=[u0(), x(), xy(), y()]) = ClosedPolygonalPath(ensure_no_repeated_locations(vertices))

PolygonalPath = Union{OpenPolygonalPath, ClosedPolygonalPath}

polygonal_path(vertices::Locs=[u0(), x(), xy(), y(), u0()]) =
  coincident_path_location(vertices[1], vertices[end]) ?
    closed_polygonal_path(vertices[1:end-1]) :
    open_polygonal_path(vertices)
polygonal_path(v::Loc, vs...) = polygonal_path([v, vs...])

ensure_no_repeated_locations(locs) =
    begin
        @assert (locs[1] != locs[end])
        locs
    end

path_start(path::PolygonalPath) = path.vertices[1]
path_end(path::OpenPolygonalPath) = path.vertices[end]
path_end(path::ClosedPolygonalPath) = path.vertices[1]

join_paths(p1::OpenPolygonalPath, p2::OpenPolygonalPath) =
  coincident_path_location(path_end(p1), path_start(p2)) ?
    polygonal_path([p1.vertices..., p2.vertices...]) :
    error("Paths are non contiguous")

# Splines

struct OpenSplinePath <: OpenPath
    vertices::Locs
    v0::Union{Bool,Vec}
    v1::Union{Bool,Vec}
end
open_spline_path(vertices=[u0(), x(), xy(), y()], v0=false, v1=false) =
    OpenSplinePath(vertices, v0, v1)

struct ClosedSplinePath <: ClosedPath
    vertices::Locs
end
closed_spline_path(vertices=[u0(), x(), xy(), y()]) = ClosedSplinePath(ensure_no_repeated_locations(vertices))

spline_path(vertices::Locs) =
  coincident_path_location(vertices[1], vertices[end]) ?
    closed_spline_path(vertices[1:end-1]) :
    open_spline_path(vertices)
spline_path(vs...) = spline_path(vs)

map_division(f, path) =
  f.(rotation_minimizing_frames(path_interpolated_vertices(path)))

map_division(f::Function, path::OpenSplinePath, n::Integer) =
  let interpolator = curve_interpolator(path.vertices),
      ps = map_division(interpolator, 0.0, 1.0, n),
      vts = map_division(t->Interpolations.gradient(interpolator, t)[1], 0.0, 1.0, n)
      #vns = map_division(t->Interpolations.hessian(interpolator, t)[1], 0.0, 1.0, n)
    f.(rotation_minimizing_frames(
         path.vertices[1],
         map(p->xyz(p[1], p[2], p[3], world_cs), ps),
         map(vy->vxyz(vy[1], vy[2], vy[3], world_cs), vts)))
  end

map_division(f::Function, path::ClosedSplinePath, n::Integer) =
  let interpolator = curve_interpolator(path.vertices),
      ps = map_division(interpolator, 0.0, 1.0, n, false),
      vts = map_division(t->Interpolations.gradient(interpolator, t)[1], 0.0, 1.0, n, false)
      #vns = map_division(t->Interpolations.hessian(interpolator, t)[1], 0.0, 1.0, n)
    f.(rotation_minimizing_frames(
         path.vertices[1],
         map(p->xyz(p[1], p[2], p[3], world_cs), ps),
         map(vy->vxyz(vy[1], vy[2], vy[3], world_cs), vts)))
  end

SplinePath = Union{OpenSplinePath, ClosedSplinePath}

path_start(path::SplinePath) = path.vertices[1]
path_end(path::OpenSplinePath) = path.vertices[end]
path_end(path::ClosedSplinePath) = path.vertices[1]

convert(::Type{ClosedSplinePath}, path::Path) =
  closed_spline_path(path_vertices(path))

#=
      fixed_normal(vn, vt) = norm(vn) < path_tolerance() ? SVector{3}(vpol(1, sph_phi(xyz(vt[1],vt[2],vt[3], world_cs))+pi/2).raw[1:3]) : vn
    map_division(
      t-> let p = interpolator(t),
              vt = Interpolations.gradient(interpolator, t)[1],
              vn = Interpolations.hessian(interpolator, t)[1],
              vy = cross(vt, vn)
            f(loc_from_o_vx_vy(
                xyz(p[1], p[2], p[3], world_cs),
                vxyz(vn[1], vn[2], vn[3], world_cs),
                vxyz(vy[1], vy[2], vy[3], world_cs)))
          end,
     0.0, 1.0, n)
 end
=#

# There is a set of operations over Paths:
# 1. translate a path a given vector
# 2. stroke a path
# 3. fill a (presumably closed) path
# 4. compute a path location given a length from the path beginning
# 5. compute a sub path from a path, a length from the path begining and a length increment
# 6. Produce the meta representation of a path

translate(path::CircularPath, v::Vec) = circular_path(path.center + v, path.radius)
translate(path::ArcPath, v::Vec) = arc_path(path.center + v, path.radius, path.start_angle, path.amplitude)
translate(path::RectangularPath, v::Vec) = rectangular_path(path.corner + v, path.dx, path.dy)
translate(path::OpenPolygonalPath, v::Vec) = open_polygonal_path(translate(path.vertices, v))
translate(path::ClosedPolygonalPath, v::Vec) = closed_polygonal_path(translate(path.vertices, v))
translate(path::OpenSplinePath, v::Vec) = open_spline_path(translate(path.vertices, v), path.v0, path.v1)
translate(path::ClosedSplinePath, v::Vec) = closed_spline_path(translate(path.vertices, v))
translate(ps::Locs, v::Vec) = map(p->p+v, ps)

scale(path::RectangularPath, s::Real, p::Loc=u0()) =
  let v = path.corner - p
    rectangular_path(p + v*s, path.dx*s, path.dy*s)
  end

path_length(path::CircularPath) = 2*pi*path.radius
path_length(path::ArcPath) = path.radius*abs(path.amplitude)
path_length(path::RectangularPath) = 2*(path.dx + path.dy)
path_length(path::OpenPolygonalPath) = path_length(path.vertices)
path_length(path::ClosedPolygonalPath) = path_length(path.vertices) + distance(path.vertices[end], path.vertices[1])
path_length(ps::Locs) =
  let p = ps[1]
      l = 0.0
    for i in 2:length(ps)
        pp = ps[i]
        l += distance(p, pp)
        p = pp
    end
    l
  end

location_at_length(path::CircularPath, d::Real) =
  loc_from_o_phi(path.center + vpol(path.radius, d/path.radius), d/path.radius + pi/2)
location_at_length(path::ArcPath, d::Real) =
  let Δα = d/path.radius,
      s = sign(path.amplitude)
    Δα <= abs(path.amplitude) + path_tolerance() ?
      loc_from_o_phi(path.center + vpol(path.radius, path.start_angle + Δα*s),
                     path.start_angle + (Δα + pi/2)*s) :
      error("Exceeded path length by ", Δα - path.amplitude)
  end
location_at_length(path::RectangularPath, d::Real) =
    let d = d % (2*(path.dx + path.dy)) # remove multiple periods
        p = path.corner
        for (delta, phi) in zip([path.dx, path.dy, path.dx, path.dy], [0, pi/2, pi, 3pi/2])
            if d - delta < path_tolerance()
                return loc_from_o_phi(add_pol(p, d, phi), phi)
            else
                p = add_pol(p, delta, phi)
                d -= delta
            end
        end
    end
location_at_length(path::OpenPolygonalPath, d::Real) =
    let p = path.vertices[1]
        for i in 2:length(path.vertices)
            pp = path.vertices[i]
            delta = distance(p, pp)
            if d - delta < path_tolerance()
                phi = pol_phi(pp - p)
                return loc_from_o_phi(add_pol(p, d, phi), phi)
            else
                p = pp
                d -= delta
            end
        end
        error("Exceeded path length by ", d)
    end
location_at_length(path::ClosedPolygonalPath, d::Real) =
  let p = path.vertices[1]
    for i in countfrom(1)
      pp = path.vertices[i%length(path.vertices)+1]
      delta = distance(p, pp)
      if d - delta < path_tolerance()
        phi = pol_phi(pp - p)
        return loc_from_o_phi(add_pol(p, d, phi), phi)
      else
        p = pp
        d -= delta
      end
    end
  end

subpath_starting_at(path::Path, d::Real) = subpath(path, d, path_length(path))
subpath_ending_at(path::Path, d::Real) = subpath(path, 0, d)

subpath(path::CircularPath, a::Real, b::Real) =
  arc_path(path.center, path.radius, a/path.radius, (b-a)/path.radius)
subpath(path::ArcPath, a::Real, b::Real) =
  b <= path_length(path) + path_tolerance() ?
      arc_path(path.center,
               path.radius,
               path.start_angle + a/path.radius*sign(path.amplitude),
               (b - a)/path.radius*sign(path.amplitude)) :
    error("Exceeded path length by ", path.amplitude - b/path.radius)
subpath(path::RectangularPath, a::Real, b::Real) =
  subpath(convert(ClosedPolygonalPath, path), a, b)
subpath(path::ClosedPolygonalPath, a::Real, b::Real) =
  subpath(convert(OpenPolygonalPath, path), a, b)
subpath(path::OpenPolygonalPath, a::Real, b::Real) =
  subpath_starting_at(subpath_ending_at(path, b), a)

subpath_starting_at(path::OpenPolygonalPath, d::Real) =
    let pts = path.vertices
        p1 = pts[1]
        for i in 2:length(pts)
            p2 = pts[i]
            delta = distance(p1, p2)
            diff = d - delta
            if diff < 0
                mp = p1 + (p2 - p1)*d/delta
                return open_polygonal_path([mp, pts[i:end]...])
            elseif abs(diff) < 1e-9
                return open_polygonal_path(i == length(pts) ? [pts[i], pts[i]] : pts[i:end])
            else
                p1 = p2
                d -= delta
            end
        end
        abs(d) < path_tolerance() ?
          path :
          error("Exceeded path length by ", d)
    end

subpath_ending_at(path::OpenPolygonalPath, d::Real) =
    let pts = path.vertices
        p1 = pts[1]
        for i in 2:length(pts)
            p2 = pts[i]
            delta = distance(p1, p2)
            diff = d - delta
            if diff < 0
                mp = p1 + (p2 - p1)*d/delta
                return open_polygonal_path([pts[1:i-1]..., mp])
            elseif abs(diff) < path_tolerance()
                return open_polygonal_path(pts[1:i])
            else
              p1 = p2
              d -= delta
            end
        end
        abs(d) < path_tolerance() ?
          path :
          error("Exceeded path length by ", d)
    end

#=
collect_vertices_length(p::Loc, vs::Locs, d::Real) =
    let pp = vs[1]
        dist = distance(p, pp)
        d <= dist ?
        push!(Vector{Loc}(), p + (pp - p)/dist*d) :
        unshift!(collect_vertices_length(pp, vs[2:end], d - dist), pp)
    end

subpath(path::OpenPolygonalPath, a::Real, b::Real) =
    let p = path.vertices[1]
        pts = []
        for i in 2:length(path.vertices)
            pp = path.vertices[i]
            delta = distance(p, pp)
            if a <= delta
                phi = pol_phi(pp - p)
                p0 = add_pol(p, d, phi)
                return open_polygonal_path(
                        unshift!(collect_vertices_length(p0,
                                                         vcat(path.vertices[i:end], path.vertices),
                                                         delta_d),
                                 p0))
            else
                p = pp
                d -= delta
            end
        end
        error("Exceeded path length")
    end

=#

meta_program(p::OpenPolygonalPath) =
    Expr(:call, :open_polygonal_path, meta_program(p.vertices))



# Path can be made of subparts
abstract type PathOp end
#struct MoveToOp <: PathOp loc::Loc end
#struct MoveOp <: PathOp vec::Vec end
#struct LineToOp <: PathOp loc::Loc end
struct LineOp <: PathOp
  vec::Vec
end
#struct CloseOp <: PathOp end
struct ArcOp <: PathOp
    radius::Real
    start_angle::Real
    amplitude::Real
end
#struct LineToXThenToYOp <: PathOp loc::Loc end
#struct LineToYThenToXOp <: PathOp loc::Loc end
struct LineXThenYOp <: PathOp
  vec::Vec
end
struct LineYThenXOp <: PathOp
  vec::Vec
end

struct PathOps <: Path
  start::Loc
  ops::Vector{<:PathOp}
  closed::Bool
end
open_path_ops(start, ops...) = PathOps(start, [ops...], false)
closed_path_ops(start, ops...) = PathOps(start, [ops...], true)

path_length(path::PathOps) =
  let len = mapreduce(path_length, +, path.ops, init=0)
    path.closed ?
      len + distance(path.start, location_at_length(path, len)) :
      len
  end
path_length(op::LineOp) = norm(op.vec)
path_length(op::ArcOp) = op.radius*abs(op.amplitude)

translate(path::PathOps, v::Vec) =
  PathOps(path.start + v, path.ops, path.closed)

#translate_op(op::MoveToOp, v) = MoveToOp(op.loc + v)
#translate_op(op::LineToOp, v) = LineToOp(op.loc + v)
#translate_op(op::LineToXThenToYOp, v) = LineToXThenToYOpp(op.loc + v)
#translate_op(op::LineToYThenToXOp, v) = LineToYThenToXOp(op.loc + v)
#translate_op(op::PathOp, v) = op

location_at_length(path::PathOps, d::Real) =
  let ops = path.ops,
      start = path.start
    for op in ops
      delta = path_length(op)
      if d < delta + path_tolerance()
        return location_at_length(op, start, d)
      else
        start = location_at_length(op, start, delta)
        d -= delta
      end
    end
    error("Exceeded path length by ", d)
  end

location_at_length(op::LineOp, start::Loc, d::Real) =
  start + op.vec*d/path_length(op)
location_at_length(op::ArcOp, start::Loc, d::Real) =
  let center = start - vpol(op.radius, op.start_angle),
      a = d/op.radius*sign(op.amplitude)
    center + vpol(op.radius, op.start_angle + a)
  end

subpath_starting_at(path::PathOps, d::Real) =
  let ops = path.ops,
      start = location_at_length(path, d)
    for i in 1:length(ops)
      op = ops[i]
      delta = path_length(op)
      if d == delta
        return PathOps(start, ops[i+1:end], false)
      elseif d < delta + path_tolerance()
        op = subpath_starting_at(op, d)
        return PathOps(start, [op, ops[i+1:end]...], false)
      else
        d -= delta
      end
    end
    error("Exceeded path length by ", d)
  end
subpath_ending_at(path::PathOps, d::Real) =
  let ops = path.ops
    for i in 1:length(ops)
      op = ops[i]
      delta = path_length(op)
      if d == delta
        return PathOps(path.start, ops[1:i], false)
      elseif d < delta + path_tolerance()
        return PathOps(path.start, [ops[1:i-1]..., subpath_ending_at(op, d)], false)
      else
        d -= delta
      end
    end
    error("Exceeded path length by ", d)
  end

subpath_starting_at(pathOp::LineOp, d::Real) =
  let len = path_length(pathOp.vec)
    LineOp(pathOp.vec*(len-d)/len)
  end
subpath_starting_at(pathOp::ArcOp, d::Real) =
  let a = d/pathOp.radius
    ArcOp(pathOp.radius, pathOp.start_angle + a, pathOp.amplitude - a)
  end

subpath_ending_at(pathOp::LineOp, d::Real) =
  let len = path_length(pathOp.vec)
    LineOp(pathOp.vec*d/len)
  end
subpath_ending_at(pathOp::ArcOp, d::Real) =
  let a = d/pathOp.radius
    ArcOp(pathOp.radius, pathOp.start_angle, a)
  end

subpath(path::Path, a::Real, b::Real) =
  subpath_starting_at(subpath_ending_at(path, b), a)

# A path sequence is a sequence of paths where the next element of the sequence
# starts at the same place where the previous element ends.

struct OpenPathSequence <: OpenPath
  paths::Vector{<:Path}
end
open_path_sequence(paths...) =
  OpenPathSequence(ensure_connected_paths([paths...]))
struct ClosedPathSequence <: ClosedPath
  paths::Vector{<:Path}
end
closed_path_sequence(paths...) =
  ClosedPathSequence(ensure_connected_paths([paths...]))

PathSequence = Union{OpenPathSequence, ClosedPathSequence}

path_sequence(paths...) =
  coincident_path_location(path_start(paths[1]), path_end(paths[end])) ?
    closed_path_sequence(paths...) :
    open_path_sequence(paths...)

ensure_connected_paths(paths) = # AML: Finish this
    paths
path_length(path::PathSequence) =
  sum(map(path_length, path.paths))
location_at_length(path::PathSequence, d::Real) =
  let paths = path.paths
    for path in paths
      delta = path_length(path)
      if d <= delta + path_tolerance()
          return location_at_length(path, d)
      else
          d -= delta
      end
    end
    error("Exceeded path length by ", d)
  end
subpath_starting_at(path::PathSequence, d::Real) =
  let subpaths = path.paths
    for i in 1:length(subpaths)
      subpath = subpaths[i]
      delta = path_length(subpath)
      if d == delta
        return OpenPathSequence(ops[i+1:end])
      elseif d < delta + path_tolerance()
        subpath = subpath_starting_at(subpath, d)
        return OpenPathSequence([subpath, subpaths[i+1:end]...])
      else
        d -= delta
      end
    end
    error("Exceeded path length by ", d)
  end
subpath_ending_at(path::PathSequence, d::Real) =
  let subpaths = path.paths
    for i in 1:length(subpaths)
      subpath = subpaths[i]
      delta = path_length(subpath)
      if d == delta
        return OpenPathSequence(subpaths[1:i])
      elseif d < delta + path_tolerance()
        return OpenPathSequence([subpaths[1:i-1]..., subpath_ending_at(subpath, d)])
      else
        d -= delta
      end
    end
    error("Exceeded path length by ", d)
  end

convert(::Type{PathOps}, path::PathSequence) =
  let start = in_world(location_at_length(path, 0)),
      ops = mapfoldl(convert_to_path_ops, vcat, path.paths, init=PathOp[])
    PathOps(start, ops, coincident_path_location(start, path_end(path)))
  end
convert_to_path_ops(path::ArcPath) =
  [ArcOp(path.radius, path.start_angle, path.amplitude)]
convert_to_path_ops(path::PolygonalPath) =
  [LineOp(in_world(p)-in_world(q)) for (p,q) in zip(path.vertices[[2:end;1]], path.vertices)]
convert_to_path_ops(path::RectangularPath) =
  let p0 = in_world(path.corner),
      p1 = in_world(add_x(path.corner, path.dx)),
      p2 = in_world(add_xy(path.corner, path.dx, path.dy)),
      p3 = in_world(add_y(path.corner, path.dy))
    [LineOp(p1-p0), LineOp(p2-p1), LineOp(p3-p2), LineOp(p0-p3)]
  end

# A path set is a set of independent paths.
struct PathSet <: Path
  paths::Vector{<:Path}
end

# Should we just use tuples instead of arrays?
path_set(paths...) =
  PathSet([paths...])

# Operations on path containers
translate(path::T, v::Vec) where T<:Union{PathSequence,PathSet} =
  T(translate.(path.paths, v))

# Convertions from/to paths
convert(::Type{OpenPath}, vs::Locs) =
  open_polygonal_path(vs)
convert(::Type{ClosedPath}, vs::Locs) =
  closed_polygonal_path(vs)
convert(::Type{Path}, vs::Locs) =
  coincident_path_location(vs[1], vs[end]) ?
    closed_polygonal_path(vs[1:end-1]) :
    open_polygonal_path(vs)
convert(::Type{ClosedPath}, p::OpenPolygonalPath) =
  closed_polygonal_path(coincident_path_location(path_start(p), path_end(p)) ?
    path_vertices(p)[1:end-1] :
    path_vertices(p))
convert(::Type{ClosedPath}, p::OpenPath) =
  if isa(p.ops[end], CloseOp) || coincident_path_location(path_start(p), path_end(p))
    closed_path(p.ops)
  else
    error("Can't convert to a Closed Path: $p")
  end
convert(::Type{ClosedPath}, ops::Vector{<:PathOp}) =
  if isa(ops[end], CloseOp)
    closed_path(ops)
  else
    error("Can't convert to a Closed Path: $ops")
  end
convert(::Type{OpenPath}, ops::Vector{<:PathOp}) = open_path(ops)
convert(::Type{Path}, ops::Vector{<:PathOp}) =
  if isa(ops[end], CloseOp)
    closed_path(ops)
  else
    open_path(ops)
  end
convert(::Type{ClosedPolygonalPath}, path::RectangularPath) =
  let p = path.corner,
      dx = path.dx,
      dy = path.dy
    closed_polygonal_path([p, add_x(p, dx), add_xy(p, dx, dy), add_y(p, dy)])
  end

path_interpolated_vertices(path, t0=0, t1=path_length(path), epsilon=collinearity_tolerance(), min_recursion=1) =
  let p0 = location_at_length(path, t0),
      p1 = location_at_length(path, t1),
      tm = (t0 + t1)/2.0,
      pm = location_at_length(path, tm)
    min_recursion < 0 && collinear_points(p0, pm, p1, epsilon) ?
      [p0, pm, p1] :
      [path_interpolated_vertices(path, t0, tm, epsilon, min_recursion - 1)...,
       path_interpolated_vertices(path, tm, t1, epsilon, min_recursion - 1)[2:end]...]
  end

path_interpolated_lengths(path, t0=0, t1=path_length(path), epsilon=collinearity_tolerance(), min_recursion=1) =
  let p0 = location_at_length(path, t0),
      p1 = location_at_length(path, t1),
      tm = (t0 + t1)/2.0,
      pm = location_at_length(path, tm)
    min_recursion < 0 && collinear_points(p0, pm, p1, epsilon) ?
      [t0, tm, t1] :
      [path_interpolated_lengths(path, t0, tm, epsilon, min_recursion - 1)...,
       path_interpolated_lengths(path, tm, t1, epsilon, min_recursion - 1)[2:end]...]
  end

convert(::Type{ClosedPolygonalPath}, path::CircularPath) =
  closed_polygonal_path(path_interpolated_vertices(path))
convert(::Type{OpenPolygonalPath}, path::CircularPath) =
  open_polygonal_path(path_interpolated_vertices(path))
convert(::Type{OpenPolygonalPath}, path::ArcPath) =
  open_polygonal_path(path_interpolated_vertices(path))
convert(::Type{OpenPolygonalPath}, path::ClosedPolygonalPath) =
  open_polygonal_path(vcat(path.vertices, [path.vertices[1]]))
convert(::Type{OpenPolygonalPath}, path::RectangularPath) =
  convert(OpenPolygonalPath, convert(ClosedPolygonalPath, path))
convert(::Type{OpenPolygonalPath}, path::OpenSplinePath) = # ERROR: ignores limit vectors
  let interpolator = curve_interpolator(path.vertices),
      fixed_normal(vn, vt) = norm(vn) < path_tolerance() ? SVector{3}(vpol(1, sph_phi(xyz(vt[1],vt[2],vt[3], world_cs))+pi/2).raw[1:3]) : vn
    open_polygonal_path(
      map_division(
        t-> let p = interpolator(t),
                vt = Interpolations.gradient(interpolator, t)[1],
                vn = fixed_normal(Interpolations.hessian(interpolator, t)[1], vt)
                vy = cross(vt, vn)
              loc_from_o_vx_vy(
                xyz(p[1], p[2], p[3], world_cs),
                vxyz(vn[1], vn[2], vn[3], world_cs),
                vxyz(vy[1], vy[2], vy[3], world_cs))
            end,
        0.0, 1.0, 64)) # HACK this must be parameterized!
    end

curve_interpolator(pts::Locs) =
    let pts = map(pts) do p
                let v = in_world(p).raw
                  SVector{3,Float64}(v[1], v[2], v[3])
                end
              end
        Interpolations.scale(
            interpolate(pts, BSpline(Cubic(Natural(OnGrid())))),
            range(0,stop=1,length=size(pts, 1)))
    end


convert(::Type{ClosedPolygonalPath}, path::ClosedPathSequence) =
  let paths = path.paths,
      vertices = []
    for path in paths
      append!(vertices, convert(OpenPolygonalPath, path).vertices[1:end-1])
    end
    closed_polygonal_path(vertices)
  end

convert(::Type{OpenPolygonalPath}, path::OpenPathSequence) =
  let paths = path.paths,
      vertices = []
    for path in paths[1:end-1]
      append!(vertices, convert(OpenPolygonalPath, path).vertices[1:end-1])
    end
    append!(vertices, convert(OpenPolygonalPath, paths[end]).vertices)
    open_polygonal_path(vertices)
  end

# It is possible to convert a PathSet to a singleton path
# by considering the first path as the outer path and all the
# others as inner paths
convert(::Type{ClosedPath}, pset::PathSet) =
  foldl(subtract_paths, pset.paths)

#### Utilities
path_vertices(path::OpenPolygonalPath) = path.vertices
path_vertices(path::ClosedPolygonalPath) = path.vertices
path_vertices(path::OpenSplinePath) = path.vertices
path_vertices(path::ClosedSplinePath) = path.vertices
path_vertices(path::Path) = path_vertices(convert(ClosedPolygonalPath, path))

subpaths(path::OpenPolygonalPath) =
  let ps = path.vertices
    map((p0, p1)->open_polygonal_path([p0, p1]), ps[1:end-1], ps[2:end])
  end
subpaths(path::ClosedPolygonalPath) =
  let ps = path.vertices
      map((p0, p1)->open_polygonal_path([p0, p1]), ps, [ps[2:end]..., ps[1]])
  end
subpaths(path::Path) = subpaths(convert(OpenPolygonalPath, path))

subtract_paths(path1::Path, path2::Path) =
  subtract_paths(convert(ClosedPolygonalPath, path1), convert(ClosedPolygonalPath, path2))
subtract_paths(path1::ClosedPolygonalPath, path2::ClosedPolygonalPath) =
  closed_polygonal_path(
    subtract_polygon_vertices(
      path_vertices(path1),
      path_vertices(path2)))

## Smoothnesss
# A smooth curve means that the curve is differentiable up to an intended order
# In practice, it means that it does not have "corners"
is_smooth_path(path::Path) = false
is_smooth_path(pts::Locs) = false

is_smooth_path(path::Union{ArcPath, CircularPath, OpenSplinePath, ClosedSplinePath}) = true
