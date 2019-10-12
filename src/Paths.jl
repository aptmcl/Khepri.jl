# Paths are an important concept for BIM (and other things)
using IntervalSets
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
       open_spline_path,
       closed_spline_path,
       path_set,
       path_sequence,
       translate,
       stroke,
       fill,
       path_length,
       curve_length,
       location_at_length,
       path_start,
       path_end,
       subpath,
       subpaths

abstract type Path end

import Base.getindex, Base.firstindex, Base.lastindex
getindex(p::Path, i::Real) = location_at_length(p, i)
firstindex(p::Path) = 0
lastindex(p::Path) = length(p)
getindex(p::Path, i::ClosedInterval) = subpath(p, i.left, i.right)

abstract type OpenPath <: Path end
abstract type ClosedPath <: Path end

struct ArcPath <: OpenPath
    center::Loc
    radius::Real
    start_angle::Real
    amplitude::Real
end
arc_path(center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi) =
    ArcPath(center, radius, start_angle, amplitude)
struct CircularPath <: ClosedPath
    center::Loc
    radius::Real
end
circular_path(Center::Loc=u0(), Radius::Real=1; center::Loc=Center, radius::Real=Radius) = CircularPath(center, radius)
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

ensure_no_repeated_locations(locs) =
    begin
        @assert (locs[1] != locs[end])
        locs
    end

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

# There is a set of operations over Paths:
# 1. translate a path a given vector
# 2. stroke a path
# 3. fill a (presumably closed) path
# 4. compute a path location given a length from the path beginning
# 5. compute a sub path from a path, a length from the path begining and a length increment
# 6. Produce the meta representation of a path

translate(path::CircularPath, v::Vec) = circular_path(path.center + v, path.radius)
translate(path::RectangularPath, v::Vec) = rectangular_path(path.corner + v, path.dx, path.dy)
translate(path::OpenPolygonalPath, v::Vec) = open_polygonal_path(translate(path.vertices, v))
translate(path::ClosedPolygonalPath, v::Vec) = closed_polygonal_path(translate(path.vertices, v))
translate(path::ArcPath, v::Vec) = arc_path(path.center + v, path.radius, path.start_angle, path.amplitude)
translate(path::OpenSplinePath, v::Vec) = open_polygonal_path(translate(path.vertices, v), path.v0, path.v1)
translate(path::ClosedSplinePath, v::Vec) = closed_polygonal_path(translate(path.vertices, v))
translate(ps::Locs, v::Vec) = map(p->p+v, ps)


curve_length(path::CircularPath) = 2*pi*path.radius
curve_length(path::RectangularPath) = 2*(path.dx + path.dy)
curve_length(path::OpenPolygonalPath) = curve_length(path.vertices)
curve_length(path::ClosedPolygonalPath) = curve_length(path.vertices) + distance(path.vertices[end], path.vertices[1])
curve_length(ps::Locs) =
  let p = ps[1]
      l = 0.0
    for i in 2:length(ps)
        pp = ps[i]
        l += distance(p, pp)
        p = pp
    end
    l
  end

# rename to path_length?

path_length = curve_length


location_at_length(path::CircularPath, d::Real) =
    loc_from_o_phi(path.center + vpol(path.radius, d/path.radius), d/path.radius+pi/2)
location_at_length(path::RectangularPath, d::Real) =
    let d = d % (2*(path.dx + path.dy)) # remove multiple periods
        p = path.corner
        for (delta, phi) in zip([path.dx, path.dy, path.dx, path.dy], [0, pi/2, pi, 3pi/2])
            if d - delta < 0
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
            if d - delta < 0
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
            if d - delta < 0
                phi = pol_phi(pp - p)
                return loc_from_o_phi(add_pol(p, d, phi), phi)
            else
                p = pp
                d -= delta
            end
        end
    end

subpath(path::CircularPath, a::Real, b::Real) =
    arc_path(path.center, path.radius, a/path.radius, (b-a)/path.radius)
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
        abs(d) < 1e-9 ?
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
            elseif abs(diff) < 1e-9
                return open_polygonal_path(pts[1:i])
            else
              p1 = p2
              d -= delta
            end
        end
        abs(d) < 1e-9 ?
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
struct LineOp <: PathOp vec::Vec end
#struct CloseOp <: PathOp end
struct ArcOp <: PathOp
    radius::Real
    start_angle::Real
    amplitude::Real
end
#struct LineToXThenToYOp <: PathOp loc::Loc end
#struct LineToYThenToXOp <: PathOp loc::Loc end
struct LineXThenYOp <: PathOp vec::Vec end
struct LineYThenXOp <: PathOp vec::Vec end

struct PathOps <: Path
    start::Loc
    ops::Vector{<:PathOp}
    closed::Bool
end
open_path_ops(start, ops...) = PathOps(start, [ops...], false)
closed_path_ops(start, ops...) = PathOps(start, [ops...], true)

length(path::PathOps) =
    let len = mapreduce(length, +, path.ops, init=0)
        path.closed ?
        len + distance(path.start, location_at_length(path, len)) :
        len
    end
length(op::LineOp) = length(op.vec)
length(op::ArcOp) = op.radius*(op.amplitude)

translate(path::PathOps, v::Vec) =
    PathOps(path.start + v, path.ops, path.closed)

#translate_op(op::MoveToOp, v) = MoveToOp(op.loc + v)
#translate_op(op::LineToOp, v) = LineToOp(op.loc + v)
#translate_op(op::LineToXThenToYOp, v) = LineToXThenToYOpp(op.loc + v)
#translate_op(op::LineToYThenToXOp, v) = LineToYThenToXOp(op.loc + v)
#translate_op(op::PathOp, v) = op

location_at_length(path::PathOps, d::Real) =
    let ops = path.ops
        start = path.start
        for op in ops
            delta = length(op)
            if d < delta
                return location_at_length(op, start, d)
            else
                start = location_at_length(op, start, delta)
                d -= delta
            end
        end
        error("Exceeded path length by ", d)
    end

location_at_length(op::LineOp, start::Loc, d::Real) =
    start + op.vec*d/length(op)
location_at_length(op::ArcOp, start::Loc, d::Real) =
    let center = start - vpol(op.radius, op.start_angle)
        a = d/op.radius
        center + vpol(op.radius, op.start_angle + a)
    end

subpath_starting_at(path::PathOps, d::Real) =
    let ops = path.ops
        start = location_at_length(path, d)
        for i in 1:length(ops)
            op = ops[i]
            delta = length(op)
            if d == delta
                return PathOps(start, ops[i+1:end], false)
            elseif d < delta
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
            delta = length(op)
            if d == delta
                return PathOps(path.start, ops[1:i], false)
            elseif d < delta
                return PathOps(path.start, [ops[1:i-1]..., subpath_ending_at(op, d)], false)
            else
                d -= delta
            end
        end
        error("Exceeded path length by ", d)
    end

subpath_starting_at(pathOp::LineOp, d::Real) =
    let len = length(pathOp.vec)
        LineOp(pathOp.vec*(len-d)/len)
    end

subpath_starting_at(pathOp::ArcOp, d::Real) =
    let a = d/pathOp.radius
        ArcOp(pathOp.radius, pathOp.start_angle + a, pathOp.amplitude - a)
    end

subpath_ending_at(pathOp::LineOp, d::Real) =
    let len = length(pathOp.vec)
        LineOp(pathOp.vec*d/len)
    end

subpath_ending_at(pathOp::ArcOp, d::Real) =
    let a = d/pathOp.radius
        ArcOp(pathOp.radius, pathOp.start_angle, a)
    end


subpath(path::PathOps, a::Real, b::Real) =
    subpath_starting_at(subpath_ending_at(path, b), a)

# A path sequence is a sequence of paths where the next element of the sequence
# starts at the same place where the previous element ends.

struct PathSequence <: Path
    paths::Vector{<:Path}
end

path_sequence(paths...) =
    PathSequence(ensure_connected_paths([paths...]))

ensure_connected_paths(paths) = # AML: Finish this
    paths


# A path set is a set of independent paths.

struct PathSet <: Path
    paths::Vector{<:Path}
end

# Should we just use tuples instead of arrays?
path_set(paths...) =
    PathSet([paths...])

# Convertions from/to paths
import Base.convert
convert(::Type{OpenPath}, vs::Locs) = open_polygonal_path(vs)
convert(::Type{ClosedPath}, vs::Locs) = closed_polygonal_path(vs)
convert(::Type{Path}, vs::Locs) =
    if vs[1] == vs[end]
        closed_polygonal_path(vs[1:end-1])
    else
        open_polygonal_path(vs)
    end
convert(::Type{ClosedPath}, p::OpenPath) =
    if isa(p.ops[end], CloseOp) || distance(path_start(p), path_end(p)) < 1e-16 #HACK Use a global
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
    let p = path.corner
        dx = path.dx
        dy = path.dy
        closed_polygonal_path([p, add_x(p, dx), add_xy(p, dx, dy), add_y(p, dy)])
    end
# The next one is looses precision. It converts from a circular to a polygonal path
# Note that the level of detail is hardwired. Maybe it will be good to make this a parameter
convert(::Type{ClosedPolygonalPath}, path::CircularPath) =
    let c = path.center
        r = path.radius
        closed_polygonal_path([c + vpol(r, phi) for phi in division(0, 2pi, 64, false)])
    end

convert(::Type{OpenPolygonalPath}, path::ClosedPolygonalPath) =
    open_polygonal_path(vcat(path.vertices, [path.vertices[1]]))
convert(::Type{OpenPolygonalPath}, path::RectangularPath) =
    convert(OpenPolygonalPath, convert(ClosedPolygonalPath, path))

#### Utilities
export path_vertices, subpaths, subtract_paths

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


subtract_paths(path1::ClosedPolygonalPath, path2::ClosedPolygonalPath) =
  closed_polygonal_path(
    subtract_polygon_vertices(
      path_vertices(path1),
      path_vertices(path2)))

path_start(path::OpenPolygonalPath) = path.vertices[1]
path_end(path::OpenPolygonalPath) = path.vertices[end]
path_start(path::ClosedPolygonalPath) = path.vertices[1]
path_end(path::ClosedPolygonalPath) = path.vertices[1]
