using IntervalSets
using Interpolations

export Shape,
       Path,
       backend,
       current_backend,
       switch_to_backend,
       void_ref,
       delete_shape,
       delete_shapes,
       delete_all_shapes,
       set_length_unit,
       collecting_shapes,
       surface_boundary,
       curve_domain,
       surface_domain,
       create_layer,
       current_layer,
       create_block,
       instantiate_block,
       reset_backend,
       connection,
       immediate_mode,
       Backend,
       @deffamily,
       @defproxy,
       dimension,
       force_creation,
       subpath,
       subpath_starting_at,
       subpath_ending_at,
       bounding_box

#Backends are types parameterized by a key identifying the backend (e.g., AutoCAD) and by the type of reference they use
abstract type Backend{K,R} end
Base.show(io::IO, b::Backend{K,R}) where {K,R} = print(io, "Backend($(K))")

#References can be (single) native references or union or substraction of References
#Unions and subtractions are needed because actual backends frequently fail those operations
abstract type GenericRef{K,T} end

struct EmptyRef{K,T} <: GenericRef{K,T} end
struct UniversalRef{K,T} <: GenericRef{K,T} end

struct NativeRef{K,T} <: GenericRef{K,T}
  value::T
end
struct UnionRef{K,T} <: GenericRef{K,T}
  values::Tuple{Vararg{GenericRef{K,T}}}
end
struct SubtractionRef{K,T} <: GenericRef{K,T}
  value::GenericRef{K,T}
  values::Tuple{Vararg{GenericRef{K,T}}}
end

ensure_ref(b::Backend{K,T}, v::GenericRef{K,T}) where {K,T} = v
ensure_ref(b::Backend{K,T}, v::T) where {K,T} = NativeRef{K,T}(v)
ensure_ref(b::Backend{K,T}, v::Vector{T}) where {K,T} =
  length(v) == 1 ?
    NativeRef{K,T}(v[1]) :
    UnionRef{K,T}(([NativeRef{K,T}(vi) for vi in v]...,))

# currying
map_ref(b::Backend{K,T}, f::Function) where {K,T} = r -> map_ref(b, f, r)

map_ref(b::Backend{K,T}, f::Function, r::NativeRef{K,T}) where {K,T} = ensure_ref(b, f(r.value))
map_ref(b::Backend{K,T}, f::Function, r::UnionRef{K,T}) where {K,T} = UnionRef{K,T}(map(map_ref(b, f), r.values))
map_ref(b::Backend{K,T}, f::Function, r::SubtractionRef{K,T}) where {K,T} = SubtractionRef{K,T}(map_ref(b, f, r.value), map(map_ref(b, f), r.values))

# currying
collect_ref(b::Backend{K,T}) where {K,T} = r -> collect_ref(b, r)

collect_ref(b::Backend{K,T}, r::NativeRef{K,T}) where {K,T} = [r.value]
collect_ref(b::Backend{K,T}, r::UnionRef{K,T}) where {K,T} = mapreduce(collect_ref(b), vcat, r.values, init=[])
collect_ref(b::Backend{K,T}, r::SubtractionRef{K,T}) where {K,T} = vcat(collect_ref(b, r.value), mapreduce(collect_ref(b), vcat, r.values, init=[]))

# Boolean algebra laws
unite_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UniversalRef{K,T}) where {K,T} = r1
unite_ref(b::Backend{K,T}, r0::UniversalRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r0

#To avoid ambiguity
unite_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::UnionRef{K,T}) where {K,T} =
  unite_ref(b, unite_refs(b, r0), unite_refs(b, r1))
unite_ref(b::Backend{K,T}, r0::EmptyRef{K,T}, r1::EmptyRef{K,T}) where {K,T} = r0
unite_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::EmptyRef{K,T}) where {K,T} = r0
unite_ref(b::Backend{K,T}, r0::EmptyRef{K,T}, r1::UnionRef{K,T}) where {K,T} = r1
unite_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::EmptyRef{K,T}) where {K,T} = r0
unite_ref(b::Backend{K,T}, r0::EmptyRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r1

unite_refs(b::Backend{K,T}, r::UnionRef{K,T}) where {K,T} =
  foldr((r0,r1)->unite_ref(b,r0,r1), r.values, init=EmptyRef{K,T}())
unite_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::GenericRef{K,T}) where {K,T} =
  unite_ref(b, unite_refs(b, r0), r1)
unite_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UnionRef{K,T}) where {K,T} =
  unite_ref(b, r0, unite_refs(b, r1))

intersect_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UniversalRef{K,T}) where {K,T} = r0
intersect_ref(b::Backend{K,T}, r0::UniversalRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r1
intersect_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::EmptyRef{K,T}) where {K,T} = r1
intersect_ref(b::Backend{K,T}, r0::EmptyRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r0
intersect_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UnionRef{K,T}) where {K,T} =
  intersect_ref(b, r0, unite_refs(b, r1))
intersect_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::GenericRef{K,T}) where {K,T} =
  intersect_ref(b, unite_refs(b, r0), r1)

#To avoid ambiguity
subtract_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::UnionRef{K,T}) where {K,T} =
  subtract_ref(b, unite_refs(b, r0), unite_refs(b, r1))
subtract_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UniversalRef{K,T}) where {K,T} = EmptyRef{K,T}()
subtract_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::EmptyRef{K,T}) where {K,T} = r0
subtract_ref(b::Backend{K,T}, r0::EmptyRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r0
subtract_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UnionRef{K,T}) where {K,T} =
  subtract_ref(b, r0, unite_refs(b, r1))
subtract_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::GenericRef{K,T}) where {K,T} =
  subtract_ref(b, unite_refs(b, r0), r1)

# References need to be created, deleted, and recreated, depending on the way the backend works
# For example, each time a shape is consumed, it becomes deleted and might need to be recreated
mutable struct LazyRef{K,R}
  backend::Backend{K,R}
  value::GenericRef{K,R}
  created::Int
  deleted::Int
end

LazyRef(backend::Backend{K,R}) where {K,R} = LazyRef{K,R}(backend, void_ref(backend), 0, 0)
LazyRef(backend::Backend{K,R}, v::GenericRef{K,R}) where {K,R} = LazyRef{K,R}(backend, v, 1, 0)

abstract type Proxy end

backend(s::Proxy) = s.ref.backend
realized(s::Proxy) = s.ref.created == s.ref.deleted + 1
mark_deleted(s::Proxy) = realized(s) ? s.ref.deleted += 1 : error("Inconsistent creation and deletion")
ref(s::Proxy) =
  if s.ref.created == s.ref.deleted
    s.ref.value = ensure_ref(s.ref.backend, realize(s.ref.backend, s))
    s.ref.created += 1
    s.ref.value
  elseif s.ref.created == s.ref.deleted + 1
    s.ref.value
  else
    error("Inconsistent creation and deletion")
  end

# We can also use a shape as a surrogate for another shape

ensure_ref(b::Backend{K,T}, v::Proxy) where {K,T} = ref(v)



#This is a dangerous operation. I'm not sure it should exist.
set_ref!(s::Proxy, value) = s.ref.value = value

abstract type Shape <: Proxy end
Shapes = Vector{<:Shape}

map_ref(f::Function, s::Shape) = map_ref(s.ref.backend, f, ref(s))
collect_ref(s::Shape) = collect_ref(s.ref.backend, ref(s))
collect_ref(ss::Shapes) = mapreduce(collect_ref, vcat, ss, init=[])


immediate_mode = Parameter(true)
in_shape_collection = Parameter(false)
collected_shapes = Parameter(Shape[])
collecting_shapes(fn) =
    with(collected_shapes, Shape[]) do
        with(in_shape_collection, true) do
            fn()
        end
        collected_shapes()
    end

create(s::Shape) =
    begin
        immediate_mode() && ref(s)
        in_shape_collection() && push!(collected_shapes(), s)
        s
    end

force_creation(s::Shape) =
    begin
        ref(s)
        s
    end

replace_in(expr::Expr, replacements) =
    if expr.head == :.
        Expr(expr.head,
             replace_in(expr.args[1], replacements), expr.args[2])
    elseif expr.head == :quote
        expr
    else
        Expr(expr.head,
             map(arg -> replace_in(arg, replacements), expr.args) ...)
    end
replace_in(expr::Symbol, replacements) =
    get(replacements, expr, esc(expr))
replace_in(expr::Any, replacements) =
    expr

showit(s, a) = begin
    print(s)
    print(":")
    println(a)
    a
end

# The undefined backend
struct Undefined_Backend <: Backend{Int,Int} end
connection(b::Undefined_Backend) = throw(UndefinedBackendException())
void_ref(b::Undefined_Backend) = EmptyRef{Int,Int}()
const current_backend = Parameter{Backend}(Undefined_Backend())

# Side-effect full operations need to have a backend selected and will generate an exception if there is none

struct UndefinedBackendException <: Exception end
Base.show(io::IO, e::UndefinedBackendException) = print(io, "No current backend.")

# Many functions default the backend to the current_backend and throw an error if there is none.
# We will simplify their definition with a macro:
# @defop delete_all_shapes()
# that expands into
# delete_all_shapes(backend::Backend=current_backend()) = throw(UndefinedBackendException())
# Note that according to Julia semantics the previous definition actually generates two different ones:
# delete_all_shapes() = delete_all_shapes(current_backend())
# delete_all_shapes(backend::Backend) = throw(UndefinedBackendException())
# Hopefully, backends will specially the function for each specific backend

macro defop(name_params)
    name, params = name_params.args[1], name_params.args[2:end]
    quote
        export $(esc(name))
        $(esc(name))($(map(esc,params)...), backend::Backend=current_backend()) =
            throw(UndefinedBackendException())
    end
end

backend(backend::Backend) = switch_to_backend(current_backend(), backend)

switch_to_backend(from::Backend, to::Backend) = current_backend(to)

@defop current_backend_name()
@defop delete_all_shapes()
@defop set_length_unit(unit::String)
@defop reset_backend()
@defop save_as(pathname::String, format::String)

macro defproxy(name, parent, fields...)
  name_str = string(name)
  struct_name = esc(Symbol(string(map(uppercasefirst,split(name_str,'_'))...)))
  field_names = map(field -> field.args[1].args[1], fields)
  field_types = map(field -> field.args[1].args[2], fields)
  field_inits = map(field -> field.args[2], fields)
  field_renames = map(esc ∘ Symbol ∘ uppercasefirst ∘ string, field_names)
  field_replacements = Dict(zip(field_names, field_renames))
  struct_fields = map((name,typ) -> :($(name) :: $(typ)), field_names, field_types)
#  opt_params = map((name,typ,init) -> :($(name) :: $(typ) = $(init)), field_renames, field_types, field_inits)
#  key_params = map((name,typ,rename) -> :($(name) :: $(typ) = $(rename)), field_names, field_types, field_renames)
#  mk_param(name,typ) = Expr(:kw, Expr(:(::), name, typ))
  mk_param(name,typ,init) = Expr(:kw, Expr(:(::), name, typ), init)
  opt_params = map(mk_param, field_renames, field_types, map(init -> replace_in(init, field_replacements), field_inits))
  key_params = map(mk_param, field_names, field_types, field_renames)
  constructor_name = esc(name)
  predicate_name = esc(Symbol("is_", name_str))
  selector_names = map(field_name -> esc(Symbol(name_str, "_", string(field_name))), field_names)
  quote
    export $(constructor_name), $(struct_name), $(predicate_name), $(selector_names...)
    struct $struct_name <: $parent
      ref::LazyRef
      $(struct_fields...)
    end
    $(constructor_name)($(opt_params...); $(key_params...), backend::Backend=current_backend(), ref::LazyRef=LazyRef(backend)) =
      create($(struct_name)(ref, $(field_names...)))
    $(predicate_name)(v::$(struct_name)) = true
    $(predicate_name)(v::Any) = false
#    $(map((selector_name, field_name) -> :($(selector_name)(v::$(struct_name)) = v.$(field_name)),
#          selector_names, field_names)...)
    Khepri.meta_program(v::$(struct_name)) =
        Expr(:call, $(Expr(:quote, name)), $(map(field_name -> :(meta_program(v.$(field_name))), field_names)...))
  end
end

abstract type Shape0D <: Shape end
abstract type Shape1D <: Shape end
abstract type Shape2D <: Shape end
abstract type Shape3D <: Shape end

is_curve(s::Shape) = false
is_surface(s::Shape) = false
is_solid(s::Shape) = false

is_curve(s::Shape1D) = true
is_surface(s::Shape2D) = true
is_solid(s::Shape3D) = true

Shapes1D = Vector{<:Any}

@defproxy(empty_shape, Shape0D)
@defproxy(universal_shape, Shape3D)
@defproxy(point, Shape0D, position::Loc=u0())
@defproxy(line, Shape1D, vertices::Locs=[u0(), ux()])
line(v0, v1, vs...) = line([v0, v1, vs...])
@defproxy(closed_line, Shape1D, vertices::Locs=[u0(), ux(), uy()])
closed_line(v0, v1, vs...) = closed_line([v0, v1, vs...])
@defproxy(spline, Shape1D, points::Locs=[u0(), ux(), uy()], v0::Union{Bool,Vec}=false, v1::Union{Bool,Vec}=false,
          interpolator::Any=LazyParameter(Any, () -> curve_interpolator(points)))
curve_interpolator(pts::Locs) =
    let pts = map(p -> in_world(p).raw, pts)
        Interpolations.scale(
            interpolate(pts,
                        BSpline(Cubic(Natural())),
                        OnGrid()),
            range(0,stop=1,length=length(pts)))
    end

# This needs to be improved to return a proper frame
evaluate(s::Spline, t::Real) = xyz(s.interpolator()[t], world_cs)

#(def-base-shape 1D-shape (spline* [pts : (Listof Loc) (list (u0) (ux) (uy))] [v0 : (U Boolean Vec) #f] [v1 : (U Boolean Vec) #f]))

@defproxy(closed_spline, Shape1D, points::Locs=[u0(), ux(), uy()])
closed_spline(v0, v1, vs...) = closed_spline([v0, v1, vs...])
@defproxy(circle, Shape1D, center::Loc=u0(), radius::Real=1)
@defproxy(arc, Shape1D, center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(elliptic_arc, Shape1D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(ellipse, Shape1D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1)
@defproxy(polygon, Shape1D, vertices::Locs=[u0(), ux(), uy()])
polygon(v0, v1, vs...) = polygon([v0, v1, vs...])
@defproxy(regular_polygon, Shape1D, edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, inscribed::Bool=false)
@defproxy(rectangle, Shape1D, c::Loc=u0(), dx::Real=1, dy::Real=1)
rectangle(p, q) =
  let v = in_cs(q - p, p.cs)
    rectangle(p, v.x, v.y)
  end
@defproxy(surface_circle, Shape2D, center::Loc=u0(), radius::Real=1)
@defproxy(surface_arc, Shape2D, center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(surface_elliptic_arc, Shape2D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(surface_ellipse, Shape2D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1)
@defproxy(surface_polygon, Shape2D, vertices::Locs=[u0(), ux(), uy()])
surface_polygon(v0, v1, vs...) = surface_polygon([v0, v1, vs...])
@defproxy(surface_regular_polygon, Shape2D, edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, inscribed::Bool=false)
@defproxy(surface_rectangle, Shape2D, c::Loc=u0(), dx::Real=1, dy::Real=1)
@defproxy(surface, Shape2D, frontier::Shapes1D=[circle()])
surface(c0::Shape, cs...) = surface([c0, cs...])
#To be removed
surface_from = surface

surface_boundary(s::Shape2D, backend::Backend=current_backend()) =
    backend_surface_boundary(backend, s)

curve_domain(s::Shape1D, backend::Backend=current_backend()) =
    backend_curve_domain(backend, s)
map_division(f::Function, s::Shape1D, n::Int, backend::Backend=current_backend()) =
    backend_map_division(backend, f, s, n)


surface_domain(s::Shape2D, backend::Backend=current_backend()) =
    backend_surface_domain(backend, s)
map_division(f::Function, s::Shape2D, nu::Int, nv::Int, backend::Backend=current_backend()) =
    backend_map_division(backend, f, s, nu, nv)

@defproxy(text, Shape0D, str::String="", c::Loc=u0(), h::Real=1)

export text_centered
text_centered(str::String="", c::Loc=u0(), h::Real=1) =
  text(str, add_xy(c, -length(str)*h*0.85/2, -h/2), h)

@defproxy(sphere, Shape3D, center::Loc=u0(), radius::Real=1)
@defproxy(torus, Shape3D, center::Loc=u0(), re::Real=1, ri::Real=1/2)
@defproxy(cuboid, Shape3D,
  b0::Loc=u0(),        b1::Loc=add_x(b0,1), b2::Loc=add_y(b1,1), b3::Loc=add_x(b2,-1),
  t0::Loc=add_z(b0,1), t1::Loc=add_x(t0,1), t2::Loc=add_y(t1,1), t3::Loc=add_x(t2,-1))

@defproxy(regular_pyramid_frustum, Shape3D, edges::Integer=4, cb::Loc=u0(), rb::Real=1, angle::Real=0, h::Real=1, rt::Real=1, inscribed::Bool=false)
regular_pyramid_frustum(edges::Integer, cb::Loc, rb::Real, angle::Real, ct::Loc, rt::Real=1, inscribed::Bool=false) =
  let (c, h) = position_and_height(cb, ct)
    regular_pyramid_frustum(edges, c, rb, angle, h, rt, inscribed)
  end

@defproxy(regular_pyramid, Shape3D, edges::Integer=3, cb::Loc=u0(), rb::Real=1, angle::Real=0, h::Real=1, inscribed::Bool=false)
regular_pyramid(edges::Integer, cb::Loc, rb::Real, a::Real, ct::Loc, inscribed::Bool=false) =
  let (c, h) = position_and_height(cb, ct)
    regular_pyramid(edges, c, rb, a, h, inscribed)
  end

@defproxy(irregular_pyramid_fustrum, Shape3D, cbs::Locs=[ux(), uy(), uxy()], cts::Locs=[uxz(), uyz(), uxyz()])
@defproxy(irregular_pyramid, Shape3D, cbs::Locs=[ux(), uy(), uxy()], ct::Loc=uz())

@defproxy(regular_prism, Shape3D, edges::Integer=3, cb::Loc=u0(), r::Real=1, angle::Real=0, h::Real=1, inscribed::Bool=false)
regular_prism(edges::Integer, cb::Loc, r::Real, angle::Real, ct::Loc, inscribed::Bool=false) =
  let (c, h) = position_and_height(cb, ct)
    regular_prism(edges, c, r, angle, h, inscribed)
  end
@defproxy(irregular_prism, Shape3D, cbs::Locs=[ux(), uy(), uxy()], v::Vec=vz(1))
irregular_prism(cbs::Locs, h::Real) =
  irregular_prism(cbs, vz(h))

@defproxy(right_cuboid, Shape3D, cb::Loc=u0(), width::Real=1, height::Real=1, h::Real=1, angle::Real=0)
right_cuboid(cb::Loc, width::Real, height::Real, ct::Loc, angle::Real=0; backend::Backend=current_backend()) =
  let (c, h) = position_and_height(cb, ct)
    right_cuboid(c, width, height, h, angle, backend=backend)
  end
@defproxy(box, Shape3D, c::Loc=u0(), dx::Real=1, dy::Real=1, dz::Real=1)
box(c0, c1) =
  let v = in_cs(c1, c0)-c0
    box(c0, v.x, v.y, v.z)
  end
@defproxy(cone, Shape3D, cb::Loc=u0(), r::Real=1, h::Real=1)
cone(cb, r, ct) =
  let (c, h) = position_and_height(cb, ct)
    cone(c, r, h)
  end
@defproxy(cone_frustum, Shape3D, cb::Loc=u0(), rb::Real=1, h::Real=1, rt::Real=1)
cone_frustum(cb, rb, ct, rt) =
  let (c, h) = position_and_height(cb, ct)
    cone_frustum(c, rb, h, rt)
  end
@defproxy(cylinder, Shape3D, cb::Loc=u0(), r::Real=1, h::Real=1)
cylinder(cb, r, ct) =
  let (c, h) = position_and_height(cb, ct)
    cylinder(c, r, h)
  end

@defproxy(extrusion, Shape3D, profile::Shape=point(), v::Vec=vz(1))
extrusion(profile, h::Real) =
  extrusion(profile, vz(h))

@defproxy(sweep, Shape3D, path::Shape1D=circle(), profile::Shape=point(), rotation::Real=0, scale::Real=1)
@defproxy(revolve, Shape3D, profile::Shape=point(), p::Loc=u0(), n::Vec=vz(1,p.cs), start_angle::Real=0, amplitude::Real=2*pi)
@defproxy(loft, Shape3D, profiles::Shapes=[], rails::Shapes=[], ruled::Bool=false, closed::Bool=false)
loft_ruled(profiles::Shapes=[]) = loft(profiles, [], true, false)
export loft_ruled

@defproxy(move, Shape3D, shape::Shape=point(), v::Vec=vx())
@defproxy(scale, Shape3D, shape::Shape=point(), s::Real=1, p::Loc=u0())
@defproxy(rotate, Shape3D, shape::Shape=point(), angle::Real=0, p::Loc=u0(), v::Vec=vz(1,p.cs))
@defproxy(transform, Shape3D, shape::Shape=point(), xform::Loc=u0())

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
       open_polygonal_path,
       closed_polygonal_path,
       path_set,
       translate,
       stroke,
       fill,
       curve_length,
       location_at_length,
       subpath

abstract type Path end

import Base.getindex, Base.endof
getindex(p::Path, i::Real) = location_at_length(p, i)
endof(p::Path) = length(p)
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

# There is a set of operations over Paths:
# 1. translate a path a given vector
# 2. stroke a path
# 3. fill a (presumably closed) path
# 4. compute a path location given a length from the path beginning
# 5. compute a sub path from a path, a length from the path begining and a length increment
# 6. Produce the meta representation of a path

translate(path::CircularPath, v::Vec) = circular_path(path.center + v, path.radius)
translate(path::RectangularPath, v::Vec) = rectangular_path(path.corner + v, path.dx, path.dy)
translate(path::OpenPolygonalPath, v::Vec) = open_polygonal_path(map(p->p+v, path.vertices))
translate(path::ClosedPolygonalPath, v::Vec) = closed_polygonal_path(map(p->p+v, path.vertices))
translate(path::ArcPath, v::Vec) = arc_path(path.center + v, path.radius, path.start_angle, path.amplitude)


# We can also translate some shapes
translate(s::Line, v::Vec) = line(map(p -> p+v, s.vertices))
translate(s::Polygon, v::Vec) = polygon(map(p -> p+v, s.vertices))
translate(s::Circle, v::Vec) = circle(s.center+v, s.radius)
translate(s::Text, v::Vec) = text(s.str, s.c+v, s.h)

# We can translate arrays of Shapes
translate(ss::Shapes, v::Vec) = translate.(ss, v)


# We will also need to compute a bounding rectangle
bounding_rectangle(s::Union{Line, Polygon}) =
    bounding_rectangle(s.vertices)

bounding_rectangle(pts::Locs) =
    let min_p = pts[1]
        max_p = min_p
        for i in 2:length(pts)
            min_p = min_loc(min_p, pts[i])
            max_p = max_loc(max_p, pts[i])
        end
        [min_p, max_p]
    end

bounding_rectangle(ss::Shapes) =
    bounding_rectangle(mapreduce(bounding_rectangle, vcat, ss))


curve_length(path::CircularPath) = 2*pi*path.radius
curve_length(path::RectangularPath) = 2*(path.dx + path.dy)
curve_length(path::OpenPolygonalPath) = length(path.vertices)
curve_length(path::ClosedPolygonalPath) = length(path.vertices) + distance(path.vertices[end], path.vertices[1])
curve_length(s::Shape) = length(convert(Path, s))
curve_length(ps::Vector{<:Loc}) =
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
        error("Exceeded path length")
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
            if d == delta
                return open_polygonal_path(pts[i:end])
            elseif d < delta
                mp = p1 + (p2 - p1)*d/delta
                return open_polygonal_path([mp, pts[i:end]...])
            else
                p1 = p2
                d -= delta
            end
        end
        error("Exceeded path length")
    end

subpath_ending_at(path::OpenPolygonalPath, d::Real) =
    let pts = path.vertices
        p1 = pts[1]
        for i in 2:length(pts)
            p2 = pts[i]
            delta = distance(p1, p2)
            if d == delta
                return open_polygonal_path(pts[1:i])
            elseif d < delta
                mp = p1 + (p2 - p1)*d/delta
                return open_polygonal_path([pts[1:i-1]..., mp])
            else
                p1 = p2
                d -= delta
            end
        end
        error("Exceeded path length")
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
stroke(path::Path, backend::Backend=current_backend()) = backend_stroke(backend, path)
# We also need a colored stroke (and probably, something that changes line thickness)
stroke(path::Path, color::RGB, backend::Backend=current_backend()) = backend_stroke_color(backend, path, color)
# By default, we ignore the color
backend_stroke_color(backend::Backend, path::Path, color::RGB) = backend_stroke(backend, path)

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
        error("Exceeded path length")
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
        error("Exceeded path length")
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
        error("Exceeded path length")
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



# The default implementation for stroking segmented path in the backend relies on two
# dedicated functions backend_stroke_arc and backend_stroke_line
stroke(path::PathOps, backend::Backend=current_backend()) = backend_stroke(backend, path)

backend_stroke(b::Backend, path::PathOps) =
    begin
        start, curr, refs = path.start, path.start, []
        for op in path.ops
            start, curr, refs = backend_stroke_op(b, op, start, curr, refs)
        end
        if path.closed
            push!(refs, backend_stroke_line(b, [curr, start]))
        end
        backend_stroke_unite(b, refs)
    end
#=
backend_stroke_op(b::Backend, op::MoveToOp, start::Loc, curr::Loc, refs) =
    (op.loc, op.loc, refs)
backend_stroke_op(b::Backend, op::MoveOp, start::Loc, curr::Loc, refs) =
    (start, curr + op.vec, refs)
backend_stroke_op(b::Backend, op::LineToOp, start::Loc, curr::Loc, refs) =
    (start, op.loc, push!(refs, backend_stroke_line(b, [curr, op.loc])))
=#
backend_stroke_op(b::Backend, op::LineOp, start::Loc, curr::Loc, refs) =
    (start, curr + op.vec, push!(refs, backend_stroke_line(b, [curr, curr + op.vec])))
#backend_stroke_op(b::Backend, op::CloseOp, start::Loc, curr::Loc, refs) =
#    (start, start, push!(refs, backend_stroke_line(b, [curr, start])))
backend_stroke_op(b::Backend, op::ArcOp, start::Loc, curr::Loc, refs) =
    let center = curr - vpol(op.radius, op.start_angle)
        (start,
         center + vpol(op.radius, op.start_angle + op.amplitude),
         push!(refs, backend_stroke_arc(b, center, op.radius, op.start_angle, op.amplitude)))
     end
#=
backend_stroke_op(b::Backend, op::LineXThenYOp, start::Loc, curr::Loc, refs) =
    (start,
     start + op.vec,
     push!(refs, backend_stroke_line(b, [curr, curr + vec_in(op.vec, curr.cs).x, curr + op.vec])))

backend_stroke_op(b::Backend, op::LineYThenXOp, start::Loc, curr::Loc, refs) =
    (start,
     start + op.vec,
     push!(refs, backend_stroke_line(b, [curr, curr + vec_in(op.vec, curr.cs).y, curr + op.vec])))
backend_stroke_op(b::Backend, op::LineToXThenToYOp, start::Loc, curr::Loc, refs) =
    (start, op.loc, push!(refs, backend_stroke_line(b, [curr, xy(curr.x, loc_in(op.loc, curr.cs).x, curr.cs), op.loc])))
backend_stroke_op(b::Backend, op::LineToYThenToXOp, start::Loc, curr::Loc, refs) =
    (start, op.loc, push!(refs, backend_stroke_line(b, [curr, xy(curr.x, loc_in(op.loc, curr.cs).y, curr.cs), op.loc])))
=#

# A path set is a set of independent paths.

struct PathSet <: Path
    paths::Vector{<:Path}
end

# Should we just use tuples instead of arrays?
path_set(paths...) =
    PathSet([paths...])

backend_stroke(b::Backend, path::PathSet) =
    for p in path.paths
        backend_stroke(b, p)
    end

# The default implementation for filling segmented path in the backend relies on
# a dedicated function backend_fill_curves

import Base.fill
fill(path, backend=current_backend()) = backend_fill(backend, path)
backend_fill(b, path) = backend_fill_curves(b, backend_stroke(b, path))




# Convertions from/to paths
import Base.convert
convert(::Type{ClosedPath}, s::Rectangle) = closed_path([MoveToOp(s.c), RectOp(vxy(s.dx, s.dy))])
convert(::Type{ClosedPath}, s::Circle) = closed_path([CircleOp(s.center, s.radius)])
convert(::Type{Path}, s::Line) = convert(OpenPath, s.vertices)
convert(::Type{OpenPath}, vs::Locs) = open_polygonal_path(vs)
convert(::Type{ClosedPath}, vs::Locs) = closed_polygonal_path(vs)
convert(::Type{Path}, vs::Locs) =
    if vs[1] == vs[end]
        closed_polygonal_path(vs[1:end-1])
    else
        open_polygonal_path(vs)
    end

convert(::Type{ClosedPath}, s::Polygon) =
    let vs = polygon_vertices(s)
        closed_path(vcat([MoveToOp(vs[1])], map(LineToOp, vs[2:end]), [CloseOp()]))
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
convert(::Type{OpenPolygonalPath}, path::ClosedPolygonalPath) =
    open_polygonal_path(vcat(path.vertices, [path.vertices[1]]))
convert(::Type{OpenPolygonalPath}, path::RectangularPath) =
    convert(OpenPolygonalPath, convert(ClosedPolygonalPath, path))

#### Utilities

path_vertices(path::OpenPolygonalPath) = path.vertices
path_vertices(path::ClosedPolygonalPath) = path.vertices
path_vertices(path::RectangularPath) = path_vertices(convert(ClosedPolygonalPath, path))

export path_vertices

#####################################################################
export curve_domain, surface_domain, frame_at
surface_domain(s::SurfaceRectangle) = (0, s.dx, 0, s.dy)
surface_domain(s::SurfaceCircle) = (0, s.radius, 0, 2pi)
surface_domain(s::SurfaceArc) = (0, s.radius, s.start_angle, s.amplitude)


frame_at(c::Shape1D, t::Real) = backend_frame_at(backend(c), c, t)
frame_at(s::Shape2D, u::Real, v::Real) = backend_frame_at(backend(s), s, u, v)

#Some specific cases can be handled in an uniform way without the backend
frame_at(s::SurfaceRectangle, u::Real, v::Real) = add_xy(s.c, u, v)
frame_at(s::SurfaceCircle, u::Real, v::Real) = add_pol(s.center, u, v)




#####################################################################
# BIM
abstract type Measure <: Proxy end

@defproxy(level, Measure, height::Real=0.0)
create(s::Measure) = s

default_level = Parameter{Level}(level())
default_level_to_level_height = Parameter{Real}(3)
upper_level(lvl, height=default_level_to_level_height()) = level(lvl.height + height, backend=backend(lvl))

#default implementation
realize(b::Backend, s::Level) = s.height

export default_level, default_level_to_level_height

#=
@defproxy(polygonal_mass, Shape3D, points::Locs, height::Real)
@defproxy(rectangular_mass, Shape3D, center::Loc, width::Real, len::Real, height::Real)

@defproxy(column, Shape3D, center::Loc, bottom_level::Any, top_level::Any, family::Any)
=#


#=

We need to provide defaults for a lot of things. For example, we want to specify
a wall that goes through a path without having to specify the kind of wall or its
thickness and height.

This means that, apart from the wall's path, all other wall features will come
from defaults. The base height will be determined by the current level and the
wall height by the current level-to-level height. Finally, the wall thickness,
constituent parts, thermal characteristics, and so on will come from the wall
defaults.  In the case of wall, we will assume that current_wall_defaults() is
a parameter that contains a set of a wall parameters.  As an example of use, we
might have:

current_wall_defaults(wall_defaults(thickness=10))

Another option is the definition of different defaults:

thick_wall_defaults = wall_defaults(thickness=10)
thin_wall_defaults = wall_defaults(thickness=5)

which then can be make current:

current_wall_defaults(thin_wall_defaults)

In most cases, the defaults are not just one value, but a bunch of them. For a
beam, we might have:

standard_beam = beam_defaults(width=10, height=20)
current_beam_defaults(standard_beam)

Another useful feature is the ability to adapt defaults. For example:

current_beam_defaults(beam_with(standard_beam, height=20))

Finally, defaults can be created for anything. For example, in a building, we
might want to define a bunch of parameters that are relevant. The syntax is as
follows:

@defaults(building,
    width::Real=20,
    length::Real=30,
    height::Real=50)

In order to access these defaults, we can use the following:

current_building_defaults().width

In some cases, defaults are supported by the backend itself. For example, in
Revit, a wall can be specified using a family. In order to realize the wall
defaults in the current backend, we need to map from the wall parameters to the
corresponding BIM family parameters. This mapping must be described in a
different structure.

For example, a beam element might have a section with a given width and height
but, in Revit, a beam element such as "...\\Structural Framing\\Wood\\M_Timber.rfa"
has, as parameters, the dimensions b and d.  This means that we need a map, such
as Dict(:width => "b", :height => "d")))). So, for a Revit family, we might use:

RevitFamily(
    "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
    Dict(:width => "b", :height => "d"))))

However, the same beam might have a different mapping in a different backend.
This means that we need another mapping to support different backends. One
possibility is to use something similar to:

backend_family(
    revit => RevitFamily(
        "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
        Dict(:width => "b", :height => "d")),
    archicad => ArchiCADFamily(
        "BeamElement",
        Dict(:width => "size_x", :height => "size_y")),
    autocad => AutoCADFamily())

Then, we need an operation that instantiates a family. This can be done on two different
levels: (1) from a backend-specific family (e.g., RevitFamily), for example:

beam_family = RevitFamily(
    "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
    Dict(:width => "b", :height => "d"))

current_beam_defaults(beam_family_instance(beam_family, width=10, height=20)

or from a generic backend family, for example:

beam_family = backend_family(
    revit => RevitFamily(
        "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
        Dict(:width => "b", :height => "d")),
    archicad => ArchiCADFamily(
        "BeamElement"
        Dict(:width => "size_x", :height => "size_y")),
    autocad => AutoCADFamily())

current_beam_defaults(beam_family_instance(beam_family, width=10, height=20)

In this last case, the generic family will use the current_backend value to identify
which family to use.

Another important feature is the use of a delegation-based implementation for
family instances. This means that we might do

current_beam_defaults(beam_family_instance(current_beam_defaults(), width=20)

to instantiate a family that uses, by default, the same parameter values used by
another family instance.

=#

abstract type Family <: Proxy end
abstract type FamilyInstance <: Family end

family(f::Family) = f
family(f::FamilyInstance) = f.family

macro deffamily(name, parent, fields...)
  name_str = string(name)
  abstract_name = esc(Symbol(string))
  struct_name = esc(Symbol(string(map(uppercasefirst,split(name_str,'_'))...)))
  field_names = map(field -> field.args[1].args[1], fields)
  field_types = map(field -> field.args[1].args[2], fields)
  field_inits = map(field -> field.args[2], fields)
  field_renames = map(esc ∘ Symbol ∘ uppercasefirst ∘ string, field_names)
  field_replacements = Dict(zip(field_names, field_renames))
  struct_fields = map((name,typ) -> :($(name) :: $(typ)), field_names, field_types)
#  opt_params = map((name,typ,init) -> :($(name) :: $(typ) = $(init)), field_renames, field_types, field_inits)
#  key_params = map((name,typ,rename) -> :($(name) :: $(typ) = $(rename)), field_names, field_types, field_renames)
#  mk_param(name,typ) = Expr(:kw, Expr(:(::), name, typ))
  mk_param(name,typ,init) = Expr(:kw, Expr(:(::), name, typ), init)
  opt_params = map(mk_param, field_renames, field_types, map(init -> replace_in(init, field_replacements), field_inits))
  key_params = map(mk_param, field_names, field_types, field_renames)
  instance_params = map(mk_param, field_names, field_types, map(name -> :(family.$(name)), field_names))
  constructor_name = esc(name)
  instance_name = esc(Symbol(name_str, "_element")) #"_instance")) beam_family_element or beam_family_instance?
  default_name = esc(Symbol("default_", name_str))
  predicate_name = esc(Symbol("is_", name_str))
  selector_names = map(field_name -> esc(Symbol(name_str, "_", string(field_name))), field_names)
  quote
    export $(constructor_name), $(instance_name), $(default_name), $(predicate_name), $(struct_name)
    struct $struct_name <: $parent
      $(struct_fields...)
      based_on::Any #Family
      ref::Parameter{Int}
    end
    $(constructor_name)($(opt_params...);
                        $(key_params...),
                        based_on=nothing) =
      $(struct_name)($(field_names...), based_on, Parameter(-1))
    $(instance_name)(family:: Family #=$(struct_name)=#; $(instance_params...), based_on=family) =
      $(struct_name)($(field_names...), based_on, Parameter(-1))
    $(default_name) = Parameter($(constructor_name)())
    $(predicate_name)(v::$(struct_name)) = true
    $(predicate_name)(v::Any) = false
#    $(map((selector_name, field_name) -> :($(selector_name)(v::$(struct_name)) = v.$(field_name)),
#          selector_names, field_names)...)
    Khepri.meta_program(v::$(struct_name)) =
        Expr(:call, $(Expr(:quote, name)), $(map(field_name -> :(meta_program(v.$(field_name))), field_names)...))
  end
end

ref(family::Family) = family.ref()==-1 ? family.ref(backend_get_family(current_backend(), family)) : family.ref()

@deffamily(slab_family, Family,
    thickness::Real=0.2,
    coating_thickness::Real=0.0)

@defproxy(slab, Shape3D, contour::ClosedPath=rectangular_path(),
          level::Level=default_level(), family::SlabFamily=default_slab_family(),
          openings::Vector{ClosedPath}=ClosedPath[])
slab(contour; level::Level=default_level(), family::SlabFamily=default_slab_family()) =
    slab(convert(ClosedPath, contour), level, family)

# Default implementation: dispatch on the slab elements
realize(b::Backend, s::Slab) =
    realize_slab(b, s.contour, s.level, s.family)

realize_slab(b::Backend, contour::ClosedPath, level::Level, family::SlabFamily) =
    let base = vz(level.height + family.coating_thickness - family.thickness),
        thickness = family.coating_thickness + family.thickness
        # Change this to a better named protocol?
        backend_slab(b, translate(contour, base), thickness)
    end

#
export add_slab_opening
add_slab_opening(s::Slab=required(), contour::ClosedPath=circular_path()) =
    let b = backend(s)
        push!(s.openings, contour)
        if realized(s)
            set_ref!(s, realize_slab_openings(b, s, ref(s), [contour]))
        end
        s
    end

realize_slab_openings(b::Backend, s::Slab, s_ref, openings) =
    let s_base_height = s.level.height,
        s_thickness = s.family.thickness
        for opening in openings
            op_path = translate(opening, vz(s_base_height-1.1*s_thickness))
            op_ref = ensure_ref(b, backend_slab(b, op_path, s_thickness*1.2))
            s_ref = ensure_ref(b, subtract_ref(b, s_ref, op_ref))
        end
        s_ref
    end

#=
Should we eliminate this?
The rational is that an opening is not an object (like a door or a window)
However, it might be interesting to have a computational object to store properties of the opening
@defproxy(slab_opening, Shape3D, slab::Slab=required(), contour::ClosedPath=rectangular_path())
# Default implementation
realize(b::Backend, s::SlabOpening) =
    let base = vz(s.slab.level.height + s.slab.family.coating_thickness - s.slab.family.thickness - 1)
        thickness = s.slab.family.coating_thickness + s.slab.family.thickness + 1
        opening_ref = ensure_ref(b, backend_slab(b, translate(s.contour, base), thickness))
        # This is a dangerous side effect. Is this really the correct approach?
        set_ref!(s.slab, ensure_ref(b, subtract_ref(b, ref(s.slab), opening_ref)))
        opening_ref
    end
=#

# Roof

@deffamily(roof_family, Family,
    thickness::Real=0.2,
    coating_thickness::Real=0.0)

@defproxy(roof, Shape3D, contour::ClosedPath=rectangular_path(), level::Level=default_level(), family::RoofFamily=default_roof_family())

# Panel

@deffamily(panel_family, Family,
    thickness::Real=0.02)

@defproxy(panel, Shape3D, vertices::Locs=[], level::Any=default_level(), family::Any=default_panel_family())

realize(b::Backend, s::Panel) =
    let p1 = s.vertices[1],
        p2 = s.vertices[2],
        p3 = s.vertices[3],
        n = vz(s.family.thickness, cs_from_o_vx_vy(p1, p2-p1, p3-p1))
        ref(irregular_prism(map(p -> in_world(p - n), s.vertices),
                            in_world(n*2)))
    end

#=

A wall contains doors and windows

=#

# Wall

@deffamily(wall_family, Family,
    thickness::Real=0.2)

@defproxy(wall, Shape3D, path::Path=rectangular_path(),
          bottom_level::Level=default_level(),
          top_level::Level=upper_level(bottom_level),
          family::WallFamily=default_wall_family(),
          doors::Shapes=Shape[], windows::Shapes=Shape[])
wall(path::Any;
     bottom_level::Level=default_level(),
     top_level::Level=upper_level(bottom_level),
     family::WallFamily=default_wall_family()) =
    wall(convert(Path, path), bottom_level, top_level, family)
wall(p0::Loc, p1::Loc;
     bottom_level::Level=default_level(),
     top_level::Level=upper_level(bottom_level),
     family::WallFamily=default_wall_family()) =
    wall([p0, p1], bottom_level=bottom_level, top_level=top_level, family=family)

# Door

@deffamily(door_family, Family,
    width::Real=1.0,
    height::Real=2.0,
    thickness::Real=0.05)

@defproxy(door, Shape3D, wall::Wall=required(), loc::Loc=u0(), flip_x::Bool=false, flip_y::Bool=false, family::DoorFamily=default_door_family())

# Window

@deffamily(window_family, Family,
    width::Real=1.0,
    height::Real=2.0,
    thickness::Real=0.05)

@defproxy(window, Shape3D, wall::Wall=required(), loc::Loc=u0(), flip_x::Bool=false, flip_y::Bool=false, family::WindowFamily=default_window_family())

# Default implementation
realize(b::Backend, w::Wall) =
    realize_wall_openings(b, w, realize_wall_no_openings(b, w), [w.doors..., w.windows...])

realize_wall_no_openings(b::Backend, w::Wall) =
    let w_base_height = w.bottom_level.height,
        w_height = w.top_level.height - w_base_height,
        w_path = translate(w.path, vz(w_base_height))
        w_thickness = w.family.thickness
        ensure_ref(b, backend_wall(b, w_path, w_height, w_thickness))
    end

realize_wall_openings(b::Backend, w::Wall, w_ref, openings) =
    let w_base_height = w.bottom_level.height,
        w_height = w.top_level.height - w_base_height,
        w_path = translate(w.path, vz(w_base_height))
        w_thickness = w.family.thickness
        for opening in openings
            w_ref = realize_wall_opening(b, w_ref, w_path, w_thickness, opening)
            realize(b, opening)
        end
        w_ref
    end

realize_wall_opening(b::Backend, w_ref, w_path, w_thickness, op) =
    let op_base_height = op.loc.y
        op_height = op.family.height
        op_thickness = op.family.thickness
        op_path = translate(subpath(w_path, op.loc.x, op.loc.x + op.family.width), vz(op_base_height))
        op_ref = ensure_ref(b, backend_wall(b, op_path, op_height, w_thickness*1.1))
        ensure_ref(b, subtract_ref(b, w_ref, op_ref))
    end

realize(b::Backend, s::Door) =
  let base_height = s.wall.bottom_level.height + s.loc.y,
      height = s.family.height,
      subpath = translate(subpath(s.wall.path, s.loc.x, s.loc.x + s.family.width), vz(base_height))
      # we emulate a door using a small wall
      backend_wall(b, subpath, height, s.family.thickness)
  end

realize(b::Backend, s::Window) =
  let base_height = s.wall.bottom_level.height + s.loc.y,
      height = s.family.height,
      subpath = translate(subpath(s.wall.path, s.loc.x, s.loc.x + s.family.width), vz(base_height))
      # we emulate a window using a small wall
      backend_wall(b, subpath, height, s.family.thickness)
  end

##

export add_door
add_door(w::Wall=required(), loc::Loc=u0(), family::DoorFamily=default_door_family()) =
    let b = backend(w)
        d = door(w, loc, family=family)
        push!(w.doors, d)
        if realized(w)
            set_ref!(w, realize_wall_openings(b, w, ref(w), [d]))
        end
        w
    end

#
export add_window
add_window(w::Wall=required(), loc::Loc=u0(), family::WindowFamily=default_window_family()) =
    let b = backend(w)
        d = window(w, loc, family=family)
        push!(w.windows, d)
        if realized(w)
            set_ref!(w, realize_wall_openings(b, w, ref(w), [d]))
        end
        w
    end

#
# We need to redefine the default method (maybe add an option to the macro to avoid defining the meta_program)
# This needs to be fixed for windows
meta_program(w::Wall) =
    if isempty(w.doors)
        Expr(:call, :wall,
             meta_program(w.path),
             meta_program(w.bottom_level),
             meta_program(w.top_level),
             meta_program(w.family))
    else
        let door = w.doors[1]
            Expr(:call, :add_door,
                 meta_program(wall(w.path, w.bottom_level, w.top_level, w.family, w.doors[2:end], w.windows)),
                 meta_program(door.loc),
                 meta_program(door.family))
        end
    end

# Beam
# Beams are mainly horizontal elements. A beam has its top axis aligned with a line defined by two points
@deffamily(beam_family, Family,
    width::Real=1.0,
    height::Real=2.0)

@defproxy(beam, Shape3D, cb::Loc=u0(), h::Real=1, angle::Real=0, family::BeamFamily=default_beam_family())
beam(cb::Loc, ct::Loc, Angle::Real=0, Family::BeamFamily=default_beam_family(); angle::Real=Angle, family::BeamFamily=Family) =
    let (c, h) = position_and_height(cb, ct)
      beam(c, h, angle, family)
    end

# Column
# Columns are mainly vertical elements. A column has its center axis aligned with a line defined by two points

@deffamily(column_family, Family,
    width::Real=1.0,
    height::Real=2.0)

@defproxy(column, Shape3D, cb::Loc=u0(), h::Real=1, angle::Real=0, family::ColumnFamily=default_column_family())
column(cb::Loc, ct::Loc, Angle::Real=0, Family::ColumnFamily=default_column_family(); angle::Real=Angle, family::ColumnFamily=Family) =
    let (c, h) = position_and_height(cb, ct)
      column(c, h, angle, family)
    end


# Tables and chairs

@deffamily(table_family, Family,
    length::Real=1.6,
    width::Real=0.9,
    height::Real=0.75,
    top_thickness::Real=0.05,
    leg_thickness::Real=0.05)

@deffamily(chair_family, Family,
    length::Real=0.4,
    width::Real=0.4,
    height::Real=1.0,
    seat_height::Real=0.5,
    thickness::Real=0.05)

@deffamily(table_chair_family, Family,
    table_family::TableFamily=default_table_family(),
    chair_family::ChairFamily=default_chair_family(),
    chairs_top::Int=1,
    chairs_bottom::Int=1,
    chairs_right::Int=2,
    chairs_left::Int=2,
    spacing::Real=0.7)

@defproxy(table, Shape3D, loc::Loc=u0(), angle::Real=0, level::Level=default_level(), family::TableFamily=default_table_family())

realize(b::Backend, s::Table) =
    backend_rectangular_table(b, add_z(s.loc, s.level.height), s.angle, s.family)

@defproxy(chair, Shape3D, loc::Loc=u0(), angle::Real=0, level::Level=default_level(), family::ChairFamily=default_chair_family())

realize(b::Backend, s::Chair) =
    backend_chair(b, add_z(s.loc, s.level.height), s.angle, s.family)

@defproxy(table_and_chairs, Shape3D, loc::Loc=u0(), angle::Real=0, level::Level=default_level(), family::TableChairFamily=default_table_chair_family())

realize(b::Backend, s::TableAndChairs) =
    backend_rectangular_table_and_chairs(b, add_z(s.loc, s.level.height), s.angle, s.family)

# Lights

@defproxy(spotlight, Shape3D, loc::Loc=z(3), dir::Vec=vz(-1), hotspot::Real=pi/4, falloff::Real=pi/3)

realize(b::Backend, s::Spotlight) =
    backend_spotlight(b, s.loc, s.dir, s.hotspot, s.falloff)

@defproxy(ieslight, Shape3D, file::String=required(), loc::Loc=z(3), dir::Vec=vz(-1), alpha::Real=0, beta::Real=0, gamma::Real=0)

realize(b::Backend, s::Ieslight) =
    backend_ieslight(b, s.file, s.loc, s.dir, s.alpha, s.beta, s.gamma)


#################################

#We need to fix this confusion between virtual and non virtual stuff

virtual = identity
#=
@defproxy(truss_node, Shape3D, p::Loc, family::Any)
@defproxy(truss_bar, Shape3D, p0::Loc, p1::Loc, angle::Real, family::Any)
=#

import Base.union
export union, intersection, subtraction

@defproxy(union_shape, Shape3D, shapes::Shapes=[])
union(shapes::Shapes) = union_shape(shapes)
union(shape::Shape, shapes...) = union_shape([shape, shapes...])

@defproxy(intersection_shape, Shape3D, shapes::Shapes=[])
intersection(shapes::Shapes) = intersection_shape(shapes)
intersection(shape::Shape, shapes...) = intersection_shape([shape, shapes...])

@defproxy(subtraction_shape2D, Shape2D, shape::Shape=surface_circle(), shapes::Shapes=[])
@defproxy(subtraction_shape3D, Shape3D, shape::Shape=surface_circle(), shapes::Shapes=[])
subtraction(shape::Shape2D, shapes...) = subtraction_shape2D(shape, [shapes...])
subtraction(shape::Shape3D, shapes...) = subtraction_shape3D(shape, [shapes...])

@defproxy(slice, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))

@defproxy(mirror, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))
@defproxy(union_mirror, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))

@defproxy(surface_grid, Shape2D, points::AbstractMatrix{<:Loc}=zeros(Loc,(2,2)), closed_u::Bool=false, closed_v::Bool=false,
          interpolator::Any=LazyParameter(Any, () -> surface_interpolator(points)))
surface_interpolator(pts::AbstractMatrix{<:Loc}) =
    let pts = map(p -> in_world(p).raw, pts)
        Interpolations.scale(
            interpolate(pts,
                        BSpline(Cubic(Natural())),
                        OnGrid()),
            range(0,stop=1,length=size(pts, 1)),
            range(0,stop=1,length=size(pts, 2)))
    end

# This needs to be improved to return a proper frame
evaluate(s::SurfaceGrid, u::Real, v::Real) = xyz(s.interpolator()[u,v], world_cs)


@defproxy(thicken, Shape3D, shape::Shape=surface_circle(), thickness::Real=1)

# Blocks

@defproxy(block, Shape, name::String="Block", shapes::Shapes = [circle()])
@defproxy(block_instance, Shape, block::Block=required(), loc::Loc=u0(), scale::Real=1.0)

################################################################################

#Backends might use different communication mechanisms, e.g., sockets, COM, RMI, etc

#We start with socket-based communication
struct Socket_Backend{K,T} <: Backend{K,T}
  connection::LazyParameter{TCPSocket}
end

connection(b::Socket_Backend{K,T}) where {K,T} = b.connection()
reset_backend(b::Socket_Backend) = reset(b.connection)

bounding_box(shape::Shape) =
  bounding_box([shape])

bounding_box(shapes::Shapes=[]) =
  if isempty(shapes)
    [u0(), u0()]
  else
    backend_bounding_box(backend(shapes[1]), shapes)
  end

delete_shape(shape::Shape) =
  delete_shapes([shape])

delete_shapes(shapes::Shapes=[]) =
  if ! isempty(shapes)
    to_delete = filter(realized, shapes)
    backend_delete_shapes(backend(shapes[1]), to_delete)
    foreach(mark_deleted, to_delete)
  end

and_delete_shape(r::Any, shape::Shape) =
  begin
    delete_shape(shape)
    r
  end

and_delete_shapes(r::Any, shapes::Shapes) =
  begin
    delete_shapes(shapes)
    r
  end

and_mark_deleted(r::Any, shape::Shape) =
    begin
        mark_deleted(shape)
        r
    end

realize_and_delete_shapes(shape::Shape, shapes::Shapes) =
    and_delete_shapes(ref(shape), shapes)

# Common implementations for realize function

realize(b::Backend, s::UnionShape) =
    unite_refs(b, map(ref, s.shapes))

realize(b::Backend, s::SubtractionShape2D) =
    subtract_ref(b, ref(s.shape), unite_refs(b, map(ref, s.shapes)))
realize(b::Backend, s::SubtractionShape3D) =
    subtract_ref(b, ref(s.shape), unite_refs(b, map(ref, s.shapes)))


realize(b::Backend, s::Loft) =
    if all(is_point, s.profiles)
      backend_loft_points(b, s.profiles, s.rails, s.ruled, s.closed)
    elseif all(is_curve, s.profiles)
      backend_loft_curves(b, s.profiles, s.rails, s.ruled, s.closed)
    elseif all(is_surface, s.profiles)
      backend_loft_surfaces(b, s.profiles, s.rails, s.ruled, s.closed)
    elseif length(s.profiles) == 2
      let (p, sh) = if is_point(s.profiles[1])
                     (s.profiles[1], s.profiles[2])
                   elseif is_point(s.profiles[2])
                     (s.profiles[2], s.profiles[1])
                   else
                     error("Cross sections are neither points nor curves nor surfaces")
                   end
        if is_curve(sh)
          backend_loft_curve_point(b, sh, p)
        elseif is_surface(sh)
          backend_loft_surface_point(b, sh, p)
        else
          error("Can't loft the shapes")
        end
      end
    else
      error("Cross sections are neither points nor curves nor surfaces")
    end

function startSketchup(port)
  ENV["ROSETTAPORT"] = port
  args = "C:\\Users\\aml\\Dropbox\\AML\\Projects\\rosetta\\sketchup\\rosetta.rb"
  println(args)
  run(`cmd /C Sketchup -RubyStartup $args`)
  #Start listening for Sketchup
  listener = listen(port)
  connection = listener.accept()
  readline(connection) == "connected" ? connection : error("Could not connect!")
end

# CAD
@defop all_shapes()
@defop all_shapes_in_layer(layer)

# BIM
@defop all_levels()
@defop all_walls()
@defop all_walls_at_level(level)

@defop disable_update()
@defop enable_update()
import Base.view
@defop view(camera::Loc, target::Loc, lens::Real)
@defop get_view()
@defop zoom_extents()
@defop view_top()

angle_of_view(size, focal_length) = 2atan(size/2focal_length)

function dolly_effect(camera, target, lens, new_camera)
  cur_dist = distance(camera, target)
  new_dist = distance(new_camera, target)
  new_lens = lens*new_dist/cur_dist
  view(new_camera, target, new_lens)
end

dolly_effect_pull_back(delta) = begin
  camera, target, lens = get_view()
  d = distance(camera, target)
  new_camera = target + (camera-target)*(d+delta)/d
  dolly_effect(camera, target, lens, new_camera)
end

@defop prompt_position(prompt::String="Select position")


"""

"""
