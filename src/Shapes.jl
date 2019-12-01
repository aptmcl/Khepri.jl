using Interpolations

export Shape,
       Path,
       backend,
       backend_name,
       current_backend,
       has_current_backend,
       switch_to_backend,
       void_ref,
       delete_shape, delete_shapes,
       delete_all_shapes,
       set_length_unit,
       collecting_shapes,
       collected_shapes,
       with_transaction,
       surface_boundary,
       curve_domain,
       surface_domain,
       get_layer,
       create_layer,
       get_or_create_layer,
       current_layer,
       get_material,
       create_material,
       get_or_create_material,
       current_material,
       create_block,
       instantiate_block,
       reset_backend,
       connection,
       Backend,
       @deffamily,
       @defproxy,
       dimension,
       force_creation,
       subpath,
       subpath_starting_at,
       subpath_ending_at,
       bounding_box,
       capture_shape, capture_shapes,
       captured_shape, captured_shapes,
       revolve


#Backends are types parameterized by a key identifying the backend (e.g., AutoCAD) and by the type of reference they use

abstract type Backend{K,R} end
Base.show(io::IO, b::Backend{K,R}) where {K,R} = print(io, backend_name(b))

backend_name(b::Backend{K,R}) where {K,R} = "Backend"

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
ensure_ref(b::Backend{K,T}, v::Vector{<:GenericRef{K, T}}) where {K,T} =
  length(v) == 1 ?
    v[1] :
    UnionRef{K,T}((v...,))

# currying
map_ref(b::Backend{K,T}, f::Function) where {K,T} = r -> map_ref(b, f, r)

map_ref(b::Backend{K,T}, f::Function, r::NativeRef{K,T}) where {K,T} = ensure_ref(b, f(r.value))
map_ref(b::Backend{K,T}, f::Function, r::UnionRef{K,T}) where {K,T} = UnionRef{K,T}(map(map_ref(b, f), r.values))
map_ref(b::Backend{K,T}, f::Function, r::SubtractionRef{K,T}) where {K,T} = SubtractionRef{K,T}(map_ref(b, f, r.value), map(map_ref(b, f), r.values))

# currying
collect_ref(b::Backend{K,T}) where {K,T} = r -> collect_ref(b, r)

collect_ref(b::Backend{K,T}, r::EmptyRef{K,T}) where {K,T} = []
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
# This is so stupid. We need call-next-method.
really_mark_deleted(s::Proxy) = really_mark_deleted(s.ref)
really_mark_deleted(ref::LazyRef) = ref.deleted = ref.created
really_mark_deleted(s::Any) = nothing
mark_deleted(s::Proxy) = really_mark_deleted(s)
# We also need to propagate this to all dependencies
mark_deleted(ss::Array{<:Proxy}) = foreach(mark_deleted, ss)
mark_deleted(s::Any) = nothing
marked_deleted(s::Proxy) = s.ref.deleted == s.ref.created

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
Base.show(io::IO, s::Shape) =
    print(io, "$(typeof(s))(...)")

Shapes = Vector{<:Shape}

map_ref(f::Function, s::Shape) = map_ref(s.ref.backend, f, ref(s))
collect_ref(s::Shape) = collect_ref(s.ref.backend, ref(s))
collect_ref(ss::Shapes) = mapreduce(collect_ref, vcat, ss, init=[])

#=
Whenever a shape is created, it might be replaced by another shape (e.g., when
it is necessary to avoid collisions), it might be eagerly realized in its backend,
depending on the kind of shape and on the kind of backend (and/or its current state).
Another possibility is for the shape to be saved in some container.
It might also be necessary to record the control flow that caused the shape to be created.
This means that we need to control what happens immediately after a shape is initialized.
The protocol after_init takes care of that.
=#

after_init(a::Any) = a
after_init(s::Shape) = maybe_trace(maybe_collect(maybe_realize(maybe_replace(s))))

#=
Shapes might need to be replaced with other shapes. For example, when elements
join at precise locations, it might be easier, from the programming point of
view, to just duplicate the joining element. However, from the constructive
point of view, just one joining element needs to be produced.
=#

maybe_replace(s::Any) = s

#=
Backends might need to immediately realize a shape while supporting further modifications
e.g., using boolean operations. Others, however, cannot do that and can only realize
shapes by request, presumably, when they have complete information about them.
A middle term might be a backend that supports both modes.
=#

delay_realize(b::Backend, s::Shape) = s
force_realize(b::Backend, s::Shape) = (ref(s); s)

maybe_realize(s::Shape, b::Backend=backend(s)) = maybe_realize(b, s)

#=
Even if a backend is eager, it might be necessary to temporarily delay the
realization of shapes, particularly, when the construction is incremental.
=#

delaying_realize = Parameter(false)
maybe_realize(b::Backend, s::Shape) =
  delaying_realize() ?
    delay_realize(b, s) :
    force_realize(b, s)

abstract type LazyBackend{K,T} <: Backend{K,T} end
maybe_realize(b::LazyBackend, s::Shape) = delay_realize(b, s)
delay_realize(b::LazyBackend, s::Shape) = (push!(b.shapes, s); s)

with_transaction(fn) =
  maybe_realize(with(fn, delaying_realize, true))

#=
Frequently, we need to collect all shapes that are created:
=#

in_shape_collection = Parameter(false)
collected_shapes = Parameter(Shape[])
collect_shape!(s::Shape) = (push!(collected_shapes(), s); s)
collecting_shapes(fn) =
    with(collected_shapes, Shape[]) do
        with(in_shape_collection, true) do
            fn()
        end
        collected_shapes()
    end
maybe_collect(s::Shape) = (in_shape_collection() && collect_shape!(s); s)


######################################################
#Traceability
traceability = Parameter(false)
trace_depth = Parameter(1000)
excluded_modules = Parameter([Base, Base.CoreLogging, Khepri])
# We a dict from shapes to file locations
# and a dict from file locations to shapes
shape_to_file_locations = IdDict()
file_location_to_shapes = Dict()

export traceability, trace_depth, excluded_modules, clear_trace!, shape_source, source_shapes

shape_source(s) = get(shape_to_file_locations, s, [])
source_shapes(file, line) = get(file_location_to_shapes, (file, line), [])

clear_trace!() =
  begin
    empty!(shape_to_file_locations)
    empty!(file_location_to_shapes)
  end
#=
We do not care about frames that are unrelated to the application.
=#
interesting_locations(frames) =
  let locations = [],
      max_depth = min(trace_depth(), length(frames)-0)#14)
    for i in 2:max_depth
      let frame = frames[i],
          linfo = frame.linfo
        if linfo isa Core.CodeInfo ||
           (linfo isa Core.MethodInstance &&
            ! (linfo.def.module in excluded_modules()))
          push!(locations, (frame.file, frame.line))
        end
      end
    end
    locations
  end

trace!(s) =
  let frames = stacktrace(),
      locations = interesting_locations(frames)
    shape_to_file_locations[s] = locations
    for location in locations
      file_location_to_shapes[location] = Shape[get(file_location_to_shapes, location, [])..., s]
    end
    s
  end

maybe_trace(s) = (traceability() && trace!(s); s)

######################################################

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
struct UndefinedBackend <: Backend{Int,Int} end
connection(b::UndefinedBackend) = throw(UndefinedBackendException())
void_ref(b::UndefinedBackend) = EmptyRef{Int,Int}()
const undefined_backend = UndefinedBackend()

const current_backend = Parameter{Backend}(undefined_backend)
has_current_backend() = current_backend() != undefined_backend

# Side-effect full operations need to have a backend selected and will generate an exception if there is none

struct UndefinedBackendException <: Exception end
Base.show(io::IO, e::UndefinedBackendException) = print(io, "No current backend.")

realize(::UndefinedBackend, ::Shape) = throw(UndefinedBackendException())

# Many functions default the backend to the current_backend and throw an error if there is none.
# We will simplify their definition with a macro:
# @defop delete_all_shapes()
# that expands into
# delete_all_shapes(backend::Backend=current_backend()) = throw(UndefinedBackendException())
# Note that according to Julia semantics the previous definition actually generates two different ones:
# delete_all_shapes() = delete_all_shapes(current_backend())
# delete_all_shapes(backend::Backend) = throw(UndefinedBackendException())
# Hopefully, backends will specialize the function for each specific backend

macro defop(name_params)
    name, params = name_params.args[1], name_params.args[2:end]
    quote
        export $(esc(name))
        $(esc(name))($(map(esc,params)...), backend::Backend=current_backend()) =
            throw(UndefinedBackendException())
    end
end

macro defshapeop(name_params)
    name, params = name_params.args[1], name_params.args[2:end]
    quote
        export $(esc(name))
        $(esc(name))(s::Shape, $(map(esc,params)...), b::Backend=backend(s)) =
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

struct WrongTypeForParam <: Exception
  param::Symbol
  value::Any
  expected_type::Type
end
Base.showerror(io::IO, e::WrongTypeForParam) =
  print(io, "$(e.param) expected a $(e.expected_type) but got $(e.value) of type $(typeof(e.value))")


macro defproxy(name, parent, fields...)
  name_str = string(name)
  struct_name = esc(Symbol(string(map(uppercasefirst,split(name_str,'_'))...)))
  field_names = map(field -> field.args[1].args[1], fields)
  field_types = map(field -> field.args[1].args[2], fields)
  field_inits = map(field -> field.args[2], fields)
#  field_renames = map(esc ∘ Symbol ∘ uppercasefirst ∘ string, field_names)
  field_renames = map(Symbol ∘ string, field_names)
  field_replacements = Dict(zip(field_names, field_renames))
  struct_fields = map((name,typ) -> :($(name) :: $(typ)), field_names, field_types)
#  opt_params = map((name,typ,init) -> :($(name) :: $(typ) = $(init)), field_renames, field_types, field_inits)
#  key_params = map((name,typ,rename) -> :($(name) :: $(typ) = $(rename)), field_names, field_types, field_renames)
#  mk_param(name,typ) = Expr(:kw, Expr(:(::), name, typ))
  mk_param(name,typ,init) = Expr(:kw, name, init) #Expr(:kw, Expr(:(::), name, typ), init)
  opt_params = map(mk_param, field_renames, field_types, map(init -> replace_in(init, field_replacements), field_inits))
  key_params = map(mk_param, field_names, field_types, field_renames)
  constructor_name = esc(name)
  predicate_name = esc(Symbol("is_", name_str))
  #mk_convert(name,typ) = :(isa($(esc(name)), $(typ)) ? $(esc(name)) : throw(WrongTypeForParam($(QuoteNode(name)), $(esc(name)), $(typ))))
  mk_convert(name,typ) = :($(esc(name)))
  field_converts = map(mk_convert, field_names, field_types)
  selector_names = map(field_name -> esc(Symbol(name_str, "_", string(field_name))), field_names)
  quote
    export $(constructor_name), $(struct_name), $(predicate_name) #, $(selector_names...)
    struct $struct_name <: $parent
      ref::LazyRef
      $(struct_fields...)
    end
    @noinline $(constructor_name)($(opt_params...); $(key_params...), backend::Backend=current_backend(), ref::LazyRef=LazyRef(backend)) =
      after_init($(struct_name)(ref, $(field_converts...)))
    $(predicate_name)(v::$(struct_name)) = true
    $(predicate_name)(v::Any) = false
    $(map((selector_name, field_name) -> :($(selector_name)(v::$(struct_name)) = v.$(field_name)),
          selector_names, field_names)...)
    Khepri.mark_deleted(v::$(struct_name)) =
      if ! marked_deleted(v)
        really_mark_deleted(v)
        $(map(field_name -> :(mark_deleted(v.$(field_name))), field_names)...)
      end
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

# HACK: Fix element type
Shapes0D = Vector{<:Any}
Shapes1D = Vector{<:Any}
Shapes2D = Vector{<:Any}

@defproxy(empty_shape, Shape0D)
@defproxy(universal_shape, Shape3D)
@defproxy(point, Shape0D, position::Loc=u0())
@defproxy(line, Shape1D, vertices::Locs=[u0(), ux()])
line(v0::Loc, v1::Loc, vs...) = line([v0, v1, vs...])
@defproxy(closed_line, Shape1D, vertices::Locs=[u0(), ux(), uy()])
closed_line(v0::Loc, v1::Loc, vs...) = closed_line([v0, v1, vs...])
@defproxy(spline, Shape1D, points::Locs=[u0(), ux(), uy()], v0::Union{Bool,Vec}=false, v1::Union{Bool,Vec}=false,
          interpolator::Parameter{Any}=Parameter{Any}(missing))
spline(v0::Loc, v1::Loc, vs...) = spline([v0, v1, vs...])

#=
evaluate(s::Spline, t::Real) =
  let interpolator = s.interpolator
    if ismissing(interpolator())
      interpolator(curve_interpolator(s.points))
    end
    let p = interpolator()(t),
        vt = Interpolations.gradient(interpolator(), t)[1],
        vn = Interpolations.hessian(interpolator(), t)[1]
      loc_from_o_vx_vy(
        xyz(p[1], p[2], p[3], world_cs),
        vxyz(vt[1], vt[2], vt[3], world_cs),
        vxyz(vn[1], vn[2], vn[3], world_cs))
    end
  end

=#
evaluate(s::Spline, t::Real) =
  let interpolator = s.interpolator
    if ismissing(interpolator())
      interpolator(curve_interpolator(s.points))
    end
    let p = interpolator()(t),
        vt = Interpolations.gradient(interpolator(), t)[1],
        vn = Interpolations.hessian(interpolator(), t)[1],
        vy = cross(vt, vn)
      loc_from_o_vx_vy(
        xyz(p[1], p[2], p[3], world_cs),
        vxyz(vn[1], vn[2], vn[3], world_cs),
        vxyz(vy[1], vy[2], vy[3], world_cs))
    end
  end

curve_domain(s::Spline) = (0.0, 1.0)
frame_at(s::Spline, t::Real) = evaluate(s, t)
map_division(f::Function, s::Spline, n::Int, backend::Backend=backend(s)) =
  backend_map_division(backend, f, s, n)
#=HACK, THIS IS NOT READY, YET. COMPARE WITH THE BACKEND VERSION!!!!!!
  let (t1, t2) = curve_domain(s)
    map_division(t1, t2, n) do t
        f(frame_at(s, t))
    end
  end
=#
#(def-base-shape 1D-shape (spline* [pts : (Listof Loc) (list (u0) (ux) (uy))] [v0 : (U Boolean Vec) #f] [v1 : (U Boolean Vec) #f]))

@defproxy(closed_spline, Shape1D, points::Locs=[u0(), ux(), uy()])
closed_spline(v0, v1, vs...) = closed_spline([v0, v1, vs...])
@defproxy(circle, Shape1D, center::Loc=u0(), radius::Real=1)
@defproxy(arc, Shape1D, center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(elliptic_arc, Shape1D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(ellipse, Shape1D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1)
@defproxy(polygon, Shape1D, vertices::Locs=[u0(), ux(), uy()])
polygon(v0, v1, vs...) = polygon([v0, v1, vs...])
@defproxy(regular_polygon, Shape1D, edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, inscribed::Bool=true)
@defproxy(rectangle, Shape1D, corner::Loc=u0(), dx::Real=1, dy::Real=1)
rectangle(p::Loc, q::Loc) =
  let v = in_cs(q - p, p.cs)
    rectangle(p, v.x, v.y)
  end
@defproxy(surface_circle, Shape2D, center::Loc=u0(), radius::Real=1)
@defproxy(surface_arc, Shape2D, center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(surface_elliptic_arc, Shape2D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(surface_ellipse, Shape2D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1)
@defproxy(surface_polygon, Shape2D, vertices::Locs=[u0(), ux(), uy()])
surface_polygon(v0, v1, vs...) = surface_polygon([v0, v1, vs...])
@defproxy(surface_regular_polygon, Shape2D, edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, inscribed::Bool=true)
@defproxy(surface_rectangle, Shape2D, corner::Loc=u0(), dx::Real=1, dy::Real=1)
surface_rectangle(p::Loc, q::Loc) =
  let v = in_cs(q - p, p.cs)
    surface_rectangle(p, v.x, v.y)
  end
@defproxy(surface, Shape2D, frontier::Shapes1D=[circle()])
surface(c0::Shape, cs...) = surface([c0, cs...])
#To be removed
surface_from = surface

@defproxy(surface_path, Shape2D, path::ClosedPath=[circular_path()])
realize(b::Backend, s::SurfacePath) = backend_fill(b, s.path)

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

@defproxy(text, Shape0D, str::String="", corner::Loc=u0(), height::Real=1)

export text_centered
text_centered(str::String="", center::Loc=u0(), height::Real=1) =
  text(str, add_xy(center, -length(str)*height*0.85/2, -height/2), height)

# This is for unknown shapes (they are opaque, the only thing you can do with then
# might be just delete them)
@defproxy(unknown, Shape3D, baseref::Any=required())


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
regular_pyramid(edges::Integer, cb::Loc, rb::Real, angle::Real, ct::Loc, inscribed::Bool=false) =
  let (c, h) = position_and_height(cb, ct)
    regular_pyramid(edges, c, rb, angle, h, inscribed)
  end

@defproxy(irregular_pyramid_frustum, Shape3D, bs::Locs=[ux(), uy(), uxy()], ts::Locs=[uxz(), uyz(), uxyz()])
@defproxy(irregular_pyramid, Shape3D, bs::Locs=[ux(), uy(), uxy()], t::Loc=uz())

@defproxy(regular_prism, Shape3D, edges::Integer=3, cb::Loc=u0(), r::Real=1, angle::Real=0, h::Real=1, inscribed::Bool=false)
regular_prism(edges::Integer, cb::Loc, r::Real, angle::Real, ct::Loc, inscribed::Bool=false) =
  let (c, h) = position_and_height(cb, ct)
    regular_prism(edges, c, r, angle, h, inscribed)
  end
@defproxy(irregular_prism, Shape3D, bs::Locs=[ux(), uy(), uxy()], v::Vec=vz(1))
irregular_prism(bs::Locs, h::Real) =
  irregular_prism(bs, vz(h))

@defproxy(right_cuboid, Shape3D, cb::Loc=u0(), width::Real=1, height::Real=1, h::Real=1, angle::Real=0)
right_cuboid(cb::Loc, width::Real, height::Real, ct::Loc, angle::Real=0; backend::Backend=current_backend()) =
  let (c, h) = position_and_height(cb, ct)
    right_cuboid(c, width, height, h, angle, backend=backend)
  end
@defproxy(box, Shape3D, c::Loc=u0(), dx::Real=1, dy::Real=dx, dz::Real=dy)
box(c0::Loc, c1::Loc) =
  let v = in_cs(c1, c0)-c0
    box(c0, v.x, v.y, v.z)
  end
@defproxy(cone, Shape3D, cb::Loc=u0(), r::Real=1, h::Real=1)
cone(cb::Loc, r::Real, ct::Loc) =
  let (c, h) = position_and_height(cb, ct)
    cone(c, r, h)
  end
@defproxy(cone_frustum, Shape3D, cb::Loc=u0(), rb::Real=1, h::Real=1, rt::Real=1)
cone_frustum(cb::Loc, rb::Real, ct::Loc, rt::Real) =
  let (c, h) = position_and_height(cb, ct)
    cone_frustum(c, rb, h, rt)
  end
@defproxy(cylinder, Shape3D, cb::Loc=u0(), r::Real=1, h::Real=1)
cylinder(cb::Loc, r::Real, ct::Loc) =
  let (c, h) = position_and_height(cb, ct)
    cylinder(c, r, h)
  end

@defproxy(extrusion, Shape3D, profile::Shape=point(), v::Vec=vz(1))
extrusion(profile, h::Real) =
  extrusion(profile, vz(h))

realize(b::Backend, s::Extrusion) =
  backend_extrusion(backend(s), s.profile, s.v)

backend_extrusion(b::Backend, p::Point, v::Vec) =
  realize_and_delete_shapes(line([p.position, p.position + v], backend=b), [p])

@defproxy(sweep, Shape3D, path::Shape1D=circle(), profile::Shape=point(), rotation::Real=0, scale::Real=1)

realize(b::Backend, s::Sweep) =
  backend_sweep(backend(s), s.path, s.profile, s.rotation, s.scale)

@defproxy(revolve_point, Shape1D, profile::Shape0D=point(), p::Loc=u0(), n::Vec=vz(1,p.cs), start_angle::Real=0, amplitude::Real=2*pi)
@defproxy(revolve_curve, Shape2D, profile::Shape1D=line(), p::Loc=u0(), n::Vec=vz(1,p.cs), start_angle::Real=0, amplitude::Real=2*pi)
@defproxy(revolve_surface, Shape3D, profile::Shape2D=circle(), p::Loc=u0(), n::Vec=vz(1,p.cs), start_angle::Real=0, amplitude::Real=2*pi)
revolve(profile::Shape=point(x(1)), p::Loc=u0(), n::Vec=vz(1,p.cs), start_angle::Real=0, amplitude::Real=2*pi) =
  if is_point(profile)
    revolve_point(profile, p, n, start_angle, amplitude)
  elseif is_curve(profile)
    revolve_curve(profile, p, n, start_angle, amplitude)
  elseif is_surface(profile)
    revolve_surface(profile, p, n, start_angle, amplitude)
  elseif is_union_shape(profile)
    union(map(s->revolve(s, p, n, start_angle, amplitude), profile.shapes))
  elseif is_empty_shape(profile)
    profile
  else
    error("Profile is neither a point nor a curve nor a surface")
  end

realize(b::Backend, s::RevolvePoint) =
  backend_revolve_point(b, s.profile, s.p, s.n, s.start_angle, s.amplitude)
realize(b::Backend, s::RevolveCurve) =
  backend_revolve_curve(b, s.profile, s.p, s.n, s.start_angle, s.amplitude)
realize(b::Backend, s::RevolveSurface) =
  backend_revolve_surface(b, s.profile, s.p, s.n, s.start_angle, s.amplitude)

@defproxy(loft_points, Shape1D, profiles::Shapes0D=Shape[], rails::Shapes=Shape[], ruled::Bool=false, closed::Bool=false)
@defproxy(loft_curves, Shape2D, profiles::Shapes1D=Shape[], rails::Shapes=Shape[], ruled::Bool=false, closed::Bool=false)
@defproxy(loft_surfaces, Shape3D, profiles::Shapes2D=Shape[], rails::Shapes=Shape[], ruled::Bool=false, closed::Bool=false)
@defproxy(loft_curve_point, Shape2D, profile::Shape1D=circle(), point::Shape0D=point(z(1)))
@defproxy(loft_surface_point, Shape3D, profile::Shape2D=surface_circle(), point::Shapes=point(z(1)))

loft(profiles::Shapes=Shape[], rails::Shapes=Shape[], ruled::Bool=false, closed::Bool=false) =
  if all(is_point, profiles)
    loft_points(profiles, rails, ruled, closed)
  elseif all(is_curve, profiles)
    loft_curves(profiles, rails, ruled, closed)
  elseif all(is_surface, profiles)
    loft_surfaces(profiles, rails, ruled, closed)
  elseif length(profiles) == 2
    let (p, sh) = if is_point(profiles[1])
                    (profiles[1], profiles[2])
                  elseif is_point(profiles[2])
                    (profiles[2], profiles[1])
                  else
                    error("Cross sections are neither points nor curves nor surfaces")
                  end
      if is_curve(sh)
        loft_curve_point(sh, p)
      elseif is_surface(sh)
        loft_surface_point(sh, p)
      else
        error("Can't loft the shapes")
      end
    end
  else
    error("Cross sections are neither points nor curves nor surfaces")
  end

loft_ruled(profiles::Shapes=Shape[]) = loft(profiles, Shape[], true, false)
export loft, loft_ruled

realize(b::Backend, s::LoftPoints) = backend_loft_points(backend(s), s.profiles, s.rails, s.ruled, s.closed)
realize(b::Backend, s::LoftCurves) = backend_loft_curves(backend(s), s.profiles, s.rails, s.ruled, s.closed)
realize(b::Backend, s::LoftSurfaces) = backend_loft_surfaces(backend(s), s.profiles, s.rails, s.ruled, s.closed)
realize(b::Backend, s::LoftCurvePoint) = backend_loft_curve_point(backend(s), s.profile, s.point)
realize(b::Backend, s::LoftSurfacePoint) = backend_loft_surface_point(backend(s), s.profile, s.point)

backend_loft_points(b::Backend, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  let f = (ruled ? (closed ? polygon : line) : (closed ? closed_spline : spline))
    and_delete_shapes(ref(f(map(point_position, profiles), backend=b)),
                      vcat(profiles, rails))
  end

@defproxy(move, Shape3D, shape::Shape=point(), v::Vec=vx())
@defproxy(scale, Shape3D, shape::Shape=point(), s::Real=1, p::Loc=u0())
@defproxy(rotate, Shape3D, shape::Shape=point(), angle::Real=0, p::Loc=u0(), v::Vec=vz(1,p.cs))
@defproxy(transform, Shape3D, shape::Shape=point(), xform::Loc=u0())

#####################################################################

# We can also translate some shapes
translate(s::Line, v::Vec) = line(map(p -> p+v, s.vertices))
translate(s::Polygon, v::Vec) = polygon(map(p -> p+v, s.vertices))
translate(s::Circle, v::Vec) = circle(s.center+v, s.radius)
translate(s::Text, v::Vec) = text(s.str, s.c+v, s.h)

# We can translate arrays of Shapes
translate(ss::Shapes, v::Vec) = translate.(ss, v)

# We can compute the length of shapes as long as we can convert them
curve_length(s::Shape) = curve_length(convert(Path, s))

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


#####################################################################
## We might also be insterested in seeing paths
# This is more for debuging

stroke(path::Path, backend::Backend=current_backend()) = backend_stroke(backend, path)
# We also need a colored stroke (and probably, something that changes line thickness)
stroke(path::Path, color::RGB, backend::Backend=current_backend()) = backend_stroke_color(backend, path, color)
# By default, we ignore the color
backend_stroke_color(backend::Backend, path::Path, color::RGB) = backend_stroke(backend, path)

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

backend_stroke(b::Backend, path::PathSet) =
    for p in path.paths
        backend_stroke(b, p)
    end

# The default implementation for filling segmented path in the backend relies on
# a dedicated function backend_fill_curves

import Base.fill
fill(path, backend=current_backend()) = backend_fill(backend, path)
backend_fill(b, path) = backend_fill_curves(b, backend_stroke(b, path))

backend_stroke(b::Backend, path::Union{OpenPathSequence,ClosedPathSequence}) =
    backend_stroke_unite(b, map(path->backend_stroke(b, path), path.paths))
backend_fill(b::Backend, path::ClosedPathSequence) =
    backend_fill_curves(b, map(path->backend_stroke(b, path), path.paths))


#####################################################################
## Conversions

convert(::Type{ClosedPath}, s::Rectangle) =
  closed_path([MoveToOp(s.corner), RectOp(vxy(s.dx, s.dy))])
convert(::Type{ClosedPath}, s::Circle) =
  closed_path([CircleOp(s.center, s.radius)])
convert(::Type{Path}, s::Line) = convert(OpenPath, s.vertices)

convert(::Type{ClosedPath}, s::Polygon) =
    let vs = polygon_vertices(s)
        closed_path(vcat([MoveToOp(vs[1])], map(LineToOp, vs[2:end]), [CloseOp()]))
    end
#####################################################################

export curve_domain, surface_domain, frame_at
surface_domain(s::SurfaceRectangle) = (0, s.dx, 0, s.dy)
surface_domain(s::SurfaceCircle) = (0, s.radius, 0, 2pi)
surface_domain(s::SurfaceArc) = (0, s.radius, s.start_angle, s.amplitude)


frame_at(c::Shape1D, t::Real) = backend_frame_at(backend(c), c, t)
frame_at(s::Shape2D, u::Real, v::Real) = backend_frame_at(backend(s), s, u, v)

#Some specific cases can be handled in an uniform way without the backend
frame_at(s::SurfaceRectangle, u::Real, v::Real) = add_xy(s.corner, u, v)
frame_at(s::SurfaceCircle, u::Real, v::Real) = add_pol(s.center, u, v)




#####################################################################
# BIM
#=
Building Information Modeling is more than just shapes.
Each BIM element is connected to other BIM elements. A wall is connected to its
windows and doors. A floor is connected to its walls, and stairs. Stairs connect
different floors, etc.

This graph of objects needs to be build programmatically and that is one problem.
The other problem is portability, as we want this to work in backends that are
not BIM tools.

We will start with the simplest of BIM elements, namely, levels, slabs, walls,
windows and doors.
=#

abstract type BIMElement <: Proxy end
abstract type Measure <: BIMElement end
BIMElements = Vector{<:BIMElement}

@defproxy(level, Measure, height::Real=0, elements::BIMElements=BIMElement[])
levels_cache = Dict{Real,Level}()
maybe_replace(level::Level) = get!(levels_cache, level.height, level)

current_levels() = values(level_cache)
default_level = Parameter{Level}(level())
default_level_to_level_height = Parameter{Real}(3)
upper_level(lvl=default_level(), height=default_level_to_level_height()) = level(lvl.height + height, backend=backend(lvl))
Base.:(==)(l1::Level, l2::Level) = l1.height == l2.height

#default implementation
realize(b::Backend, s::Level) = s.height

export all_levels, default_level, default_level_to_level_height, upper_level

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
a parameter that contains a set of wall parameters.  As an example of use, we
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
    Dict(:width => "b", :height => "d"))

To make things more interesting, some families might require instantiation on
different levels. For example, a circular window family in Revit needs to be
loaded using RVTLoadFamily, which requires the name of the family, then it needs
to be instantiated with RVTFamilyElement, which requires the radius of the window,
and finally needs to be inserted on the wall, using RVTInsertWindow, which requires
the opening angle. This might be different for different families, so we need a
flexible way of using the parameters. One hipotesis is to specify those different
moments as different dictionaries.

RevitFamily(
  "C:\\ProgramData\\Autodesk\\RVT 2019\\Libraries\\US Metric\\Windows\\CIRCULAR WINDOW.rfa",
  ("radius_window" => :radius),
  ("angle_window" => :opening_angle))

If no parameters are needed on a particular phases, we might use an empty dictionary
to describe that phase. Given the typical Revit families, defaults parameters seem
to be the best approach here. This means that the previous beam family might be
equivalent to

RevitFamily(
  "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
  ("b" => :width, "d" => :height),
  ())

However, the same beam might have a different mapping in a different backend.
This means that we need another mapping to support different backends. One
possibility is to use something similar to:

backend_family(
    revit => RevitFamily(
        "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
        ("b" => :width, "d" => :height)),
    archicad => ArchiCADFamily(
        "BeamElement",
        ("size_x" => :width, "size_y" => :height)),
    autocad => AutoCADFamily())

Then, we need an operation that instantiates a family. This can be done on two different
levels: (1) from a backend-specific family (e.g., RevitFamily), for example:

beam_family = RevitFamily(
    "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
    ("b" => :width, "d" => :height))

current_beam_defaults(beam_family_instance(beam_family, width=10, height=20)

or from a generic backend family, for example:

beam_family = backend_family(
    revit => RevitFamily(
        "C:\\ProgramData\\Autodesk\\RVT 2017\\Libraries\\US Metric\\Structural Framing\\Wood\\M_Timber.rfa",
        ("b" => :width, "d" => :height)),
    archicad => ArchiCADFamily(
        "BeamElement"
        ("size_x" => :width, "size_y" => :height)),
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

#HACK Using Dict instead of IdDict just because get is not fully implemented for IdDict
#FIXME after updating to Julia 1.2
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
      based_on::Union{Family, Nothing}
      implemented_as::IdDict{<:Backend, <:Family}
      ref::Parameter{Any}
    end
    $(constructor_name)($(opt_params...);
                        $(key_params...),
                        based_on=nothing,
                        implemented_as=IdDict{Backend, Family}()) =
      $(struct_name)($(field_names...), based_on, implemented_as, Parameter{Any}(nothing))
    $(instance_name)(family:: Family, implemented_as=copy(family.implemented_as); $(instance_params...)) =
      $(struct_name)($(field_names...), family, implemented_as, Parameter{Any}(nothing))
    $(default_name) = Parameter($(constructor_name)())
    $(predicate_name)(v::$(struct_name)) = true
    $(predicate_name)(v::Any) = false
#    $(map((selector_name, field_name) -> :($(selector_name)(v::$(struct_name)) = v.$(field_name)),
#          selector_names, field_names)...)
    Khepri.meta_program(v::$(struct_name)) =
        Expr(:call, $(Expr(:quote, name)), $(map(field_name -> :(meta_program(v.$(field_name))), field_names)...))
  end
end

# When dispatching a BIM operation to a backend, we also need to dispatch the family

backend_family(b::Backend, family::Family) =
# replace this with next fragment after updating Julia
  if haskey(family.implemented_as, b)
    family.implemented_as[b]
  else
    family.based_on == nothing ? # this is not a family_element (nor a derivation of a family_element)
      error("Family $(family) is missing the implementation for backend $(b)") :
      backend_family(b, family.based_on)
  end
#=  get(family.implemented_as, b) do
    family.based_on == nothing ? # this is not a family_element (nor a derivation of a family_element)
      error("Family $(family) is missing the implementation for backend $(b)") :
      backend_family(b, family.based_on)
  end
=#

copy_struct(s::T) where T = T([getfield(s, k) for k ∈ fieldnames(T)]...)

# Backends will install their own families on top of the default families, e.g.,
# set_backend_family(default_beam_family(), revit, revit_beam_family)
set_backend_family(family::Family, backend::Backend, backend_family::Family) =
    family.implemented_as[backend]=backend_family

# Finally, we can implement a generic backend caching mechanism for families

realize(b::Backend, f::Family) =
  if f.ref() == nothing
    f.ref(backend_get_family_ref(b, f, backend_family(b, f)))
  else
    f.ref()
  end

export backend_family, set_backend_family

@deffamily(slab_family, Family,
    thickness::Real=0.2,
    coating_thickness::Real=0.0)

slab_family_elevation(b::Backend, family::SlabFamily) =
  family.coating_thickness - family.thickness
slab_family_thickness(b::Backend, family::SlabFamily) =
  family.coating_thickness + family.thickness

@defproxy(slab, Shape3D, contour::ClosedPath=rectangular_path(),
          level::Level=default_level(), family::SlabFamily=default_slab_family(),
          openings::Vector{<:ClosedPath}=ClosedPath[])

# Default implementation: dispatch on the slab elements
realize(b::Backend, s::Slab) =
    realize_slab(b, s.contour, s.openings, s.level, s.family)

realize_slab(b::Backend, contour::ClosedPath, holes::Vector{<:ClosedPath}, level::Level, family::Family) =
    let base = vz(level.height + slab_family_elevation(b, family)),
        thickness = slab_family_thickness(b, family)
        # Change this to a better named protocol?
        backend_slab(b, translate(contour, base), map(c -> translate(c, base), holes), thickness, family)
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
        s_thickness = slab_family_thickness(b, s.family)
        for opening in openings
            op_path = translate(opening, vz(s_base_height-1.1*s_thickness))
            op_ref = ensure_ref(b, backend_slab(b, op_path, s_thickness*1.2))
            s_ref = ensure_ref(b, subtract_ref(b, s_ref, op_ref))
        end
        s_ref
    end

# Roof

@deffamily(roof_family, Family,
    thickness::Real=0.2,
    coating_thickness::Real=0.0)

slab_family_elevation(b::Backend, family::RoofFamily) = 0
slab_family_thickness(b::Backend, family::RoofFamily) =
  family.coating_thickness + family.thickness

@defproxy(roof, Shape3D, contour::ClosedPath=rectangular_path(),
          level::Level=default_level(), family::RoofFamily=default_roof_family(),
          openings::Vector{<:ClosedPath}=ClosedPath[])

realize(b::Backend, s::Roof) =
    realize_slab(b, s.contour, s.openings, s.level, s.family)

# Panel

@deffamily(panel_family, Family,
    thickness::Real=0.02)

@defproxy(panel, Shape3D, vertices::Locs=Loc[], level::Any=default_level(), family::Any=default_panel_family())

#TODO Pass the provided backend
realize(b::Backend, s::Panel) =
  let #p1 = s.vertices[1],
      #p2 = s.vertices[2],
      #p3 = s.vertices[3],
      #n = vz(s.family.thickness, cs_from_o_vx_vy(p1, p2-p1, p3-p1))
      verts = in_world.(s.vertices)
      n = vertices_normal(verts)*(s.family.thickness/2)
    ref(irregular_prism(map(p -> p - n, verts), n*2))
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
          offset::Real=0.0,
          doors::Shapes=Shape[], windows::Shapes=Shape[])
wall(p0::Loc, p1::Loc;
     bottom_level::Level=default_level(),
     top_level::Level=upper_level(bottom_level),
     family::WallFamily=default_wall_family(),
     offset::Real=0.0) =
  wall([p0, p1],
       bottom_level=bottom_level,
       top_level=top_level,
       family=family,
       offset=offset)

#=
Walls can be joined. That is very important because the wall needs to have
uniform thickness along the entire path.
=#
export join_walls
join_walls(wall1, wall2) =
  if wall1.bottom_level != wall2.bottom_level
    error("Walls with different bottom levels")
  elseif wall1.top_level != wall2.top_level
    error("Walls with different top levels")
  elseif wall1.family != wall2.family
    error("Walls with different families")
  elseif wall1.offset != wall2.offset
    error("Walls with different offsets")
  else
    let w = wall(join_paths(wall1.path, wall2.path),
                 wall1.bottom_level, wall1.top_level,
                 wall1.family, wall1.offset),
        len = path_length(wall1.path)
      for (es,l) in ((wall1.doors, 0), (wall2.doors, len))
        for e in es
          add_door(w, e.loc+vx(l), e.family)
        end
      end
      for (es,l) in ((wall1.windows, 0), (wall2.windows, len))
        for e in es
          add_window(w, e.loc+vx(l), e.family)
        end
      end
      for w in (wall1, wall2)
        delete_shapes(w.doors)
        delete_shapes(w.windows)
        delete_shape(w)
      end
      w
    end
  end

join_walls(walls...) =
  reduce(join_walls, walls)



# Right and Left considering observer looking along with curve direction
r_thickness(w::Wall) = (+1+w.offset)/2*w.family.thickness
l_thickness(w::Wall) = (-1+w.offset)/2*w.family.thickness
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
      w_path = translate(w.path, vz(w_base_height)),
      r_thickness = r_thickness(w),
      l_thickness = l_thickness(w)
    ensure_ref(b, backend_wall(b, w_path, w_height, l_thickness, r_thickness, w.family))
  end

realize_wall_openings(b::Backend, w::Wall, w_ref, openings) =
  let w_base_height = w.bottom_level.height,
      w_height = w.top_level.height - w_base_height,
      w_path = translate(w.path, vz(w_base_height)),
      r_thickness = r_thickness(w),
      l_thickness = l_thickness(w)
    for opening in openings
      w_ref = realize_wall_opening(b, w_ref, w_path, l_thickness, r_thickness, opening, w.family)
      ref(opening)
    end
    w_ref
  end

realize_wall_opening(b::Backend, w_ref, w_path, l_thickness, r_thickness, op, family) =
  let op_base_height = op.loc.y,
      op_height = op.family.height,
      op_path = translate(subpath(w_path, op.loc.x, op.loc.x + op.family.width), vz(op_base_height)),
      op_ref = ensure_ref(b, backend_wall(b, op_path, op_height, l_thickness, r_thickness, family))
    ensure_ref(b, subtract_ref(b, w_ref, op_ref))
  end

realize(b::Backend, s::Union{Door, Window}) =
  let base_height = s.wall.bottom_level.height + s.loc.y,
      height = s.family.height,
      subpath = translate(subpath(s.wall.path, s.loc.x, s.loc.x + s.family.width), vz(base_height))
      backend_wall_element(b, s, subpath, height, s.family.thickness, s.family)
  end

backend_wall_element(b::Backend, s::Union{Door, Window}, path, height, thickness, family) =
  # we emulate doors and windows using a small wall
  backend_wall(b::Backend, path, height, thickness/2, thickness/2, family)

##

export add_door
add_door(w::Wall=required(), loc::Loc=u0(), family::DoorFamily=default_door_family()) =
  backend_add_door(backend(w), w, loc, family)

backend_add_door(b::Backend, w::Wall, loc::Loc, family::DoorFamily) =
    let d = door(w, loc, family=family)
        push!(w.doors, d)
        if realized(w)
            set_ref!(w, realize_wall_openings(b, w, ref(w), [d]))
        end
        w
    end

#
export add_window
add_window(w::Wall=required(), loc::Loc=u0(), family::WindowFamily=default_window_family()) =
  backend_add_window(backend(w), w, loc, family)

backend_add_window(b::Backend, w::Wall, loc::Loc, family::WindowFamily) =
  let d = window(w, loc, family=family)
    push!(w.windows, d)
    if realized(w)
        set_ref!(w, realize_wall_openings(b, w, ref(w), [d]))
    end
    w
  end

#=
A curtain wall is a special kind of wall that is made of a frame with windows.
=#

@deffamily(curtain_wall_frame_family, Family,
  width::Real=0.1,
  depth::Real=0.1,
  depth_offset::Real=0.25)

@deffamily(curtain_wall_family, Family,
  n_curtain_panels::Int=3,
  panel::PanelFamily=panel_family(thickness=0.05),
  boundary_frame::CurtainWallFrameFamily=
    curtain_wall_frame_family(width=0.1,depth=0.1,depth_offset=0.25),
  mullion_frame::CurtainWallFrameFamily=
    curtain_wall_frame_family(width=0.08,depth=0.09,depth_offset=0.2),
  transom_frame::CurtainWallFrameFamily=
    curtain_wall_frame_family(width=0.06,depth=0.1,depth_offset=0.11))

@defproxy(curtain_wall, Shape3D,
          path::Path=rectangular_path(),
          bottom_level::Level=default_level(),
          top_level::Level=upper_level(bottom_level),
          family::CurtainWallFamily=default_curtain_wall_family(),
          offset::Real=0.0)
curtain_wall(p0::Loc, p1::Loc;
     bottom_level::Level=default_level(),
     top_level::Level=upper_level(bottom_level),
     family::CurtainWallFamily=default_curtain_wall_family(),
     offset::Real=0.0) =
  curtain_wall([p0, p1], bottom_level=bottom_level, top_level=top_level,
         family=family, offset=offset)

realize(b::Backend, s::CurtainWall) =
  let th = s.family.panel.thickness,
      bfw = s.family.boundary_frame.width,
      bfd = s.family.boundary_frame.depth,
      bfdo = s.family.boundary_frame.depth_offset,
      mfw = s.family.mullion_frame.width,
      mfd = s.family.mullion_frame.depth,
      mdfo = s.family.mullion_frame.depth_offset,
      tfw = s.family.transom_frame.width,
      tfd = s.family.transom_frame.depth,
      tfdo = s.family.transom_frame.depth_offset,
      #pts = map(t->in_world(location_at_length(s.path, t)),
      #          division(0, path_length, s.family.n_curtain_panels)),
      path = s.path, #open_polygonal_path(pts),
      path_length = path_length(path),
      bottom = level_height(s.bottom_level),
      top = level_height(s.top_level),
      height = top - bottom,
      refs = []
    push!(refs, backend_curtain_wall(b, s, subpath(path, bfw, path_length-bfw), bottom+bfw, height-2*bfw, th, :panel))
    push!(refs, backend_curtain_wall(b, s, path, bottom, bfw, bfd, :boundary_frame))
    push!(refs, backend_curtain_wall(b, s, path, top-bfw, bfw, bfd, :boundary_frame))
    push!(refs, backend_curtain_wall(b, s, subpath(path, 0, bfw), bottom+bfw, height-2*bfw, bfd, :boundary_frame))
    push!(refs, backend_curtain_wall(b, s, subpath(path, path_length-bfw, path_length), bottom+bfw, height-2*bfw, bfd, :boundary_frame))
    push!(refs, backend_curtain_wall(b, s, subpath(path, bfw, path_length-bfw), bottom+height/2-tfw/2, tfw, tfd, :transom_frame))
    let n = s.family.n_curtain_panels
      for i in 1:n-1
        l = path_length/n*i
        push!(refs, backend_curtain_wall(b, s, subpath(path, l-mfw/2, l+mfw/2), bottom+bfw, height-2*bfw, mfd, :mullion_frame))
      end
    end
    [ensure_ref(b,r) for r in refs]
  end

backend_curtain_wall(b::Backend, s, path::Path, bottom::Real, height::Real, thickness::Real, kind::Symbol) =
  backend_wall(b, translate(path, vz(bottom)), height, thickness/2, thickness/2, getproperty(s.family, kind))

#
# We need to redefine the default method (maybe add an option to the macro to avoid defining the meta_program)
# This needs to be fixed for windows
#=
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
=#

# Beam
# Beams are mainly horizontal elements. By default, a beam is aligned along its top axis
@deffamily(beam_family, Family,
#  width::Real=1.0,
#  height::Real=2.0,
  profile::ClosedPath=top_aligned_rectangular_profile(1, 2))
#beam_family(Width::Real=1.0, Height::Real=2.0; width=Width, height=Height) =
#  beam_family(rectangular_path(xy(-width/2,-height), width, height))

@defproxy(beam, Shape3D, cb::Loc=u0(), h::Real=1, angle::Real=0, family::BeamFamily=default_beam_family())
beam(cb::Loc, ct::Loc, Angle::Real=0, Family::BeamFamily=default_beam_family(); angle::Real=Angle, family::BeamFamily=Family) =
    let (c, h) = position_and_height(cb, ct)
      beam(c, h, angle, family)
    end

# Column
# Columns are mainly vertical elements. A column has its center axis aligned with a line defined by two points

@deffamily(column_family, Family,
    #width::Real=1.0,
    #height::Real=2.0,
  profile::ClosedPath=rectangular_profile(1, 2))

@defproxy(free_column, Shape3D, cb::Loc=u0(), h::Real=1, angle::Real=0, family::ColumnFamily=default_column_family())
free_column(cb::Loc, ct::Loc, Angle::Real=0, Family::ColumnFamily=default_column_family(); angle::Real=Angle, family::ColumnFamily=Family) =
    let (c, h) = position_and_height(cb, ct)
      free_column(c, h, angle, family)
    end

@defproxy(column, Shape3D, cb::Loc=u0(), angle::Real=0,
  bottom_level::Level=default_level(), top_level::Level=upper_level(bottom_level),
  family::ColumnFamily=default_column_family())

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

@defproxy(pointlight, Shape3D, loc::Loc=z(3), color::RGB=rgb(255,255,255), range::Real=10, intensity::Real=4, level::Level=default_level())

realize(b::Backend, s::Pointlight) =
    backend_pointlight(b, add_z(s.loc, s.level.height), s.color, s.range, s.intensity)

@defproxy(spotlight, Shape3D, loc::Loc=z(3), dir::Vec=vz(-1), hotspot::Real=pi/4, falloff::Real=pi/3)

realize(b::Backend, s::Spotlight) =
    backend_spotlight(b, s.loc, s.dir, s.hotspot, s.falloff)

@defproxy(ieslight, Shape3D, file::String=required(), loc::Loc=z(3), dir::Vec=vz(-1), alpha::Real=0, beta::Real=0, gamma::Real=0)

realize(b::Backend, s::Ieslight) =
    backend_ieslight(b, s.file, s.loc, s.dir, s.alpha, s.beta, s.gamma)


#################################

@deffamily(truss_node_family, Family,
    radius::Real=0.2,
    support::Any=false) #(Option node_support)

@deffamily(truss_bar_family, Family,
    radius::Real=0.03,
    section::Any=false,
    material::Any=false,
    created::Parameter{Bool}=Parameter(false)) # HACK: This should be merged with the lazy creation of families

@defproxy(truss_node, Shape3D, p::Loc=u0(), family::TrussNodeFamily=default_truss_node_family())
@defproxy(truss_bar, Shape3D, p0::Loc=u0(), p1::Loc=u0(), angle::Real=0, family::TrussBarFamily=default_truss_bar_family())

realize(b::Backend, s::TrussNode) = sphere(s.p, s.family.radius)
realize(b::Backend, s::TrussBar) = cylinder(s.p0, s.family.radius, s.p1)

import Base.union
export union, intersection, subtraction

@defproxy(union_shape, Shape3D, shapes::Shapes=Shape[])
union(shapes::Shapes) =
  let non_empty_shapes = filter(s -> !is_empty_shape(s), shapes)
    non_empty_shapes == [] ?
      empty_shape() :
      union_shape(non_empty_shapes)
  end

union(shape::Shape, shapes...) = union([shape, shapes...])

@defproxy(intersection_shape, Shape3D, shapes::Shapes=Shape[])
intersection(shapes::Shapes) = intersection_shape(shapes)
intersection(shape::Shape, shapes...) = intersection_shape([shape, shapes...])

@defproxy(subtraction_shape2D, Shape2D, shape::Shape=surface_circle(), shapes::Shapes=Shape[])
@defproxy(subtraction_shape3D, Shape3D, shape::Shape=surface_sphere(), shapes::Shapes=Shape[])
subtraction(shape::Shape2D, shapes...) = subtraction_shape2D(shape, [shapes...])
subtraction(shape::Shape3D, shapes...) = subtraction_shape3D(shape, [shapes...])

@defproxy(slice, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))

@defproxy(mirror, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))
@defproxy(union_mirror, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))

@defproxy(surface_grid, Shape2D, points::AbstractMatrix{<:Loc}=zeros(Loc,(2,2)), closed_u::Bool=false, closed_v::Bool=false,
          interpolator::Parameter{Any}=Parameter{Any}(missing))
surface_interpolator(pts::AbstractMatrix{<:Loc}) =
    let pts = map(pts) do p
                let v = in_world(p).raw
                  SVector{3,Float64}(v[1], v[2], v[3])
                end
              end
        Interpolations.scale(
            interpolate(pts, BSpline(Cubic(Natural(OnGrid())))),
            range(0,stop=1,length=size(pts, 1)),
            range(0,stop=1,length=size(pts, 2)))
    end

# For interpolator to work, we need this:

convert(::Type{AbstractMatrix{<:Loc}}, pts::Vector{<:Vector{<:Loc}}) =
  permutedims(hcat(pts...))

evaluate(s::SurfaceGrid, u::Real, v::Real) =
  let interpolator = s.interpolator
    if ismissing(interpolator())
      interpolator(surface_interpolator(s.points))
    end
    let p = interpolator()(u,v)
        v = Interpolations.gradient(interpolator(), u, v)
      loc_from_o_vx_vy(
        xyz(p[1], p[2], p[3], world_cs),
        vxyz(v[1][1], v[1][2], v[1][3], world_cs),
        vxyz(v[2][1], v[2][2], v[2][3], world_cs))
    end
  end

surface_domain(s::SurfaceGrid) = (0.0, 1.0, 0.0, 1.0)
frame_at(s::SurfaceGrid, u::Real, v::Real) = evaluate(s, u, v)
map_division(f::Function, s::SurfaceGrid, nu::Int, nv::Int, backend::Backend=current_backend()) =
  let (u1, u2, v1, v2) = surface_domain(s)
    map_division(u1, u2, nu) do u
      map_division(v1, v2, nv) do v
        f(frame_at(s, u, v))
      end
    end
  end


@defproxy(thicken, Shape3D, shape::Shape=surface_circle(), thickness::Real=1)

# Blocks

@defproxy(block, Shape, name::String="Block", shapes::Shapes = Shape[])
@defproxy(block_instance, Shape, block::Block=required(), loc::Loc=u0(), scale::Real=1.0)

################################################################################

#Backends might use different communication mechanisms, e.g., sockets, COM, RMI, etc

#We start with socket-based communication
struct SocketBackend{K,T} <: Backend{K,T}
  connection::LazyParameter{TCPSocket}
  remote::NamedTuple
end

# To simplify remote calls
macro remote(b, call)
  let op = call.args[1],
      args = map(esc, call.args[2:end]),
      b = esc(b)
    :(call_remote(getfield(getfield($(b), :remote), $(QuoteNode(op))), $(b).connection(), $(args...)))
  end
end

macro get_remote(b, op)
  let b = esc(b)
    :(getfield(getfield($(b), :remote), $(QuoteNode(op))))
  end
end


SocketBackend{K,T}(c::LazyParameter{TCPSocket}) where {K,T} =
  SocketBackend{K,T}(c, NamedTuple{}())

connection(b::SocketBackend{K,T}) where {K,T} = b.connection()

reset_backend(b::SocketBackend) =
  begin
    for f in b.remote
      reset_opcode(f)
    end
    close(b.connection())
    reset(b.connection)
  end

#One less dynamic option is to use a file-based backend. To that end, we implement
#the IOBuffer_Backend

struct IOBufferBackend{K,T} <: Backend{K,T}
  out::IOBuffer
end
connection(backend::IOBufferBackend) = backend.out


################################################################################
bounding_box(shape::Shape) =
  bounding_box([shape])

bounding_box(shapes::Shapes=Shape[]) =
  if isempty(shapes)
    [u0(), u0()]
  else
    backend_bounding_box(backend(shapes[1]), shapes)
  end

delete_shape(shape::Shape) =
  delete_shapes([shape])

delete_shapes(shapes::Shapes=Shape[]) =
  if ! isempty(shapes)
    to_delete = filter(realized, shapes)
    backend_delete_shapes(backend(shapes[1]), to_delete)
    foreach(mark_deleted, shapes)
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

and_mark_deleted(r::Any, shape) =
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
@defop delete_all_shapes_in_layer(layer)

# BIM
@defop all_levels()
@defop all_walls()
@defop all_walls_at_level(level)

@defop disable_update()
@defop enable_update()
@defop set_view(camera::Loc, target::Loc, lens::Real)
@defop get_view()
@defop zoom_extents()
@defop view_top()
@defop get_layer(name::String)
@defop create_layer(name::String)
@defop current_layer()
@defop current_layer(layer)
@defop set_layer_active(layer, status)
@defop switch_to_layer(layer)
@defop get_material(name::String)
@defop create_material(name::String)
@defop current_material()
@defop current_material(material)

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

@defop select_position(prompt::String="Select a position")
@defop select_positions(prompt::String="Select positions")
@defop select_point(prompt::String="Select a point")
@defop select_points(prompt::String="Select points")
@defop select_curve(prompt::String="Select a curve")
@defop select_curves(prompt::String="Select curves")
@defop select_surface(prompt::String="Select a surface")
@defop select_surfaces(prompt::String="Select surfaces")
@defop select_solid(prompt::String="Select a solid")
@defop select_solids(prompt::String="Select solids")
@defop select_shape(prompt::String="Select a shape")
@defop select_shapes(prompt::String="Select shapes")
@defop highlight_shapes(shapes::Shapes)
@defshapeop highlight_shape()
@defshapeop register_for_changes()
@defshapeop unregister_for_changes()
@defshapeop waiting_for_changes()
@defop changed_shape(shapes::Shapes)

capture_shape(s=select_shape("Select shape to be captured")) =
  if s != nothing
    generate_captured_shape(s, backend(s))
  end

capture_shapes(ss=select_shapes("Select shapes to be captured")) =
    generate_captured_shapes(ss, backend(ss[1]))

register_for_changes(shapes::Shapes) =
  map(shapes) do shape
    register_for_changes(shape, backend(shape))
  end

unregister_for_changes(shapes::Shapes) =
  map(shapes) do shape
    unregister_for_changes(shape, backend(shape))
  end

waiting_for_changes(shapes::Shapes) =
  waiting_for_changes(shapes[1], backend(shapes[1]))

export on_change
on_change(f, shape::Shape) = on_change(f, [shape])
on_change(f, shapes) =
  let registered = register_for_changes(shapes)
    try
      while waiting_for_changes(shapes)
        let changed = changed_shape(shapes)
          f()
        end
      end
    finally
      unregister_for_changes(registered)
    end
  end

#
export with_shape_dependency
with_shape_dependency(f, ss) =
    let shapes = collecting_shapes() do
                    f()
                end
        on_change(ss) do
            try
                delete_shapes(shapes)
            catch e
            end
            shapes = collecting_shapes() do
                f()
            end
        end
    end

#
export internalize_shape, internalize_shapes

internalize_shape(s=select_shape("Select shape to be internalized")) =
  if s != nothing
    println(meta_program(s))
  end

internalize_shapes(ss=select_shapes("Select shapes to be internalized")) =
    println(meta_program(ss))

# Later, this will be used to create images.
export to_render
to_render(f, name) =
  begin
    delete_all_shapes()
    f()
  end

# Seletion

select_one_with_prompt(prompt::String, b::Backend, f::Function) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = f(connection(b), prompt)
      length(ans) > 0 ? shape_from_ref(ans[1], b) : nothing
    end
  end

select_many_with_prompt(prompt::String, b::Backend, f::Function) =
  begin
    @info "$(prompt) on the $(b) backend."
    map(id -> shape_from_ref(id, b), f(connection(b), prompt))
  end

export render_view
render_view(name::String) =
  let path = prepare_for_saving_file(render_pathname(name))
    render_view(path, current_backend())
    path
  end
