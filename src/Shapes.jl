export Shape,
       Shapes,
       Path,
       backend,
       new_backend,
       backend_name,
       current_backend,
       has_current_backend,
       switch_to_backend,
       void_ref,
       delete_shape, delete_shapes,
       delete_all_shapes,
       set_length_unit,
       is_collecting_shape,
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

show(io::IO, b::Backend{K,R}) where {K,R} = print(io, backend_name(b))

backend_name(b::Backend{K,R}) where {K,R} = typeof(b)

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
# currying
unite_ref(b::Backend{K,T}) where {K,T} = (r0::GenericRef{K,T}, r1::GenericRef{K,T}) -> unite_ref(b, r0, r1)

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

# currying
intersect_ref(b::Backend{K,T}) where {K,T} = (r0::GenericRef{K,T}, r1::GenericRef{K,T}) -> intersect_ref(b, r0, r1)

intersect_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UniversalRef{K,T}) where {K,T} = r0
intersect_ref(b::Backend{K,T}, r0::UniversalRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r1
intersect_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::EmptyRef{K,T}) where {K,T} = r1
intersect_ref(b::Backend{K,T}, r0::EmptyRef{K,T}, r1::GenericRef{K,T}) where {K,T} = r0
intersect_ref(b::Backend{K,T}, r0::GenericRef{K,T}, r1::UnionRef{K,T}) where {K,T} =
  intersect_ref(b, r0, unite_refs(b, r1))
intersect_ref(b::Backend{K,T}, r0::UnionRef{K,T}, r1::GenericRef{K,T}) where {K,T} =
  intersect_ref(b, unite_refs(b, r0), r1)

#To avoid ambiguity
# currying
subtract_ref(b::Backend{K,T}) where {K,T} = (r0::GenericRef{K,T}, r1::GenericRef{K,T}) -> subtract_ref(b, r0, r1)

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
show(io::IO, s::Shape) =
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
delay_realize(b::LazyBackend, s::Shape) = save_shape!(b, s)

# By default, save_shape! assumes there is a field in the backend to store shapes
save_shape!(b::Backend, s::Shape) = (push!(b.shapes, s); s)


with_transaction(fn) =
  maybe_realize(with(fn, delaying_realize, true))

#=
Frequently, we need to collect all shapes that are created:
=#

# HACK: Replace in_shape_collection with is_collecting_shapes
in_shape_collection = Parameter(false)
is_collecting_shapes = in_shape_collection
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
# The undefined backend
struct UndefinedBackend <: Backend{Int,Int} end
connection(b::UndefinedBackend) = throw(UndefinedBackendException())
void_ref(b::UndefinedBackend) = EmptyRef{Int,Int}()
const undefined_backend = UndefinedBackend()

const current_backend = Parameter{Backend}(undefined_backend)
has_current_backend() = current_backend() != undefined_backend

# Side-effect full operations need to have a backend selected and will generate an exception if there is none

struct UndefinedBackendException <: Exception end
show(io::IO, e::UndefinedBackendException) = print(io, "No current backend.")

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

macro defopnamed(name_params)
    name, params = name_params.args[1], name_params.args[2:end]
    par_names = map(par -> par.args[1].args[1], params)
    par_types = map(par -> par.args[1].args[2], params)
    par_inits = map(par -> par.args[2], params)
    backend_call = Symbol("backend_", name)
    quote
        export $(name)
        $(process_named_params(:($(name)($(params...), backend::Backend=current_backend()) =
            $(backend_call)(backend, $(par_names...)))))
        $(esc(:($(backend_call)(backend::Backend, $(map((name, typ)->Expr(:(::), name, typ), par_names, par_types)...)) =
            UndefinedBackendException())))
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

new_backend(b::Backend = current_backend()) = backend(b)

struct WrongTypeForParam <: Exception
  param::Symbol
  value::Any
  expected_type::Type
end
Base.showerror(io::IO, e::WrongTypeForParam) =
  print(io, "$(e.param) expected a $(e.expected_type) but got $(e.value) of type $(typeof(e.value))")

macro defproxy(name_typename, parent, fields...)
  (name, typename) = name_typename isa Symbol ?
    (name_typename, Symbol(string(map(uppercasefirst,split(string(name_typename),'_'))...))) :
    name_typename.args
  name_str = string(name)
  struct_name = esc(typename)
  field_names = map(field -> field.args[1].args[1], fields)
  field_types = map(field -> esc(field.args[1].args[2]), fields)
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
    # we don't need to convert anything because Julia already does that with the default constructor
    # and, by the same idea, we don't need to define parameter types.
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

# This might be usable, so
export @defproxy, realize, Shape0D, Shape1D, Shape2D, Shape3D

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
=#
#=
curve_domain(s::Spline) = (0.0, 1.0)
frame_at(s::Spline, t::Real, backend::Backend=backend(s)) = evaluate(s, t)
=#
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
realize(b::Backend, s::Polygon) =
  backend_polygon(b, s.vertices)
@defproxy(regular_polygon, Shape1D, edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, inscribed::Bool=true)
realize(b::Backend, s::RegularPolygon) =
  backend_polygon(b, regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
@defproxy(rectangle, Shape1D, corner::Loc=u0(), dx::Real=1, dy::Real=1)
rectangle(p::Loc, q::Loc) =
  let v = in_cs(q - p, p.cs)
    rectangle(p, v.x, v.y)
  end
realize(b::Backend, s::Rectangle) =
  backend_polygon(b, [
    s.corner,
    add_x(s.corner, s.dx),
    add_xy(s.corner, s.dx, s.dy),
    add_y(s.corner, s.dy)])
@defproxy(surface_circle, Shape2D, center::Loc=u0(), radius::Real=1)
@defproxy(surface_arc, Shape2D, center::Loc=u0(), radius::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(surface_elliptic_arc, Shape2D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1, start_angle::Real=0, amplitude::Real=pi)
@defproxy(surface_ellipse, Shape2D, center::Loc=u0(), radius_x::Real=1, radius_y::Real=1)
@defproxy(surface_polygon, Shape2D, vertices::Locs=[u0(), ux(), uy()])
surface_polygon(v0, v1, vs...) = surface_polygon([v0, v1, vs...])
realize(b::Backend, s::SurfacePolygon) =
  backend_surface_polygon(b, s.vertices)
@defproxy(surface_regular_polygon, Shape2D, edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, inscribed::Bool=true)
realize(b::Backend, s::SurfaceRegularPolygon) =
  backend_surface_polygon(b, regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
@defproxy(surface_rectangle, Shape2D, corner::Loc=u0(), dx::Real=1, dy::Real=1)
surface_rectangle(p::Loc, q::Loc) =
  let v = in_cs(q - p, p.cs)
    surface_rectangle(p, v.x, v.y)
  end
realize(b::Backend, s::SurfaceRectangle) =
  let c = s.corner,
      dx = s.dx,
      dy = s.dy
    backend_surface_polygon(b, [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy), c])
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






path_vertices(s::Shape1D) = path_vertices(shape_path(s))
shape_path(s::Circle) = circular_path(s.center, s.radius)
shape_path(s::Spline) = open_spline_path(s.points, s.v0, s.v1)
shape_path(s::ClosedSpline) = closed_spline_path(s.points)




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
realize(b::Backend, s::Cuboid) =
  backend_pyramid_frustum(b, [s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3])

@defproxy(regular_pyramid_frustum, Shape3D, edges::Integer=4, cb::Loc=u0(), rb::Real=1, angle::Real=0, h::Real=1, rt::Real=1, inscribed::Bool=true)
regular_pyramid_frustum(edges::Integer, cb::Loc, rb::Real, angle::Real, ct::Loc, rt::Real=1, inscribed::Bool=true) =
  let (c, h) = position_and_height(cb, ct)
    regular_pyramid_frustum(edges, c, rb, angle, h, rt, inscribed)
  end
realize(b::Backend, s::RegularPyramidFrustum) =
  backend_pyramid_frustum(
    b,
    regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
    regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed))
backend_pyramid_frustum(b::Backend, bot_vs::Locs, top_vs::Locs) =
  let refs = [backend_surface_polygon(b, reverse(bot_vs))]
    push!(refs, backend_surface_polygon(b, top_vs))
    for (v1, v2, v3, v4) in zip(bot_vs, circshift(bot_vs, -1), circshift(top_vs, -1), top_vs)
      push!(refs, backend_surface_polygon(b, [v1, v2, v3, v4]))
    end
    refs
  end

@defproxy(regular_pyramid, Shape3D, edges::Integer=3, cb::Loc=u0(), rb::Real=1, angle::Real=0, h::Real=1, inscribed::Bool=true)
regular_pyramid(edges::Integer, cb::Loc, rb::Real, angle::Real, ct::Loc, inscribed::Bool=true) =
  let (c, h) = position_and_height(cb, ct)
    regular_pyramid(edges, c, rb, angle, h, inscribed)
  end
realize(b::Backend, s::RegularPyramid) =
  backend_pyramid(
    b,
    regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
    add_z(s.cb, s.h))

backend_pyramid(b::Backend, bot_vs::Locs, top::Loc) =
  let refs = [backend_surface_polygon(b, reverse(bot_vs))]
    for (v1, v2, v3) in zip(bot_vs, circshift(bot_vs, -1), repeated(top))
      push!(refs, backend_surface_polygon(b, [v1, v2, v3]))
    end
    refs
  end

@defproxy(irregular_pyramid_frustum, Shape3D, bs::Locs=[ux(), uy(), uxy()], ts::Locs=[uxz(), uyz(), uxyz()])
realize(b::Backend, s::IrregularPyramidFrustum) =
  backend_pyramid_frustum(b, s.bs, s.ts)

@defproxy(irregular_pyramid, Shape3D, bs::Locs=[ux(), uy(), uxy()], t::Loc=uz())
realize(b::Backend, s::IrregularPyramid) =
  backend_pyramid(b, s.bs, s.t)

@defproxy(regular_prism, Shape3D, edges::Integer=3, cb::Loc=u0(), r::Real=1, angle::Real=0, h::Real=1, inscribed::Bool=true)
regular_prism(edges::Integer, cb::Loc, r::Real, angle::Real, ct::Loc, inscribed::Bool=true) =
  let (c, h) = position_and_height(cb, ct)
    regular_prism(edges, c, r, angle, h, inscribed)
  end
realize(b::Backend, s::RegularPrism) =
  let ps = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
    backend_pyramid_frustum(b, ps, map(p -> add_z(p, s.h), ps))
  end

@defproxy(irregular_prism, Shape3D, bs::Locs=[ux(), uy(), uxy()], v::Vec=vz(1))
irregular_prism(bs::Locs, h::Real) =
  irregular_prism(bs, vz(h))
realize(b::Backend, s::IrregularPrism) =
  backend_pyramid_frustum(b, s.bs, map(p -> (p + s.v), s.bs))

@defproxy(right_cuboid, Shape3D, cb::Loc=u0(), width::Real=1, height::Real=1, h::Real=1)
right_cuboid(cb::Loc, width::Real, height::Real, ct::Loc, angle::Real=0; backend::Backend=current_backend()) =
  let (c, h) = position_and_height(cb, ct),
      o = angle == 0 ? c : loc_from_o_phi(c, angle)
    right_cuboid(o, width, height, h, backend=backend)
  end
realize(b::Backend, s::RightCuboid) =
  backend_right_cuboid(b, s.cb, s.width, s.height, s.h, nothing)

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
realize(b::Backend, s::Cylinder) =
  backend_cylinder(b, s.cb, s.r, s.h)

@defproxy(extrusion, Shape3D, profile::Shape=point(), v::Vec=vz(1))
extrusion(profile, h::Real) =
  extrusion(profile, vz(h))

realize(b::Backend, s::Extrusion) =
  backend_extrusion(backend(s), s.profile, s.v)

backend_extrusion(b::Backend, p::Point, v::Vec) =
  realize_and_delete_shapes(line([p.position, p.position + v], backend=b), [p])

@defproxy(sweep, Shape3D, path::Union{Shape1D, Path}=circle(), profile::Union{Shape,Path}=point(), rotation::Real=0, scale::Real=1)

realize(b::Backend, s::Sweep) =
  backend_sweep(backend(s), s.path, s.profile, s.rotation, s.scale)

backend_sweep(b::Backend, path::Union{Shape,Path}, profile::Union{Shape,Path}, rotation::Real, scale::Real) =
  let vertices = in_world.(path_vertices(profile)),
      frames = map_division(identity, path, 100) #rotation_minimizing_frames(path_frames(path))
    backend_surface_grid(
      b,
      [xyz(cx(p), cy(p), cz(p), frame.cs) for p in vertices, frame in frames],
      is_closed_path(profile),
      is_closed_path(path),
      is_smooth_path(profile),
      is_smooth_path(path))
  end

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

backend_stroke(b::Backend, path::RectangularPath) =
  let c = path.corner,
      dx = path.dx,
      dy = path.dy
    backend_stroke_line(b, (c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy), c))
  end

backend_stroke(b::Backend, path::OpenPolygonalPath) =
	backend_stroke_line(b, path.vertices)

backend_stroke(b::Backend, path::ClosedPolygonalPath) =
  backend_stroke_line(b, [path.vertices...,path.vertices[1]])

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

backend_stroke(b::Backend, m::Mesh) =
  let vs = m.vertices
    for face in m.faces
      backend_stroke_line(b, vs[face.+1]) #1-indexed
    end
  end
backend_fill(b::Backend, m::Mesh) =
  backend_surface_mesh(b, m.vertices, m.faces)

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
## Paths can be used to generate surfaces and solids

@defproxy(sweep_path, Shape3D, path::Path=polygonal_path(), profile::Path=circular_path(), rotation::Real=0, scale::Real=1)


#####################################################################
export curve_domain, surface_domain, frame_at, curve_length, frame_at_length
surface_domain(s::SurfaceRectangle) = (0, s.dx, 0, s.dy)
surface_domain(s::SurfaceCircle) = (0, s.radius, 0, 2pi)
surface_domain(s::SurfaceArc) = (0, s.radius, s.start_angle, s.amplitude)

curve_length(c::Shape1D) = backend_curve_length(backend(c), c)
frame_at(c::Shape1D, t::Real) = backend_frame_at(backend(c), c, t)
frame_at_length(c::Shape1D, t::Real) = backend_frame_at_length(backend(c), c, t)
frame_at(s::Shape2D, u::Real, v::Real) = backend_frame_at(backend(s), s, u, v)

#Some specific cases can be handled in an uniform way without the backend
frame_at(s::SurfaceRectangle, u::Real, v::Real) = add_xy(s.corner, u, v)
frame_at(s::SurfaceCircle, u::Real, v::Real) = add_pol(s.center, u, v)

export union, intersection, subtraction
#=
We do some pre-filtering to deal with the presence of empty shapes or to simplify one-arg cases.
=#

@defproxy(union_shape, Shape3D, shapes::Shapes=Shape[])
union(shapes::Shapes) =
  let non_empty_shapes = filter(s -> !is_empty_shape(s), shapes),
      count_non_empty_shapes = length(non_empty_shapes)
    count_non_empty_shapes == 0 ? empty_shape() :
    count_non_empty_shapes == 1 ? non_empty_shapes[1] :
    union_shape(non_empty_shapes)
  end

union(shape::Shape, shapes...) = union([shape, shapes...])

@defproxy(intersection_shape, Shape3D, shapes::Shapes=Shape[])
intersection(shapes::Shapes) = intersection_shape(shapes)
intersection(shape::Shape, shapes...) =
  is_empty_shape(shape) || any(is_empty_shape, shapes) ? empty_shape() :
  shapes == [] ? shape : intersection_shape([shape, shapes...])

@defproxy(subtraction_shape2D, Shape2D, shape::Shape=surface_circle(), shapes::Shapes=Shape[])
@defproxy(subtraction_shape3D, Shape3D, shape::Shape=surface_sphere(), shapes::Shapes=Shape[])
subtraction(shape::Shape2D, shapes...) =
  is_empty_shape(shape) ? empty_shape() :
    let non_empty_shapes = filter(s -> !is_empty_shape(s), shapes),
        count_non_empty_shapes = length(non_empty_shapes)
      count_non_empty_shapes == 0 ? shape : subtraction_shape2D(shape, [non_empty_shapes...])
    end
subtraction(shape::Shape3D, shapes...) =
  is_empty_shape(shape) ? empty_shape() :
    let non_empty_shapes = filter(s -> !is_empty_shape(s), shapes),
        count_non_empty_shapes = length(non_empty_shapes)
      count_non_empty_shapes == 0 ? shape : subtraction_shape3D(shape, [non_empty_shapes...])
    end

@defproxy(slice, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))

@defproxy(mirror, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))
@defproxy(union_mirror, Shape3D, shape::Shape=sphere(), p::Loc=u0(), n::Vec=vz(1))

@defproxy(surface_grid, Shape2D, points::AbstractMatrix{<:Loc}=zeros(Loc,(2,2)),
          closed_u::Bool=false, closed_v::Bool=false,
          smooth_u::Bool=true, smooth_v::Bool=true,
          interpolator::Parameter{Any}=Parameter{Any}(missing))

#=
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
=#

# For interpolator to work, we need this:

convert(::Type{AbstractMatrix{<:Loc}}, pts::Vector{<:Vector{<:Loc}}) =
  permutedims(hcat(pts...))

#=
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
=#

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

realize(b::Backend, s::SurfaceGrid) =
  backend_surface_grid(b, s.points, s.closed_u, s.closed_v, s.smooth_u, s.smooth_v)

@defproxy(surface_mesh, Shape2D, vertices::Locs=[u0(), ux(), uy()], faces::Vector{Vector{Int}}=[[0,1,2]])
realize(b::Backend, s::SurfaceMesh) =
  backend_surface_mesh(b, s.vertices, s.faces)



@defproxy(parametric_surface, Shape2D, definition::Function=(u,v)->xyz(u,v,0),
          domain_u::Tuple{Real,Real}=(0,1), domain_v::Tuple{Real,Real}=(0,1))

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

realize(b::Backend, s::Union{SubtractionShape2D,SubtractionShape3D}) =
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
@defop disable_update()
@defop enable_update()
@defop set_view(camera::Loc, target::Loc, lens::Real=50, aperture::Real=32)
@defop get_view()
@defop set_sun(altitude::Real, azimuth::Real)
@defop add_ground_plane()
@defop zoom_extents()
@defop view_top()
@defop get_layer(name::String)
@defopnamed create_layer(name::String="Layer", active::Bool=true, color::RGB=rgb(1,1,1))
@defop current_layer()
@defop current_layer(layer)
@defop set_layer_active(layer, status)
@defop switch_to_layer(layer)
@defop get_material(name::String)
@defop create_material(name::String)
@defop current_material()
@defop current_material(material)
@defop set_normal_sky()
@defop set_overcast_sky()

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
@defshapeop highlight_shape()
@defshapeop register_for_changes()
@defshapeop unregister_for_changes()
@defshapeop waiting_for_changes()
@defop changed_shape(shapes::Shapes)

export highlight_shapes
highlight_shapes(shapes::Shapes) =
  highlight_shapes(shapes, shapes == [] ? current_backend() : backend(shapes[1]))

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
    render_view(name)
  end

# Seletion

select_one_with_prompt(prompt::String, b::Backend, f::Function) =
  let ans = select_many_with_prompt(prompt, b, f)
    length(ans) > 0 ? ans[1] : nothing
  end

select_many_with_prompt(prompt::String, b::Backend, f::Function) =
  begin
    @info "$(prompt) on the $(b) backend."
    map(id -> shape_from_ref(id, b), f(connection(b), prompt))
  end

export render_view
render_view(name::String="View") =
  let path = prepare_for_saving_file(render_pathname(name))
    render_view(path, current_backend())
    path
  end
export save_view
save_view(name::String="View") =
  let path = prepare_for_saving_file(render_pathname(name))
    save_view(path, current_backend())
    path
  end

export realistic_sky
realistic_sky(;
    date::DateTime=DateTime(2020, 9, 21, 10, 0, 0),
    latitude::Real=39,
    longitude::Real=9,
    meridian::Real=0,
    altitude::Union{Missing,Real}=missing,
    azimuth::Union{Missing,Real}=missing,
    turbidity::Real=5,
    withsun::Bool=true) =
  ismissing(altitude) ?
    backend_realistic_sky(
      current_backend(),
      date, latitude, longitude, meridian, turbidity, withsun) :
    backend_realistic_sky(
      current_backend(),
      altitude, azimuth, turbidity, withsun)

export ground
ground(level::Loc=z(0), color::RGB=rgb(0.25,0.25,0.25)) =
  backend_ground(current_backend(), level, color)


############################################################
# Analysis

abstract type Analysis end
abstract type StructuralAnalysis <: Analysis end
abstract type LightingAnalysis <: Analysis end


###########################################################
# Geometric properties

# Axis-aligned Bounding Box

# Centroid

export centroid
centroid(s::Sphere) = s.center
centroid(s::Cylinder) = add_z(s.cb, s.h/2)
