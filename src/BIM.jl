
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

A third problem is related to level of detail. Even in the same backend, we might
need to represent BIM elements at different levels of detail.

We will start with the simplest of BIM elements, namely, levels, slabs, walls,
windows and doors.
=#

abstract type BIMElement <: Proxy end
abstract type Measure <: BIMElement end
BIMElements = Vector{<:BIMElement}

# Level
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

# LOD
abstract type LOD <: BIMElement end
struct LOD100 <: LOD end
#@defproxy(lod100, BIMElement, value::Integer=500)
#lod_cache = Dict{Integer,Level}()
#maybe_replace(lod::Lod) = get!(lod_cache, lod.value, lod)

#current_lods() = values(lod_cache)
#default_lod = Parameter{Lod}(lod())
#Base.:(==)(l1::Lod, l2::Lod) = l1.value == l2.value



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
  backend_wall(b::Backend, path, height, -thickness/2, thickness/2, family)

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
  backend_wall(b, translate(path, vz(bottom)), height, -thickness/2, thickness/2, getproperty(s.family, kind))

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


###################################
# BIM
@defop all_levels()
@defop all_walls()
@defop all_walls_at_level(level)
