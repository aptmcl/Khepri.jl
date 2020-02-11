#=
A backend for POVRay
=#
export POVRay,
       povray,
       default_povray_material

# We will use MIME types to encode for POVRay
const MIMEPOVRay = MIME"text/povray"

import Base.show
# First, basic types
show(io::IO, ::MIMEPOVRay, s::String) =
  show(io, s)
show(io::IO, ::MIMEPOVRay, r::Real) =
  show(io, r)
show(io::IO, mime::MIMEPOVRay, v::Vector) =
  begin
    write(io, "<")
    for (i, e) in enumerate(v)
      if i > 1
        write(io, ", ")
      end
      show(io, mime, e)
    end
    write(io, ">")
    nothing
  end
show(io::IO, ::MIMEPOVRay, p::Union{Loc,Vec}) =
  # swap y with z to make it consistent with POVRay coordinate system
  print(io, "<$(p.x), $(p.z), $(p.y)>")
show(io::IO, ::MIMEPOVRay, c::RGB) =
  print(io, "color rgb <$(Float64(red(c))), $(Float64(green(c))), $(Float64(blue(c)))>")
show(io::IO, ::MIMEPOVRay, c::RGBA) =
  print(io, "color rgbt <$(Float64(red(c))), $(Float64(green(c))), $(Float64(blue(c))), $(Float64(1-alpha(c)))>")

write_povray_object(f::Function, io::IO, type, material, args...) =
  let mime = MIMEPOVRay()
    write(io, "$(type) {")
    for (i, arg) in enumerate(args)
      if i == 1
        write(io, "\n  ")
      elseif i > 1
        write(io, ", ")
      end
      show(io, mime, arg)
    end
    write(io, '\n')
    f()
    if ! isnothing(material)
      write_povray_material(io, material)
    end
    write(io, "}\n")
  end

write_povray_object(io::IO, type, material, args...) =
  write_povray_object(io, type, material, args...) do
  end

write_povray_param(io::IO, name::String, value::Any) =
  begin
    write(io, "  ", name, " ")
    show(io, MIMEPOVRay(), value)
    write(io, '\n')
  end

write_povray_material(io::IO, material) =
  begin
    write(io, "  ")
    show(io, MIMEPOVRay(), material)
    write(io, '\n')
  end

write_povray_camera(io::IO, camera, target, lens) =
  write_povray_object(io, "camera", nothing) do
    write_povray_param(io, "location", camera)
    write_povray_param(io, "look_at", target)
    write_povray_param(io, "angle", view_angles(lens)[1])
  end

write_povray_pointlight(io::IO, location, color) =
  write_povray_object(io, "light_source", nothing, location, color)


#=
struct POVRayPigment
  color::RGB
end

show(io::IO, mime::MIMEPOVRay, p::POVRayPigment) =
  write_povray_object(io, mime, "pigment") do
    show(io, mime, p.color)
  end

#=
NORMAL:
  normal { [NORMAL_IDENTIFIER] [NORMAL_TYPE] [NORMAL_MODIFIER...] }
NORMAL_TYPE:
  PATTERN_TYPE Amount |
  bump_map { BITMAP_TYPE "bitmap.ext" [BUMP_MAP_MODS...]}
NORMAL_MODIFIER:
  PATTERN_MODIFIER | NORMAL_LIST | normal_map { NORMAL_MAP_BODY } |
  slope_map{ SLOPE_MAP_BODY } | bump_size Amount |
  no_bump_scale Bool | accuracy Float
=#

struct POVRayNormal
  type::String
  amount::Real
end

povray_normal(; bumps) =
  POVRayNormal("bumps", bumps)

show(io::IO, mime::MIMEPOVRay, n::POVRayNormal) =
  write_povray_object(io, mime, "normal") do
    write_povray_param(io, mime, n.type, n.amount)
  end

#=
FINISH:
  finish { [FINISH_IDENTIFIER] [FINISH_ITEMS...] }
FINISH_ITEMS:
  ambient COLOR | diffuse [albedo] Amount [, Amount] | emission COLOR |
  brilliance Amount | phong [albedo] Amount | phong_size Amount | specular [albedo] Amount |
  roughness Amount | metallic [Amount] | reflection COLOR |
  crand Amount | conserve_energy BOOL_ON_OFF |
  reflection { Color_Reflecting_Min [REFLECTION_ITEMS...] } |
  subsurface { translucency COLOR } |
  irid { Irid_Amount [IRID_ITEMS...] }
REFLECTION_ITEMS:
  COLOR_REFLECTION_MAX | fresnel BOOL_ON_OFF |
  falloff FLOAT_FALLOFF | exponent FLOAT_EXPONENT |
  metallic FLOAT_METALLIC
IRID_ITEMS:
  thickness Amount | turbulence Amount
=#

struct POVRayFinish
  ambient::Union{RGB,Nothing}
  diffuse::Union{Real,Nothing}
  emission::Union{RGB,Nothing}
  brilliance::Union{Real,Nothing}
  phong::Union{Real,Nothing}
  phong_size::Union{Real,Nothing}
  specular::Union{Real,Nothing}
  roughness::Union{Real,Nothing}
  metallic::Union{Real,Nothing}
  reflection::Union{RGB,Nothing}
  crand::Union{Real,Nothing}
  conserve_energy::Union{Bool,Nothing}
  subsurface::Any
  irid::Any
end

povray_finish(; ambient::Union{RGB,Nothing}=nothing,
                diffuse::Union{Real,Nothing}=nothing,
                emission::Union{RGB,Nothing}=nothing,
                brilliance::Union{Real,Nothing}=nothing,
                phong::Union{Real,Nothing}=nothing,
                phong_size::Union{Real,Nothing}=nothing,
                specular::Union{Real,Nothing}=nothing,
                roughness::Union{Real,Nothing}=nothing,
                metallic::Union{Real,Nothing}=nothing,
                reflection::Union{RGB,Nothing}=nothing,
                crand::Union{Real,Nothing}=nothing,
                conserve_energy::Union{Bool,Nothing}=nothing,
                subsurface::Any=nothing,
                irid::Any=nothing) =
  POVRayFinish(ambient,
               diffuse,
               emission,
               brilliance,
               phong,
               phong_size,
               specular,
               roughness,
               metallic,
               reflection,
               crand,
               conserve_energy,
               subsurface,
               irid)

show(io::IO, mime::MIMEPOVRay, f::POVRayFinish) =
  write_povray_object(io, mime, "finish") do
    for name in fieldnames(typeof(f))
      let v = getfield(f, name)
        if ! isnothing(v)
          write_povray_param(io, mime, string(name), v)
        end
      end
    end
  end

struct POVRayTexture
  pygment::Union{POVRayPigment,Nothing}
  normal::Union{POVRayNormal,Nothing}
  finish::Union{POVRayFinish,Nothing}
end

struct POVRayLibraryTexture
  name::String
end

show(io::IO, mime::MIMEPOVRay, m::POVRayLibraryTexture) =
  write_povray_object(io, mime, "texture") do
    show(io, mime, m.name)
    #show(io, mime, m.interior)
  end

show(io::IO, mime::MIMEPOVRay, t::POVRayTexture) =
  write_povray_object(io, mime, "texture") do
    show(io, mime, t.pygment)
    show(io, mime, t.normal)
    show(io, mime, t.finish)
  end

struct POVRayInterior
  pygment::Union{POVRayPigment,Nothing}
  normal::Union{POVRayNormal,Nothing}
  finish::Union{POVRayFinish,Nothing}
end

struct POVRayMaterial
  name::String
  texture::POVRayTexture
  interior::Union{POVRayInterior,Nothing}
end
=#

abstract type POVRayMaterial end

struct POVRayDefinition <: POVRayMaterial
  name::String
  kind::String
  description::String
end

write_povray_definition(io::IO, d::POVRayDefinition) =
  write(io, "#declare $(d.name) =\n  $(d.kind) $(d.description)\n")

show(io::IO, ::MIMEPOVRay, d::POVRayDefinition) =
  write(io, "$(d.kind) { $(d.name) }")

struct POVRayInclude <: POVRayMaterial
  filename::AbstractString
  kind::String
  name::String
end

write_povray_definition(io::IO, m::POVRayInclude) =
  write(io, "#include \"$(m.filename)\"\n")

show(io::IO, ::MIMEPOVRay, m::POVRayInclude) =
  write(io, "$(m.kind) { $(m.name) }")

export povray_definition, povray_include
const povray_definition = POVRayDefinition
const povray_include = POVRayInclude

const texture_concrete =
  povray_definition("Concrete", "texture", """{
   pigment {
     granite turbulence 1.5 color_map {
       [0  .25 color White color Gray75] [.25  .5 color White color Gray75]
       [.5 .75 color White color Gray75] [.75 1.1 color White color Gray75]}}
  finish {
    ambient 0.2 diffuse 0.3 crand 0.03 reflection 0 }
  normal {
    dents .5 scale .5 }}""")

const default_povray_material =
  #Parameter(povray_include("stones.inc", "texture", "T_Grnt25"))
  Parameter(texture_concrete)

####################################################
# Sky models

####################################################

abstract type POVRayKey end
const POVRayId = Int
const POVRayRef = GenericRef{POVRayKey, POVRayId}
const POVRayNativeRef = NativeRef{POVRayKey, POVRayId}
const POVRayUnionRef = UnionRef{POVRayKey, POVRayId}
const POVRaySubtractionRef = SubtractionRef{POVRayKey, POVRayId}

mutable struct POVRayBackend{K,T} <: LazyBackend{K,T}
  shapes::Shapes
  shape_material::Dict{Shape,POVRayAbstractMaterial}
  sky::String
  buffer::LazyParameter{IOBuffer}
  materials::Dict
  camera::Loc
  target::Loc
  lens::Real
  sun_altitude::Real
  sun_azimuth::Real
end

const POVRay = POVRayBackend{POVRayKey, POVRayId}

# In POVRay, everytime we save a shape, we attach the default_povray_material
save_shape!(b::POVRayBackend, s::Shape) =
  begin
    push!(b.shapes, s)
    b.shape_material[s] = default_povray_material()
    s
  end

#=
The POVRay backend cannot realize shapes immediately, only when requested.
=#

void_ref(b::POVRay) = POVRayNativeRef(-1)

const povray =
  POVRay(Shape[],
         Dict{Shape,POVRayAbstractMaterial}(),
         "",
         LazyParameter(IOBuffer, IOBuffer),
         Dict(),
         xyz(10,10,10),
         xyz(0,0,0),
         35,
         90,
         0)

buffer(b::POVRay) = b.buffer()
get_material(b::POVRay, key) = get!(b.materials, key, key)
get_material(b::POVRay, s::Shape) = get_material(b, get(b.shape_material, s, default_povray_material()))

# Sky

const povray_normal_sky = """
sky_sphere{
 pigment{ gradient <0,1,0>
          color_map{
          [0.0 color rgb<1,1,1>        ]
          [0.8 color rgb<0.1,0.25,0.75>]
          [1.0 color rgb<0.1,0.25,0.75>]}
        } // end pigment
 }
 """

set_normal_sky(b::POVRay) =
  b.sky = povray_normal_sky

#

set_view(camera::Loc, target::Loc, lens::Real, b::POVRay) =
  begin
    set_view(camera, target, lens, autocad)
    b.camera = camera
    b.target = target
    b.lens = lens
  end

get_view(b::POVRay) =
  b.camera, b.target, b.lens

set_sun(altitude, azimuth, b::POVRay) =
  begin
    b.altitude = altitude
    b.azimuth = azimuth
  end

delete_all_shapes(b::POVRay) =
  begin
    delete_all_shapes(autocad)
    (empty!(b.shapes); empty!(b.materials); empty!(b.shape_material); nothing)
  end

realize(b::POVRay, s::Sphere) =
  let mat = get_material(b, s)
    write_povray_object(buffer(b), "sphere", mat, in_world(s.center), s.radius)
    void_ref(b)
  end

#=
realize(b::ACAD, s::Torus) =
  ACADTorus(connection(b), s.center, vz(1, s.center.cs), s.re, s.ri)
realize(b::ACAD, s::Cuboid) =
  ACADIrregularPyramidFrustum(connection(b), [s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3])
realize(b::ACAD, s::RegularPyramidFrustum) =
    ACADIrregularPyramidFrustum(connection(b),
                                regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                                regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed))
realize(b::ACAD, s::RegularPyramid) =
  ACADIrregularPyramid(connection(b),
                          regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                          add_z(s.cb, s.h))
realize(b::ACAD, s::IrregularPyramid) =
  ACADIrregularPyramid(connection(b), s.bs, s.t)
realize(b::ACAD, s::RegularPrism) =
  let ps = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
    ACADIrregularPyramidFrustum(connection(b),
                                   ps,
                                   map(p -> add_z(p, s.h), ps))
  end
realize(b::ACAD, s::IrregularPyramidFrustum) =
    ACADIrregularPyramidFrustum(connection(b), s.bs, s.ts)

realize(b::ACAD, s::IrregularPrism) =
  ACADIrregularPyramidFrustum(connection(b),
                              s.bs,
                              map(p -> (p + s.v), s.bs))
## FIXME: deal with the rotation angle
realize(b::ACAD, s::RightCuboid) =
  ACADCenteredBox(connection(b), s.cb, s.width, s.height, s.h)
=#

realize(b::POVRay, s::Box) =
  let buf = buffer(b),
      bot = in_world(s.c),
      top = in_world(s.c + vxyz(s.dx, s.dy, s.dz, s.c.cs)),
      mat = get_material(b, s)
    write_povray_object(buf, "box", mat, bot, top)
    void_ref(b)
  end

#=
realize(b::ACAD, s::Cone) =
  ACADCone(connection(b), add_z(s.cb, s.h), s.r, s.cb)
realize(b::ACAD, s::ConeFrustum) =
  ACADConeFrustum(connection(b), s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt)
=#

realize(b::POVRay, s::Cylinder) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_povray_object(buf, "cylinder", mat, bot, top, s.r)
    void_ref(b)
  end
#=
#=
realize(b::POVRay, s::EmptyShape) =
    EmptyRef{POVRayId}()
realize(b::POVRay, s::UniversalShape) =
    UniversalRef{POVRayId}()

realize(b::POVRay, s::Move) =
    let r = map_ref(s.shape) do r
                POVRayMove(connection(b), r, s.v)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::POVRay, s::Scale) =
    let r = map_ref(s.shape) do r
                POVRayScale(connection(b), r, s.p, s.s)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::POVRay, s::Rotate) =
    let r = map_ref(s.shape) do r
                POVRayRotate(connection(b), r, s.p, s.v, s.angle)
                r
            end
        mark_deleted(s.shape)
        r
    end

=#
=#
# BIM

realize_pyramid_fustrum(b::POVRay, top, bot, side, bot_vs::Locs, top_vs::Locs, closed=true) =
  let buf = buffer(b)
    if closed
      write_povray_polygon(buf, bot, reverse(bot_vs))
      write_povray_polygon(buf, top, top_vs)
    end
    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
      write_povray_polygon(buf, side, vs)
    end
  end

realize_polygon(b::POVRay, path::Path, acw=true) =
  realize_polygon(b, path_vertices(path), acw)

realize_polygon(b::POVRay, vs, acw=true) =
  let buf = buffer(b)
    polygon(vs, backend=autocad)
    write_povray_polygon(buf, acw ? vs : reverse(vs))
  end

write_povray_polygon(io::IO, mat, vs) =
  write_povray_object(buf, "polygon", mat, length(vs), vs...)


#=
POVRay families need to know the different kinds of materials
that go on each surface.
In some cases it might be the same material, but in others, such
as slabs, outside walls, etc, we will have different materials.
=#

const POVRayMaterialFamily = BackendMaterialFamily{POVRayMaterial}
povray_material_family(mat::POVRayMaterial) =
  POVRayMaterialFamily(mat)

const POVRaySlabFamily = BackendSlabFamily{POVRayMaterial}
povray_slab_family(top::POVRayMaterial, bot::POVRayMaterial=top, side::POVRayMaterial=bot) =
  POVRaySlabFamily(top, bot, side)

const POVRayFoofFamily = BackendRoofFamily{POVRayMaterial}
povray_roof_family(top::POVRayMaterial, bot::POVRayMaterial=top, side::POVRayMaterial=bot) =
  POVRayRoofFamily(top, bot, side)

const POVRayWallFamily = BackendWallFamily{POVRayMaterial}
povray_wall_family(right::POVRayMaterial, left::POVRayMaterial=right) =
  POVRayWallFamily(right, left)

export povray_material_family,
       povray_slab_family,
       povray_roof_family,
       povray_wall_family,
       default_povray_material
#=
set_backend_family(default_wall_family(), povray, povray_include_texture())
set_backend_family(default_slab_family(), povray,
  povray_slab_family(povray_generic_floor_20, povray_generic_ceiling_80))
set_backend_family(default_roof_family(), povray,
  povray_roof_family(povray_generic_floor_20, povray_outside_facade_30))
set_backend_family(default_beam_family(), povray,
  povray_material_family(povray_generic_metal))
set_backend_family(default_column_family(), povray,
  povray_material_family(povray_generic_metal))
set_backend_family(default_door_family(), povray,
  povray_material_family(povray_generic_furniture_50))
set_backend_family(default_panel_family(), povray,
  povray_material_family(povray_glass_material("GenericGlass", gray=0.3)))
=#
#=
create_ground_plane(shapes, material=default_povray_ground_material()) =
  if shapes == []
    error("No shapes selected for analysis. Use add-povray-shape!.")
  else
    let (p0, p1) = bounding_box(union(shapes)),
        (center, ratio) = (quad_center(p0, p1, p2, p3),
                  distance(p0, p4)/distance(p0, p2));
     ratio == 0 ?
      error("Couldn"t compute height. Use add-povray-shape!.") :
      let pts = map(p -> intermediate_loc(center, p, ratio*10), [p0, p1, p2, p3]);
         create_surface_layer(pts, 0, ground_layer(), material)
        end
       end
  end
        w = max(floor_extra_factor()*distance(p0, p1), floor_extra_width())
        with(current_layer,floor_layer()) do
          box(xyz(min(p0.x, p1.x)-w, min(p0.y, p1.y)-w, p0.z-1-floor_distance()),
              xyz(max(p0.x, p1.x)+w, max(p0.y, p1.y)+w, p0.z-0-floor_distance()))
        end
      end
    end

=#

backend_slab(b::POVRay, profile, openings, thickness, family) =
  let mattop = realize(b, family).top_material,
      matbot = realize(b, family).bottom_material,
      matside = realize(b, family).side_material,
      path = profile
    for op in openings
      path = subtract_paths(path, op)
    end
    realize_pyramid_fustrum(
      b, mattop, matbot, matside,
      path_vertices(path),
      path_vertices(translate(path, vz(thickness))))
  end

#=
#=
#FIXME define the family parameters for beams
realize(b::POVRay, s::Beam) =
    ref(right_cuboid(s.p0, 0.2, 0.2, s.p1, 0))

=#
realize(b::POVRay, s::Panel) =
  let p1 = s.vertices[1],
      p2 = s.vertices[2],
      p3 = s.vertices[3],
      n = vz(s.family.thickness/2, cs_from_o_vx_vy(p1, p2-p1, p3-p1)),
      mod = next_modifier(b, realize(b, s.family).material)
    realize_pyramid_fustrum(
        b, mod, mod, mod,
        map(p -> in_world(p - n), s.vertices),
        map(p -> in_world(p + n), s.vertices))
  end

closed_path_for_height(path, h) =
  let ps = path_vertices(path)
    closed_polygonal_path([ps..., reverse(map(p -> p+vz(h), ps))...])
  end

#=
One important restriction is that POVRay only supports _planar_ polygons.
This forces us to create multiple subpaths.
=#
realize(b::POVRay, w::Wall) =
  let w_base_height = w.bottom_level.height,
      w_height = w.top_level.height - w_base_height,
      r_thickness = r_thickness(w),
      l_thickness = l_thickness(w),
      w_path = translate(w.path, vz(w_base_height)),
      w_paths = subpaths(w_path),
      r_w_paths = subpaths(offset(w_path, -r_thickness)),
      l_w_paths = subpaths(offset(w_path, l_thickness)),
      openings = [w.doors..., w.windows...],
      prevlength = 0,
      modright = next_modifier(b, realize(b, w.family).right_material),
      modleft = next_modifier(b, realize(b, w.family).left_material)
    for (w_seg_path, r_w_path, l_w_path) in zip(w_paths, r_w_paths, l_w_paths)
      let currlength = prevlength + path_length(w_seg_path),
          c_r_w_path = closed_path_for_height(r_w_path, w_height),
          c_l_w_path = closed_path_for_height(l_w_path, w_height)
        realize_pyramid_fustrum(b, modleft, modright, modright, c_l_w_path, c_r_w_path, false)
        openings = filter(openings) do op
          if prevlength <= op.loc.x < currlength ||
             prevlength <= op.loc.x + op.family.width <= currlength # contained (at least, partially)
            let op_height = op.family.height,
                op_at_start = op.loc.x <= prevlength,
                op_at_end = op.loc.x + op.family.width >= currlength,
                op_path = subpath(w_path,
                                  max(prevlength, op.loc.x),
                                  min(currlength, op.loc.x + op.family.width)),
                r_op_path = offset(op_path, -r_thickness),
                l_op_path = offset(op_path, l_thickness),
                fixed_r_op_path =
                  open_polygonal_path([path_start(op_at_start ? r_w_path : r_op_path),
                                       path_end(op_at_end ? r_w_path : r_op_path)]),
                fixed_l_op_path =
                  open_polygonal_path([path_start(op_at_start ? l_w_path : l_op_path),
                                       path_end(op_at_end ? l_w_path : l_op_path)]),
                c_r_op_path = closed_path_for_height(translate(fixed_r_op_path, vz(op.loc.y)), op_height),
                c_l_op_path = closed_path_for_height(translate(fixed_l_op_path, vz(op.loc.y)), op_height),
                idxs = closest_vertices_indexes(path_vertices(c_r_w_path), path_vertices(c_r_op_path))
              realize_pyramid_fustrum(b, modleft, modright, modright, c_r_op_path, c_l_op_path, false)
              c_r_w_path =
                closed_polygonal_path(
                  inject_polygon_vertices_at_indexes(path_vertices(c_r_w_path), path_vertices(c_r_op_path), idxs))
              c_l_w_path =
                closed_polygonal_path(
                  inject_polygon_vertices_at_indexes(path_vertices(c_l_w_path), path_vertices(c_l_op_path), idxs))
              # preserve if not totally contained
              ! (op.loc.x >= prevlength && op.loc.x + op.family.width <= currlength)
            end
          else
            true
          end
        end
        prevlength = currlength
        realize_polygon(b, modleft, "wall", c_l_w_path, false)
        realize_polygon(b, modright, "wall", c_r_w_path, true)
      end
    end
    void_ref(b)
  end

realize(b::POVRay, w::Window) = nothing
realize(b::POVRay, w::Door) = nothing

backend_wall(b::POVRay, w_path, w_height, l_thickness, r_thickness, family) =
  error("BUM")
#=  let w_paths = subpaths(w_path),
      r_w_paths = subpaths(offset(w_path, -r_thickness)),
      l_w_paths = subpaths(offset(w_path, l_thickness)),
      modright = next_modifier(b, realize(b, s.family).material_right),
      modleft = next_modifier(b, realize(b, s.family).material_left)
    for (w_seg_path, r_w_path, l_w_path) in zip(w_paths, r_w_paths, l_w_paths)
      let c_r_w_path = closed_path_for_height(r_w_path, w_height),
          c_l_w_path = closed_path_for_height(l_w_path, w_height)
        realize_pyramid_fustrum(b, modleft, modright, modright, c_l_w_path, c_r_w_path, false)
        # This is surely wrong!
        #realize_polygon(b, family, "wall", c_l_w_path, false)
        #realize_polygon(b, family, "wall", c_r_w_path, true)
      end
    end
    void_ref(b)
  end
=#
=#

set_backend_family(default_wall_family(), povray, acad_layer_family("Wall"))
set_backend_family(default_slab_family(), povray, acad_layer_family("Slab"))
set_backend_family(default_roof_family(), povray, acad_layer_family("Roof"))
set_backend_family(default_beam_family(), povray, acad_layer_family("Beam"))
set_backend_family(default_column_family(), povray, acad_layer_family("Column"))
set_backend_family(default_door_family(), povray, acad_layer_family("Door"))
set_backend_family(default_panel_family(), povray, acad_layer_family("Panel"))

set_backend_family(default_table_family(), povray, acad_layer_family("Table"))
set_backend_family(default_chair_family(), povray, acad_layer_family("Chair"))
set_backend_family(default_table_chair_family(), povray, acad_layer_family("TableChairs"))

set_backend_family(default_curtain_wall_family(), povray, acad_layer_family("CurtainWall"))
set_backend_family(default_curtain_wall_family().panel, povray, acad_layer_family("CurtainWall-Panel"))
set_backend_family(default_curtain_wall_family().boundary_frame, povray, acad_layer_family("CurtainWall-Boundary"))
set_backend_family(default_curtain_wall_family().transom_frame, povray, acad_layer_family("CurtainWall-Transom"))
set_backend_family(default_curtain_wall_family().mullion_frame, povray, acad_layer_family("CurtainWall-Mullion"))
#current_backend(povray)

set_backend_family(default_truss_node_family(), povray, acad_layer_family("TrussNodes"))
set_backend_family(default_truss_bar_family(), povray, acad_layer_family("TrussBars"))





set_sun_direction(v::Vec, b::POVRay) =
  b.sun_direction = v


add_ground_plane(b::POVRay) =
  @warn "Not generating ground plane"

used_materials(b::POVRay) =
  unique(map(f -> realize(s.family, b), b.shapes))

####################################################

export_to_povray(path::String, b::POVRay=current_backend()) =
  let buf = b.buffer()
    # First pass, to fill material dictionary
    take!(buf)
    for s in b.shapes
      realize(b, s)
    end
    open(path, "w") do out
      # write materials
      write_povray_include(out, POVRayIncludeMaterial("colors.inc", "dummy"))
      for (k,v) in b.materials
        write_povray_include(out, k)
      end
      # write sky
      write(out, b.sky, "\n")
      # write the objects
      write(out, String(take!(buf)))
      # write the view
      write_povray_camera(out, b.camera, b.target, b.lens)
    end
  end

#=
Simulations need to be done on a temporary folder, so that we can have multiple
simulations running at the same time.
=#

povray_simulation_path() = joinpath(string(@__DIR__), "FOO.pov")
  #=
  let (path, io) = mktemp(mktempdir(tempdir(), prefix="POVRay_"))
    close(io)
    path
  end
 =#
export povray_folder
const povray_folder = Parameter("C:/Program Files/POV-Ray/v3.7/bin/")

povray_cmd(cmd::AbstractString="pvengine64") = povray_folder() * cmd

##########################################

export povray_render
povray_render() =
  let path = povray_simulation_path()
    export_to_povray(path)
    #run(`$(povray_cmd()) +w +h $path`)
    run(`$(povray_cmd()) Width=$(render_width()) Height=$(render_height()) /RENDER $path`, wait=false)
    #run(`$(povray_cmd("wxFalseColor")) $picpath`, wait=false)
  end
