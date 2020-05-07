using PlotlyJS

export Jupyter,
       jupyter

####################################################
abstract type JupyterKey end
const JupyterId = Any #PlotlyJS.GenericTrace
const JupyterRef = GenericRef{JupyterKey, JupyterId}
const JupyterNativeRef = NativeRef{JupyterKey, JupyterId}
const JupyterUnionRef = UnionRef{JupyterKey, JupyterId}
const JupyterSubtractionRef = SubtractionRef{JupyterKey, JupyterId}

const JupyterMaterial = String

mutable struct JupyterBackend{K,T} <: LazyBackend{K,T}
  shapes::Shapes
  shape_material::Dict{Shape,JupyterMaterial}
  materials::Dict{JupyterMaterial,JupyterMaterial}
  camera::Loc
  target::Loc
  lens::Real
end

const Jupyter = JupyterBackend{JupyterKey, JupyterId}
# Traits
has_boolean_ops(::Type{Jupyter}) = HasBooleanOps{false}()

const default_jupyter_material = Parameter{JupyterMaterial}("Red")

# In Jupyter, everytime we save a shape, we attach the default_jupyter_material
save_shape!(b::Jupyter, s::Shape) =
  begin
    prepend!(b.shapes, [s])
    b.shape_material[s] = default_jupyter_material()
    s
  end

#=
The Jupyter backend cannot realize shapes immediately, only when requested.
=#

void_ref(b::Jupyter) = JupyterNativeRef(-1)

const jupyter =
  Jupyter(Shape[],
         Dict{Shape,JupyterMaterial}(),
         Dict{JupyterMaterial,JupyterMaterial}(),
         xyz(10,10,10),
         xyz(0,0,0),
         35)

get_material(b::Jupyter, key) = get!(b.materials, key, key)
get_material(b::Jupyter, s::Shape) = get_material(b, get(b.shape_material, s, default_Jupyter_material()))

#

set_view(camera::Loc, target::Loc, lens::Real, b::Jupyter) =
  begin
    b.camera = camera
    b.target = target
    b.lens = lens
  end

get_view(b::Jupyter) =
  b.camera, b.target, b.lens

###################################
#
delete_all_shapes(b::Jupyter) =
  begin
    (empty!(b.shapes); empty!(b.materials); empty!(b.shape_material); nothing)
  end

backend_delete_shapes(b::Jupyter, shapes::Shapes) =
  begin
    b.shapes = filter(s->isnothing(findfirst(s1->s1===s, shapes)), b.shapes)
    for s in shapes delete!(b.shape_material, s) end
  end
#=
realize(b::Jupyter, s::Sphere) =
  let mat = get_material(b, s)
    write_Jupyter_object(buffer(b), "sphere", mat, in_world(s.center), s.radius)
    void_ref(b)
  end

realize(b::Jupyter, s::Torus) =
  let buf = buffer(b)
    write_Jupyter_object(buf, "torus", get_material(b, s), s.re, s.ri) do
      let p = in_world(s.center),
          t = s.center.cs.transform
        write_Jupyter_object(buf, "matrix", nothing,
                            [t[1,1], t[2,1], t[3,1],
                             t[1,2], t[2,2], t[3,2],
                             t[1,3], t[2,3], t[3,3],
                             p.x, p.y, p.z])
      end
    end
    void_ref(b)
  end

#=
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
  write_Jupyter_object(buffer(b), "prism", get_material(b, s), 0, norm(s.v)) do
    s.bs,
                              map(p -> (p + s.v), s.bs))
## FIXME: deal with the rotation angle
realize(b::ACAD, s::RightCuboid) =
  ACADCenteredBox(connection(b), s.cb, s.width, s.height, s.h)
=#

realize(b::Jupyter, s::Box) =
  let buf = buffer(b),
      bot = in_world(s.c),
      top = in_world(s.c + vxyz(s.dx, s.dy, s.dz, s.c.cs)),
      mat = get_material(b, s)
    write_Jupyter_object(buf, "box", mat, bot, top)
    void_ref(b)
  end

realize(b::Jupyter, s::Cone) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_Jupyter_object(buf, "cone", mat, bot, s.r, top, 0)
    void_ref(b)
  end

realize(b::Jupyter, s::ConeFrustum) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_Jupyter_object(buf, "cone", mat, bot, s.rb, top, s.rt)
    void_ref(b)
  end

realize(b::Jupyter, s::Cylinder) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_Jupyter_object(buf, "cylinder", mat, bot, top, s.r)
    void_ref(b)
  end

write_Jupyter_mesh(buf::IO, mat, points, closed_u, closed_v, smooth_u, smooth_v) =
  let si = size(points, 1),
      sj = size(points, 2),
      pts = in_world.(points),
      vcs = in_world.(add_z.(points, 1) .- points),
      idxs = quad_grid_indexes(si, sj, closed_u, closed_v)
#    let ps = reshape(permutedims(pts), :)
#      delete_all_shapes(autocad)
#      for tr in idxs
#        surface_polygon(ps[map(x->x+1,tr)], backend=autocad)
#      end
#    end
    write_Jupyter_object(buf, "mesh2", mat) do
      write_Jupyter_object(buf, "vertex_vectors", nothing, si*sj, reshape(permutedims(pts), :)...)
      # Must understand how to handle smoothness along one direction
      if smooth_u && smooth_v
        write_Jupyter_object(buf, "normal_vectors", nothing, si*sj, reshape(permutedims(vcs), :)...)
      end
      write_Jupyter_object(buf, "face_indices", nothing, length(idxs), idxs...)
      # Must understand how to handle smoothness along one direction
      #write_Jupyter_object(buf, "normal_indices", nothing, length(idxs), idxs...)
    end
  end
=#
realize(b::Jupyter, s::SurfaceGrid) =
  let mat = 1, #get_material(b, s)
      pts = map_division(in_world, s, size(s.points,1)-1, size(s.points,2)-1),
      xs = map(cx, map(r->r[1], pts)),
      ys = map(cy, pts[1]),
      zs = map(r->map(cz, r), pts)
    JupyterNativeRef(PlotlyJS.surface(x=xs, y=ys, z=zs))
  end

#=
realize(b::Jupyter, s::SweepPath) =
  let vertices = in_world.(path_vertices(s.profile)),
      frames = map_division(identity, s.path, 20),
      buf = buffer(b),
      mat = get_material(b, s)
    write_Jupyter_mesh(
      buf,
      mat,
      [xyz(cx(p), cy(p), cz(p), frame.cs) for p in vertices, frame in frames],
      is_closed_path(s.profile),
      is_closed_path(s.path),
      is_smooth_path(s.profile),
      is_smooth_path(s.path))
    void_ref(b)
  end

=#
# HACK: JUST FOR TESTING
realize(b::Jupyter, s::Thicken) =
  realize(b, s.shape)
#=






realize(b::Jupyter, s::EmptyShape) = void_ref(b)
realize(b::Jupyter, s::UniversalShape) = void_ref(b)
#=
realize(b::Jupyter, s::Move) =
    let r = map_ref(s.shape) do r
                JupyterMove(connection(b), r, s.v)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::Jupyter, s::Scale) =
    let r = map_ref(s.shape) do r
                JupyterScale(connection(b), r, s.p, s.s)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::Jupyter, s::Rotate) =
    let r = map_ref(s.shape) do r
                JupyterRotate(connection(b), r, s.p, s.v, s.angle)
                r
            end
        mark_deleted(s.shape)
        r
    end

=#

realize(b::Jupyter, s::UnionShape) =
  let shapes = filter(! is_empty_shape, s.shapes)
    length(shapes) == 1 ?
      (ref(shapes[1]); delete_shape(shapes[1])) :
      write_Jupyter_object(buffer(b), "union", get_material(b, s)) do
        for ss in shapes
          ref(ss)
          delete_shape(ss)
        end
      end
    void_ref(b)
  end


realize(b::Jupyter, s::IntersectionShape) =
  write_Jupyter_object(buffer(b), "intersection", get_material(b, s)) do
    for ss in s.shapes
      ref(ss)
      delete_shape(ss)
    end
    void_ref(b)
  end

realize(b::Jupyter, s::SubtractionShape3D) =
  write_Jupyter_object(buffer(b), "difference", get_material(b, s)) do
    ref(s.shape)
    delete_shape(s.shape)
    for ss in s.shapes
      ref(ss)
      delete_shape(ss)
    end
    void_ref(b)
  end

# BIM

realize_box(b::Jupyter, mat, p, dx, dy, dz) =
  let buf = buffer(b),
      bot = in_world(p),
      top = in_world(add_xyz(p, dx, dy, dz))
    write_Jupyter_object(buf, "box", mat, bot, top)
    void_ref(b)
  end

realize_prism(b::Jupyter, top, bot, side, path::PathSet, h::Real) =
  # PathSets require a different approach
  let buf = buffer(b),
      bot_vss = map(path_vertices, path.paths),
      top_vss = map(path_vertices, translate(path, vz(h)).paths)
    write_Jupyter_polygons(buf, bot, map(reverse, bot_vss))
    write_Jupyter_polygons(buf, top, top_vss)
    for (bot_vs, top_vs) in zip(bot_vss, top_vss)
      for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
        write_Jupyter_polygon(buf, side, vs)
      end
    end
    void_ref(b)
  end

realize_pyramid_frustum(b::Jupyter, top, bot, side, bot_vs::Locs, top_vs::Locs, closed=true) =
  let buf = buffer(b)
    if closed
      write_Jupyter_polygon(buf, bot, reverse(bot_vs))
      write_Jupyter_polygon(buf, top, top_vs)
    end
    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
      write_Jupyter_polygon(buf, side, vs)
    end
    void_ref(b)
  end

realize_polygon(b::Jupyter, mat, path::PathSet, acw=true) =
  acw ?
    write_Jupyter_polygons(buffer(b), mat, map(path_vertices, path.paths)) :
    write_Jupyter_polygons(buffer(b), mat, map(reverse ∘ path_vertices, path.paths))

realize_polygon(b::Jupyter, mat, vs::Locs, acw=true) =
  let buf = buffer(b)
    write_Jupyter_polygon(buf, mat, acw ? vs : reverse(vs))
  end

write_Jupyter_polygon(io::IO, mat, vs) =
  write_Jupyter_object(io, "polygon", mat, length(vs)+1, vs..., vs[1])

write_Jupyter_polygons(io::IO, mat, vss) =
  write_Jupyter_object(io, "polygon", mat,
    mapreduce(length, +, vss) + length(vss),
    mapreduce(vs->[vs..., vs[1]], vcat, vss)...)

# Polygons with holes need a PathSets in Jupyter

subtract_paths(b::Jupyter, c_r_w_path::PathSet, c_l_w_path::PathSet, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path.paths..., c_r_op_path),
  path_set(c_l_w_path.paths..., c_l_op_path)

subtract_paths(b::Jupyter, c_r_w_path::Path, c_l_w_path::Path, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path, c_r_op_path),
  path_set(c_l_w_path, c_l_op_path)


#=
Jupyter families need to know the different kinds of materials
that go on each surface.
In some cases it might be the same material, but in others, such
as slabs, outside walls, etc, we will have different materials.
=#

const JupyterMaterialFamily = BackendMaterialFamily{JupyterMaterial}
Jupyter_material_family(mat::JupyterMaterial) =
  JupyterMaterialFamily(mat)

const JupyterSlabFamily = BackendSlabFamily{JupyterMaterial}
Jupyter_slab_family(top::JupyterMaterial, bot::JupyterMaterial=top, side::JupyterMaterial=bot) =
  JupyterSlabFamily(top, bot, side)

const JupyterRoofFamily = BackendRoofFamily{JupyterMaterial}
Jupyter_roof_family(top::JupyterMaterial, bot::JupyterMaterial=top, side::JupyterMaterial=bot) =
  JupyterRoofFamily(top, bot, side)

const JupyterWallFamily = BackendWallFamily{JupyterMaterial}
Jupyter_wall_family(right::JupyterMaterial, left::JupyterMaterial=right) =
  JupyterWallFamily(right, left)

export Jupyter_material_family,
       Jupyter_slab_family,
       Jupyter_roof_family,
       Jupyter_wall_family,
       default_Jupyter_material

# Layers
current_layer(b::Jupyter) =
  default_Jupyter_material()

current_layer(layer, b::Jupyter) =
  default_Jupyter_material(layer)

backend_create_layer(b::Jupyter, name::String, active::Bool, color::RGB) =
  begin
    @assert active
    Jupyter_material(name, red=red(color), green=green(color), blue=blue(color))
  end

#=
create_ground_plane(shapes, material=default_Jupyter_ground_material()) =
  if shapes == []
    error("No shapes selected for analysis. Use add-Jupyter-shape!.")
  else
    let (p0, p1) = bounding_box(union(shapes)),
        (center, ratio) = (quad_center(p0, p1, p2, p3),
                  distance(p0, p4)/distance(p0, p2));
     ratio == 0 ?
      error("Couldn"t compute height. Use add-Jupyter-shape!.") :
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

#=
#FIXME define the family parameters for beams
realize(b::Jupyter, s::Beam) =
    ref(right_cuboid(s.p0, 0.2, 0.2, s.p1, 0))
=#
realize(b::Jupyter, s::Union{Door, Window}) =
  void_ref(b)

used_materials(b::Jupyter) =
  unique(map(f -> realize(s.family, b), b.shapes))

####################################################
=#
export export_to_jupyter
export_to_jupyter(b::Jupyter=current_backend()) =
  let x = 1
    PlotlyJS.plot([ref(s).value for s in b.shapes])
  end
