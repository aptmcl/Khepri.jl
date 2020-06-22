export meshcat

#=
ATTENTION!!!

Use Chrome debugger to visualize threeJS examples, break the execution and
evaluate in the console:

JSON.stringify(scene.toJSON(), undefined, 2)

=#



## Primitives
using UUIDs, MsgPack

send_meshcat(vis::Visualizer, obj) =
  write(vis.core, pack(obj))

send_setobject(vis, path, obj) =
  let msg = (type="set_object", path=path, object=obj),
      data = pack(msg)
    vis.core.tree[path].object = data
    write(vis.core, data)
    path
  end

send_settransform(vis, path, transform) =
  let msg = (type="set_transform", path=path, matrix=convert(Vector{Float32}, transform)),
      data = pack(msg)
    vis.core.tree[path].transform = data
    write(vis.core, data)
    nothing
  end

send_setproperty(vis, path, property, value) =
  let msg = (type="set_property", path=path, property=property, value=value),
      data = pack(msg)
    vis.core.tree[path].properties[property] = data
    write(vis.core, data)
    nothing
  end

send_setvisible(vis, path, value) =
  send_meshcat(vis, (type="set_property", path=path, property="visible", value=value))

#send_setproperty(connection(meshcat), "/Background/hide_background", "value", true)
#send_meshcat(connection(meshcat), (type="hide_background",))

send_delobject(vis, path) =
  send_meshcat(vis, (type="delete", path=path))

meshcat_point(p::Loc) =
  let p = in_world(p)
    [cx(p),cz(p),cy(p)]
  end

meshcat_transform(p::Loc) =
  # MeshCat swaps Y with Z
  let m = translated_cs(p.cs, p.x, p.y, p.z).transform*[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    m[:]
  end

meshcat_metadata() =
  (type="Object", version=4.5)

meshcat_material(color) =
  (uuid=string(uuid1()),
   type="MeshLambertMaterial",
    #"vertexColors" => 0,
    #"transparent" => false,
    #"opacity" => 1.0,
    #"depthTest" => true,
    #"linewidth" => 1.0,
    #"depthFunc" => 3,
    side=2,
    color="0x$(hex(color))",
    #color="0xAAAAAA",
    #"reflectivity" => 0.5,
    #depthWrite=true
    )

meshcat_line_material(color) =
  (uuid=string(uuid1()),
   type="LineBasicMaterial",
   color="0x$(hex(color))",
   #linewidth=2, Due to limitations of the OpenGL Core Profile with the WebGL renderer on most platforms linewidth will always be 1 regardless of the set value.
   #depthFunc=3,
   #depthTest=true,
   #depthWrite=true,
   #stencilWrite=false,
   #stencilWriteMask=255,
   #stencilFunc=519,
   #stencilRef=0,
   #stencilFuncMask=255,
   #stencilFail=7680,
   #stencilZFail=7680,
   #stencilZPass=7680
   )

meshcat_object(type, geom, material, p=u0(world_cs)) =
  (metadata=meshcat_metadata(),
   geometries=[geom],
   materials=[material],
   object=(uuid=string(uuid1()),
           type=type,
           geometry=geom.uuid,
           material=material.uuid,
           matrix=meshcat_transform(p)))


meshcat_buffer_geometry_attributes_position(vertices) =
  (itemSize=3,
   type="Float32Array",
   array=convert(Vector{Float32}, reduce(vcat, [meshcat_point(v) for v in vertices])))

meshcat_line(vertices, material) =
  let geom = (uuid=string(uuid1()),
              type="BufferGeometry",
              data=(
                attributes=(
                  position=meshcat_buffer_geometry_attributes_position(vertices),),))
    meshcat_object("Line", geom, material)
  end

meshcat_surface_polygon(vertices, material) =
  let geom = (uuid=string(uuid1()),
              type="BufferGeometry",
              data=(
                attributes=(
                  position=meshcat_buffer_geometry_attributes_position([vertices..., vertices[1]]),),))
    meshcat_object("Mesh", geom, material)
  end

meshcat_circle_mesh(center, radius, start_angle, amplitude, material) =
  let geom = (uuid=string(uuid1()),
              type="CircleBufferGeometry",
              radius=radius,
              segments=64,
              thetaStart=start_angle,
              thetaLength=amplitude)
    meshcat_object("Mesh", geom, material, center)
  end

#=



=#

meshcat_centered_box(p, dx, dy, dz, material) =
  let geom = (uuid=string(uuid1()),
              type="BoxGeometry",
              width=dx,
              depth=dy,
              height=dz)
    meshcat_object("Mesh", geom, material, p)
  end

meshcat_box(p, dx, dy, dz, material) =
  meshcat_centered_box(p+vxyz(dx/2, dy/2, dz/2, p.cs), dx, dy, dz, material)

#=
send_delobject(v, "/meshcat")
send_setobject(v, "/meshcat/box1", meshcat_box(xyz(1,1,1), 1, 2, 3))
send_setobject(v, "/meshcat/box2", meshcat_box(loc_from_o_vx_vy(u0(), vxy(1,1), vxy(-1,1)), 1, 2, 3))
send_delobject(v, "/meshcat/box2")
=#

meshcat_sphere(p, r, material) =
  let geom = (uuid=string(uuid1()),
              type="SphereGeometry",
              radius=r,
              widthSegments=64,
              heightSegments=64)
    meshcat_object("Mesh", geom, material, p)
  end

meshcat_centered_cone(p, r, h, material) =
  let geom = (uuid=string(uuid1()),
              type="ConeGeometry",
              radius=r,
              height=h,
              radialSegments=64)
    meshcat_object("Mesh", geom, material, p)
  end

meshcat_cone(p, r, h, material) =
  meshcat_centered_cone(p+vz(h/2, p.cs), r, h, material)

meshcat_centered_cone_frustum(p, rb, rt, h, material) =
  let geom = (uuid=string(uuid1()),
              type="CylinderGeometry",
              radiusTop=rt,
              radiusBottom=rb,
              height=h,
              radialSegments=64)
    meshcat_object("Mesh", geom, material, p)
  end

meshcat_cone_frustum(p, rb, rt, h, material) =
  meshcat_centered_cone_frustum(p+vz(h/2, p.cs), rb, rt, h, material)

meshcat_cylinder(p, r, h, material) =
  meshcat_centered_cone_frustum(p+vz(h/2, p.cs), r, r, h, material)

#=
send_setobject(v, "/meshcat/sphere1", meshcat_sphere(xyz(-2,-3,0), 1))
send_setobject(v, "/meshcat/cylinder1", meshcat_cylinder(loc_from_o_vx_vy(x(-3), vxy(1,1), vxz(-1,1)), 1, 5))
=#
meshcat_mesh(vertices, faces, material) =
  let geom = (uuid=string(uuid1()),
              type="BufferGeometry",
              data=(
                attributes=(
                  position=meshcat_buffer_geometry_attributes_position(vertices),
                  #uv=?
                  ),
                index=(
                  itemSize=3,
                  type="Uint32Array",
                  array=convert(Vector{UInt32}, collect(Iterators.flatten(faces))))))
    meshcat_object("Mesh", geom, material, u0(world_cs))
  end

####################################################
abstract type MCATKey end
const MCATId = String #MCATlyJS.GenericTrace
const MCATRef = GenericRef{MCATKey, MCATId}
const MCATNativeRef = NativeRef{MCATKey, MCATId}
const MCATUnionRef = UnionRef{MCATKey, MCATId}
const MCATSubtractionRef = SubtractionRef{MCATKey, MCATId}

abstract type MCATMaterial end

struct MCATLayer
  name
  material
  line_material
end

mcat_layer(name, color) =
  MCATLayer(name, meshcat_material(color), meshcat_line_material(color))

mutable struct MCATBackend{K,T} <: Backend{K,T}
  connection::LazyParameter{Visualizer}
  count::Int64
  layer::MCATLayer
end

const MCAT = MCATBackend{MCATKey, MCATId}

connection(b::MCAT) = b.connection()
material(b::MCAT) = b.layer.material
line_material(b::MCAT) = b.layer.line_material

create_MCAT_connection() =
  let (width, height) = render_size(),
      vis = Visualizer()
    display(render(vis))
    MeshCat.wait(vis)
    vis
  end

meshcat = MCAT(LazyParameter(Visualizer, create_MCAT_connection),
               0,
               mcat_layer("default", RGB(1,1,1)))

export display_view
display_view() =
  let vis = connection(meshcat)
    display(render(vis))
    MeshCat.wait(vis)
  end

const meshcat_root_path = "/Khepri"

next_id(b::MCATBackend{K,T}) where {K,T} =
  begin
      b.count += 1
      string(meshcat_root_path, "/", b.layer.name, "/", b.count - 1)
  end

add_object(b::MCAT, obj) =
  send_setobject(connection(b), next_id(b), obj)

has_boolean_ops(::Type{MCAT}) = HasBooleanOps{false}()
void_ref(b::MCAT) = MCATNativeRef("")

reset_backend(b::MCAT) =
  display(render(connection(b)))

new_backend(b::MCAT) =
  reset(b.connection)
#

#=
setprop!(
  connection(meshcat)["/Cameras/default/rotated/<object>"], #"Cameras","default","rotated","<object>"],
  "zoom",
  1.0)
=#

#=
function MeshCat.lower(s::Sphere)
  error("BUM")
    Dict{String, Any}(
        "uuid" => string(uuid1()),
        "type" => "SphereGeometry",
        "radius" => s.radius,
        "widthSegments" => 20,
        "heightSegments" => 20,
    )
end

add_object!(meshcat, s)
=#
#=
MeshCat.setcontrol!(
  connection(meshcat)["Background"], #"Cameras","default","rotated","<object>"],
  "hide_background")

{
    type: "set_property",
    path: "/Cameras/default/rotated/<object>",
    property: "zoom",
    value: 2.0
}
=#

set_view(camera::Loc, target::Loc, lens::Real, aperture::Real, b::MCAT) =
  let v = connection(b),
      (x1,y1,z1) = raw_point(camera),
      (x2,y2,z2) = raw_point(target)
    send_settransform(v, "/Cameras/default", [
      1, 0, 0, 1,
      0, 1, 0, 0,
      0, 0, 1, 0,
      x2, y2, z2, 1])
    send_setproperty(v, "/Cameras/default/rotated/<object>", "position", [x1-x2,z1-z2,y2-y1])
    #  b.lens = lens
  end

send_setproperty(connection(meshcat), "/Cameras/default/rotated/<object>", "position", [10,0,10])
send_setproperty(connection(meshcat), "/Cameras/default/rotated/<object>", "quaternion", [pi/2,0,0,0])



get_view(b::MCAT) =
  b.camera, b.target, b.lens

###################################
#
delete_all_shapes(b::MCAT) =
  let c = connection(b)
    send_delobject(c, meshcat_root_path)
  end

backend_delete_shapes(b::MCAT, shapes::Shapes) =
  let c = connection(b)
    for r in collect_ref(shapes)
      send_delobject(c, r)
    end
  end

backend_stroke(b::MCAT, path::CircularPath) =
  add_object(b, meshcat_line(path_frames(path), line_material(b)))

backend_stroke(b::MCAT, path::RectangularPath) =
  let c = path.corner,
      dx = path.dx,
      dy = path.dy
    add_object(
      b,
      meshcat_line([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)],
                   line_material(b)))
  end

backend_stroke(b::MCAT, path::ArcPath) =
  add_object(
    b,
    meshcat_line(
      path_frames(arc_path(path.center, path.radius, path.start_angle, path.amplitude)),
      line_material(b)))

backend_stroke(b::MCAT, path::OpenPolygonalPath) =
	add_object(b, meshcat_line(path.vertices, line_material(b)))

backend_stroke(b::MCAT, path::ClosedPolygonalPath) =
  add_object(b, meshcat_line([path.vertices...,path.vertices[1]], line_material(b)))
  #=

backend_fill(b::MCAT, path::ClosedPolygonalPath) =
  add_object(b, meshcat_surface_polygon(path.vertices, material(b)))


backend_fill(b::MCAT, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        @remote(b, SurfaceClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)]))
    end
backend_stroke(b::MCAT, path::OpenSplinePath) =
  if (path.v0 == false) && (path.v1 == false)
    #@remote(b, Spline(path.vertices))
    @remote(b, InterpSpline(
                     path.vertices,
                     path.vertices[2]-path.vertices[1],
                     path.vertices[end]-path.vertices[end-1]))
  elseif (path.v0 != false) && (path.v1 != false)
    @remote(b, InterpSpline(path.vertices, path.v0, path.v1))
  else
    @remote(b, InterpSpline(
                     path.vertices,
                     path.v0 == false ? path.vertices[2]-path.vertices[1] : path.v0,
                     path.v1 == false ? path.vertices[end-1]-path.vertices[end] : path.v1))
  end
backend_stroke(b::MCAT, path::ClosedSplinePath) =
    @remote(b, InterpClosedSpline(path.vertices))
backend_fill(b::MCAT, path::ClosedSplinePath) =
    backend_fill_curves(b, @remote(b, InterpClosedSpline(path.vertices)))

backend_fill_curves(b::MCAT, refs::MCATIds) = @remote(b, SurfaceFromCurves(refs))
backend_fill_curves(b::MCAT, ref::MCATId) = @remote(b, SurfaceFromCurves([ref]))

backend_stroke_line(b::MCAT, vs) = @remote(b, PolyLine(vs))

backend_stroke_arc(b::MCAT, center::Loc, radius::Real, start_angle::Real, amplitude::Real) =
  let end_angle = start_angle + amplitude
    @remote(b, Arc(center, vz(1, center.cs), radius, start_angle, end_angle))
  end
backend_stroke_unite(b::MCAT, refs) = @remote(b, JoinCurves(refs))

realize(b::MCAT, s::Point) =
  @remote(b, Point(s.position))
=#
realize(b::MCAT, s::Line) =
  stroke(open_polygonal_path(s.vertices), b)

realize(b::MCAT, s::Spline) =
 # This should be merged with opensplinepath
  if (s.v0 == false) && (s.v1 == false)
    add_object(b, meshcat_line(path_frames(open_spline_path(s.points)), line_material(b)))
  elseif (s.v0 != false) && (s.v1 != false)
    error("Finish this")
  else
    error("Finish this")
  end

#=
realize(b::MCAT, s::ClosedSpline) =
  @remote(b, InterpClosedSpline(s.points))
realize(b::MCAT, s::Circle) =
  @remote(b, Circle(s.center, vz(1, s.center.cs), s.radius))
realize(b::MCAT, s::Arc) =
  if s.radius == 0
    @remote(b, Point(s.center))
  elseif s.amplitude == 0
    @remote(b, Point(s.center + vpol(s.radius, s.start_angle, s.center.cs)))
  elseif abs(s.amplitude) >= 2*pi
    @remote(b, Circle(s.center, vz(1, s.center.cs), s.radius))
  else
    end_angle = s.start_angle + s.amplitude
    if end_angle > s.start_angle
      @remote(b, Arc(s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle))
    else
      @remote(b, Arc(s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle))
    end
  end

realize(b::MCAT, s::Ellipse) =
  if s.radius_x > s.radius_y
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x))
  else
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y))
  end
realize(b::MCAT, s::EllipticArc) =
  error("Finish this")

realize(b::MCAT, s::Polygon) =
  @remote(b, ClosedPolyLine(s.vertices))
realize(b::MCAT, s::RegularPolygon) =
  @remote(b, ClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))
realize(b::MCAT, s::Rectangle) =
  @remote(b, ClosedPolyLine(
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)]))
realize(b::MCAT, s::SurfaceCircle) =
  @remote(b, SurfaceCircle(s.center, vz(1, s.center.cs), s.radius))
realize(b::MCAT, s::SurfaceArc) =
    #@remote(b, SurfaceArc(s.center, vz(1, s.center.cs), s.radius, s.start_angle, s.start_angle + s.amplitude))
    if s.radius == 0
        @remote(b, Point(s.center))
    elseif s.amplitude == 0
        @remote(b, Point(s.center + vpol(s.radius, s.start_angle, s.center.cs)))
    elseif abs(s.amplitude) >= 2*pi
        @remote(b, SurfaceCircle(s.center, vz(1, s.center.cs), s.radius))
    else
        end_angle = s.start_angle + s.amplitude
        if end_angle > s.start_angle
            @remote(b, SurfaceFromCurves(
                [@remote(b, Arc(s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle)),
                 @remote(b, PolyLine([add_pol(s.center, s.radius, end_angle),
                                              add_pol(s.center, s.radius, s.start_angle)]))]))
        else
            @remote(b, SurfaceFromCurves(
                [@remote(b, Arc(s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle)),
                 @remote(b, PolyLine([add_pol(s.center, s.radius, s.start_angle),
                                              add_pol(s.center, s.radius, end_angle)]))]))
        end
    end

realize(b::MCAT, s::SurfaceEllipse) =
  if s.radius_x > s.radius_y
    @remote(b, SurfaceEllipse(s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x))
  else
    @remote(b, SurfaceEllipse(s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y))
  end


realize(b::MCAT, s::SurfacePolygon) =
  @remote(b, SurfaceClosedPolyLine(s.vertices))
realize(b::MCAT, s::SurfaceRegularPolygon) =
  @remote(b, SurfaceClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))
realize(b::MCAT, s::SurfaceRectangle) =
  @remote(b, SurfaceClosedPolyLine(
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)]))
realize(b::MCAT, s::Surface) =
  let #ids = map(r->@remote(b, NurbSurfaceFrom(r)), @remote(b, SurfaceFromCurves(collect_ref(s.frontier))))
      ids = @remote(b, SurfaceFromCurves(collect_ref(s.frontier)))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::MCAT, s::Shape2D) =
    map(c -> shape_from_ref(c, b), @remote(b, CurvesFromSurface(ref(s).value)))
=#

realize(b::MCAT, s::Sphere) =
  add_object(b, meshcat_sphere(s.center, s.radius, material(b)))
#=
realize(b::MCAT, s::Torus) =
  let buf = buffer(b)
    write_MCAT_object(buf, "torus", get_material(b, s), s.re, s.ri) do
      let p = in_world(s.center),
          t = s.center.cs.transform
        write_MCAT_object(buf, "matrix", nothing,
                            [t[1,1], t[2,1], t[3,1],
                             t[1,2], t[2,2], t[3,2],
                             t[1,3], t[2,3], t[3,3],
                             p.x, p.y, p.z])
      end
    end
    void_ref(b)
  end

#=
realize(b::MCAT, s::Cuboid) =
  MCATIrregularPyramidFrustum(connection(b), [s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3])
realize(b::MCAT, s::RegularPyramidFrustum) =
    MCATIrregularPyramidFrustum(connection(b),
                                regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                                regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed))
realize(b::MCAT, s::RegularPyramid) =
  MCATIrregularPyramid(connection(b),
                          regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                          add_z(s.cb, s.h))
realize(b::MCAT, s::IrregularPyramid) =
  MCATIrregularPyramid(connection(b), s.bs, s.t)
realize(b::MCAT, s::RegularPrism) =
  let ps = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
    MCATIrregularPyramidFrustum(connection(b),
                                   ps,
                                   map(p -> add_z(p, s.h), ps))
  end
realize(b::MCAT, s::IrregularPyramidFrustum) =
    MCATIrregularPyramidFrustum(connection(b), s.bs, s.ts)

realize(b::MCAT, s::IrregularPrism) =
  write_MCAT_object(buffer(b), "prism", get_material(b, s), 0, norm(s.v)) do
    s.bs,
                              map(p -> (p + s.v), s.bs))
## FIXME: deal with the rotation angle
=#
=#
realize(b::MCAT, s::RightCuboid) =
  add_object(b, meshcat_box(s.cb, s.width, s.height, s.h, material(b)))

realize(b::MCAT, s::Box) =
  let #mat = get_material(b, s)
    add_object(b, meshcat_box(s.c, s.dx, s.dy, s.dz, material(b)))
  end

realize(b::MCAT, s::Cone) =
  add_object(b, meshcat_cone(s.cb, s.r, s.h, material(b)))

realize(b::MCAT, s::ConeFrustum) =
  add_object(b, meshcat_cone_frustum(s.cb, s.rb, s.rt, s.h, material(b)))

realize(b::MCAT, s::Cylinder) =
  send_setobject(connection(b), next_id(b), meshcat_cylinder(s.cb, s.r, s.h, material(b)))

meshcat_faces(si, sj, closed_u, closed_v) =
  let idx(i,j) = (i-1)*sj+(j-1),
      idxs = [],
      quad(a,b,c,d) = (push!(idxs, (a, b, d)); push!(idxs, (d, b, c)))
    for i in 1:si-1
      for j in 1:sj-1
        quad(idx(i,j), idx(i+1,j), idx(i+1,j+1), idx(i,j+1))
      end
      if closed_v
        quad(idx(i,sj), idx(i+1,sj), idx(i+1,1), idx(i,1))
      end
    end
    if closed_u
      for j in 1:sj-1
        quad(idx(si,j), idx(1,j), idx(1,j+1), idx(si,j+1))
      end
      if closed_v
        quad(idx(si,sj), idx(1,sj), idx(1,1), idx(si,1))
      end
    end
    idxs
  end

realize(b::MCAT, s::SurfaceGrid) =
  let pts = in_world.(s.points),
      si = size(pts, 1),
      sj = size(pts, 2),
      idxs = meshcat_faces(si, sj, s.closed_u, s.closed_v)
    send_setobject(connection(b), next_id(b), meshcat_mesh(reshape(permutedims(pts),:), idxs, material(b))) #reshape(permutedims(pts), :), idxs))
  end

#=
realize(b::MCAT, s::SweepPath) =
  let vertices = in_world.(path_vertices(s.profile)),
      frames = map_division(identity, s.path, 20),
      buf = buffer(b),
      mat = get_material(b, s)
    write_MCAT_mesh(
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
realize(b::MCAT, s::Thicken) =
  realize(b, s.shape)

realize(b::MCAT, s::EmptyShape) = void_ref(b)

realize(b::MCAT, s::UniversalShape) = void_ref(b)

#=
#=
realize(b::MCAT, s::Move) =
    let r = map_ref(s.shape) do r
                MCATMove(connection(b), r, s.v)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::MCAT, s::Scale) =
    let r = map_ref(s.shape) do r
                MCATScale(connection(b), r, s.p, s.s)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::MCAT, s::Rotate) =
    let r = map_ref(s.shape) do r
                MCATRotate(connection(b), r, s.p, s.v, s.angle)
                r
            end
        mark_deleted(s.shape)
        r
    end

=#

realize(b::MCAT, s::UnionShape) =
  let shapes = filter(! is_empty_shape, s.shapes)
    length(shapes) == 1 ?
      (ref(shapes[1]); delete_shape(shapes[1])) :
      write_MCAT_object(buffer(b), "union", get_material(b, s)) do
        for ss in shapes
          ref(ss)
          delete_shape(ss)
        end
      end
    void_ref(b)
  end


realize(b::MCAT, s::IntersectionShape) =
  write_MCAT_object(buffer(b), "intersection", get_material(b, s)) do
    for ss in s.shapes
      ref(ss)
      delete_shape(ss)
    end
    void_ref(b)
  end

realize(b::MCAT, s::SubtractionShape3D) =
  write_MCAT_object(buffer(b), "difference", get_material(b, s)) do
    ref(s.shape)
    delete_shape(s.shape)
    for ss in s.shapes
      ref(ss)
      delete_shape(ss)
    end
    void_ref(b)
  end
=#
# BIM

realize_box(b::MCAT, mat, p, dx, dy, dz) =
  add_object(b, meshcat_box(p, dx, dy, dz, mat))

#=
realize_prism(b::MCAT, top, bot, side, path::PathSet, h::Real) =
  # PathSets require a different approach
  let buf = buffer(b),
      bot_vss = map(path_vertices, path.paths),
      top_vss = map(path_vertices, translate(path, vz(h)).paths)
    write_MCAT_polygons(buf, bot, map(reverse, bot_vss))
    write_MCAT_polygons(buf, top, top_vss)
    for (bot_vs, top_vs) in zip(bot_vss, top_vss)
      for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
        write_MCAT_polygon(buf, side, vs)
      end
    end
    void_ref(b)
  end

realize_pyramid_frustum(b::MCAT, top, bot, side, bot_vs::Locs, top_vs::Locs, closed=true) =
  let buf = buffer(b)
    if closed
      write_MCAT_polygon(buf, bot, reverse(bot_vs))
      write_MCAT_polygon(buf, top, top_vs)
    end
    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
      write_MCAT_polygon(buf, side, vs)
    end
    void_ref(b)
  end

realize_polygon(b::MCAT, mat, path::PathSet, acw=true) =
  acw ?
    write_MCAT_polygons(buffer(b), mat, map(path_vertices, path.paths)) :
    write_MCAT_polygons(buffer(b), mat, map(reverse ∘ path_vertices, path.paths))

realize_polygon(b::MCAT, mat, vs::Locs, acw=true) =
  let buf = buffer(b)
    write_MCAT_polygon(buf, mat, acw ? vs : reverse(vs))
  end

write_MCAT_polygon(io::IO, mat, vs) =
  write_MCAT_object(io, "polygon", mat, length(vs)+1, vs..., vs[1])

write_MCAT_polygons(io::IO, mat, vss) =
  write_MCAT_object(io, "polygon", mat,
    mapreduce(length, +, vss) + length(vss),
    mapreduce(vs->[vs..., vs[1]], vcat, vss)...)

# Polygons with holes need a PathSets in MCAT

subtract_paths(b::MCAT, c_r_w_path::PathSet, c_l_w_path::PathSet, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path.paths..., c_r_op_path),
  path_set(c_l_w_path.paths..., c_l_op_path)

subtract_paths(b::MCAT, c_r_w_path::Path, c_l_w_path::Path, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path, c_r_op_path),
  path_set(c_l_w_path, c_l_op_path)


#=
MCAT families need to know the different kinds of materials
that go on each surface.
In some cases it might be the same material, but in others, such
as slabs, outside walls, etc, we will have different materials.
=#

const MCATMaterialFamily = BackendMaterialFamily{MCATMaterial}
MCAT_material_family(mat::MCATMaterial) =
  MCATMaterialFamily(mat)

const MCATSlabFamily = BackendSlabFamily{MCATMaterial}
MCAT_slab_family(top::MCATMaterial, bot::MCATMaterial=top, side::MCATMaterial=bot) =
  MCATSlabFamily(top, bot, side)

const MCATRoofFamily = BackendRoofFamily{MCATMaterial}
MCAT_roof_family(top::MCATMaterial, bot::MCATMaterial=top, side::MCATMaterial=bot) =
  MCATRoofFamily(top, bot, side)

const MCATWallFamily = BackendWallFamily{MCATMaterial}
MCAT_wall_family(right::MCATMaterial, left::MCATMaterial=right) =
  MCATWallFamily(right, left)

export MCAT_material_family,
       MCAT_slab_family,
       MCAT_roof_family,
       MCAT_wall_family,
       default_MCAT_material

=#

# Layers
current_layer(b::MCAT) =
  b.layer

current_layer(layer, b::MCAT) =
  b.layer = layer

backend_create_layer(b::MCAT, name::String, active::Bool, color::RGB) =
  begin
    @assert active
    mcat_layer(name, color)
  end

#=
create_ground_plane(shapes, material=default_MCAT_ground_material()) =
  if shapes == []
    error("No shapes selected for analysis. Use add-MCAT-shape!.")
  else
    let (p0, p1) = bounding_box(union(shapes)),
        (center, ratio) = (quad_center(p0, p1, p2, p3),
                  distance(p0, p4)/distance(p0, p2));
     ratio == 0 ?
      error("Couldn"t compute height. Use add-MCAT-shape!.") :
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
realize(b::MCAT, s::Beam) =
    ref(right_cuboid(s.p0, 0.2, 0.2, s.p1, 0))
=#
realize(b::MCAT, s::Union{Door, Window}) =
  void_ref(b)

####################################################

#=
set_backend_family(default_wall_family(), meshcat,
  radiance_wall_family(radiance_generic_interior_wall_70))
set_backend_family(default_slab_family(), meshcat,
  radiance_slab_family(radiance_generic_floor_20, radiance_generic_ceiling_80))
set_backend_family(default_roof_family(), meshcat,
  radiance_roof_family(radiance_generic_floor_20, radiance_outside_fMCATe_30))
set_backend_family(default_beam_family(), meshcat,
  radiance_material_family(radiance_generic_metal))
set_backend_family(default_column_family(), meshcat,
  radiance_material_family(radiance_generic_metal))
set_backend_family(default_door_family(), meshcat,
  radiance_material_family(radiance_generic_furniture_50))
set_backend_family(default_panel_family(), meshcat,
  radiance_material_family(radiance_glass_material("GenericGlass", gray=0.3)))
=#
