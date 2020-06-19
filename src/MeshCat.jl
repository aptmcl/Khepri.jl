#=
convert(::Type{GeometryTypes.Point}, p::Khepri.Loc) =
  let p = in_world(p)
    GeometryTypes.Point(Float64(cx(p)), Float64(cy(p)), Float64(cz(p)))
  end
convert(::Type{GeometryTypes.Vec}, v::Khepri.Vec) =
  let v = in_world(v)
    GeometryTypes.Vec(Float64(Khepri.cx(v)), Float64(Khepri.cy(v)), Float64(Khepri.cz(v)))
  end
convert(::Type{GeometryTypes.Normal}, v::Khepri.Vec) =
  let v = in_world(v)
    GeometryTypes.Normal(Float64(Khepri.cx(v)), Float64(Khepri.cy(v)), Float64(Khepri.cz(v)))
  end

add_object!(vis, id, obj) =
  begin

    id
  end
=#
export meshcat

## Primitives

send_meshcat(vis::Visualizer, obj) =
  write(vis.core, pack(obj))

send_setobject(vis, path, obj) =
  let msg = (type="set_object", path=path, object=obj)
    send_meshcat(vis, msg)
    path
  end

using UUIDs, MsgPack

send_delobject(vis, path) =
  send_meshcat(vis, (type="delete", path=path))

meshcat_transform(p::Loc) =
  # MeshCat swaps Y with Z
  let m = translated_cs(p.cs, p.x, p.y, p.z).transform*[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
    m[:]
  end

meshcat_metadata() =
  (type="Object", version=4.5)

meshcat_materials(material_id) =
  [(
    #"vertexColors" => 0,
    #"transparent" => false,
    #"opacity" => 1.0,
    #"depthTest" => true,
    #"linewidth" => 1.0,
    #"depthFunc" => 3,
    side=2,
    color="0xFFFFFF",
    #"reflectivity" => 0.5,
    type="MeshLambertMaterial",
    uuid=material_id,
    depthWrite=true
    )]

meshcat_centered_box_geometries(geom_id, dx, dy, dz) =
  [(type="BoxGeometry",
    uuid=geom_id,
    width=dx,
    depth=dy,
    height=dz)]

meshcat_object(geom_id, material_id, p) =
  (geometry=geom_id,
   material=material_id,
   matrix=meshcat_transform(p),
   type="Mesh",
   uuid=string(uuid1()))

meshcat_centered_box(p, dx, dy, dz) =
  let geom_id = string(uuid1()),
      material_id = string(uuid1())
    (metadata=meshcat_metadata(),
     geometries=meshcat_centered_box_geometries(geom_id, dx, dy, dz),
     materials=meshcat_materials(material_id),
     object=meshcat_object(geom_id, material_id, p))
  end

meshcat_box(p, dx, dy, dz) =
  meshcat_centered_box(p+vxyz(dx/2, dy/2, dz/2, p.cs), dx, dy, dz)


#=
send_delobject(v, "/meshcat")
send_setobject(v, "/meshcat/box1", meshcat_box(xyz(1,1,1), 1, 2, 3))
send_setobject(v, "/meshcat/box2", meshcat_box(loc_from_o_vx_vy(u0(), vxy(1,1), vxy(-1,1)), 1, 2, 3))
send_delobject(v, "/meshcat/box2")
=#

meshcat_sphere(p, r) =
  let geom_id = string(uuid1()),
      material_id = string(uuid1())
    (metadata=meshcat_metadata(),
     geometries=[(type="SphereGeometry",
                  uuid=geom_id,
                  radius=r,
                  widthSegments=64,
                  heightSegments=64)],
     materials=meshcat_materials(material_id),
     object=meshcat_object(geom_id, material_id, p))
  end

meshcat_centered_cone(p, r, h) =
  let geom_id = string(uuid1()),
      material_id = string(uuid1())
    (metadata=meshcat_metadata(),
     geometries=[(type="ConeGeometry",
                  uuid=geom_id,
                  radius=r,
                  height=h,
                  radialSegments=64)],
     materials=meshcat_materials(material_id),
     object=meshcat_object(geom_id, material_id, p))
  end

meshcat_cone(p, r, h) =
  meshcat_centered_cone(p+vz(h/2, p.cs), r, h)

meshcat_centered_cone_frustum(p, rb, rt, h) =
  let geom_id = string(uuid1()),
      material_id = string(uuid1())
    (metadata=meshcat_metadata(),
     geometries=[(type="CylinderGeometry",
                  uuid=geom_id,
                  radiusTop=rt,
                  radiusBottom=rb,
                  height=h,
                  radialSegments=64)],
     materials=meshcat_materials(material_id),
     object=meshcat_object(geom_id, material_id, p))
  end

meshcat_cone_frustum(p, rb, rt, h) =
  meshcat_centered_cone_frustum(p+vz(h/2, p.cs), rb, rt, h)

meshcat_cylinder(p, r, h) =
  meshcat_centered_cone_frustum(p+vz(h/2, p.cs), r, r, h)

#=
send_setobject(v, "/meshcat/sphere1", meshcat_sphere(xyz(-2,-3,0), 1))
send_setobject(v, "/meshcat/cylinder1", meshcat_cylinder(loc_from_o_vx_vy(x(-3), vxy(1,1), vxz(-1,1)), 1, 5))
=#

swap_y_z(v) = [v[1],v[3],v[2]]

meshcat_mesh(vertices, faces) =
  let geom_id = string(uuid1()),
      material_id = string(uuid1()),
      pts = (itemSize=3,
             type="Float32Array",
             array=convert(Vector{Float32}, reduce(vcat, [swap_y_z(v.raw) for v in vertices]))),
      idxs = (itemSize=3,
              type="Uint32Array",
              array=convert(Vector{UInt32}, collect(Iterators.flatten(faces))))
    (metadata=meshcat_metadata(),
     geometries=[(type="BufferGeometry",
                  uuid=geom_id,
                  data=(attributes=(position=pts,
                                    #uv=?
                                    ),
                        index=idxs))],
     materials=meshcat_materials(material_id),
     object=meshcat_object(geom_id, material_id, u0()))
  end

#=
Dict{String,Any}(
  "object" => Dict{String,Any}(
    "type" => "set_object",
    "path" => "/meshcat"
    "geometries" => Any[Dict{String,Any}(
      "uuid" => "1c892160-b1fe-11ea-3542-9fc0f2d5a2ea",
      "data" => Dict{String,Any}(
        "attributes" => Dict{String,Any}(
          "position" => Dict{String,Any}(
            "array" => Float32[-0.0001, -0.0001, -0.0001, 0.0, 0.0, -0.0001, -0.0001, 0.0, -0.0001, 0.0, -0.0001, 0.0, 0.0, -0.0001, -0.0001, -0.0001, 0.0, 0.0, -0.0001, -0.0001, 0.0, 0.0999, -0.0001, -0.1, 0.05, 0.0, -0.05, 0.06666667, -0.033333335, -0.033333335, 0.0999, 0.0, -0.1, 0.0999, -0.1, -0.0001, 0.1, -0.05, -0.05, 0.1, -0.1, -0.0001, 0.1, -0.0001, -0.1, 0.05, -0.05, 0.0, 0.0999, -0.1, 0.0, -0.0001, 0.0999, -0.1, -0.05, 0.1, -0.05, -0.033333335, 0.06666667, -0.033333335, -0.0001, 0.1, -0.1, 0.0, 0.05, -0.05, 0.0, 0.0999, -0.1, -0.1, 0.0999, -0.0001, -0.1, 0.1, -0.0001, -0.05, 0.05, 0.0, -0.1, 0.0999, 0.0, 0.05, 0.05, -0.1, 0.033333335, 0.033333335, -0.06666667, -0.05, 0.0, 0.05, -0.033333335, -0.033333335, 0.06666667, -0.0001, -0.1, 0.0999, 0.0, -0.05, 0.05, 0.0, -0.1, 0.0999, -0.1, -0.0001, 0.0999, -0.1, 0.0, 0.0999, -0.05, -0.05, 0.1, -0.0001, -0.1, 0.1, -0.1, -0.0001, 0.1, 0.033333335, -0.06666667, 0.033333335, 0.05, -0.1, 0.05, -0.06666667, 0.033333335, 0.033333335, -0.1, 0.05, 0.05],
            "type" => "Float32Array",
            "itemSize" => 3)),
        "index" => Dict{String,Any}(
          "array" => MeshCat.PackedVector{UInt32}(UInt32[0x00000000, 0x00000001, 0x00000002, 0x00000000, 0x00000003, 0x00000004, 0x00000000, 0x00000004, 0x00000001, 0x00000000, 0x00000002, 0x00000005, 0x00000000, 0x00000006, 0x00000003, 0x00000000, 0x00000005, 0x00000006, 0x00000007, 0x00000008, 0x00000009, 0x00000007, 0x0000000a, 0x00000008, 0x0000000b, 0x0000000c, 0x00000009, 0x0000000b, 0x0000000d, 0x0000000c, 0x00000007, 0x0000000c, 0x0000000e, 0x00000007, 0x00000009, 0x0000000c, 0x00000004, 0x00000008, 0x00000001, 0x00000004, 0x00000009, 0x00000008, 0x0000000b, 0x0000000f, 0x00000010, 0x0000000b, 0x00000009, 0x0000000f, 0x00000004, 0x0000000f, 0x00000009, 0x00000004, 0x00000003, 0x0000000f, 0x00000011, 0x00000012, 0x00000013, 0x00000011, 0x00000014, 0x00000012, 0x00000002, 0x00000015, 0x00000013, 0x00000002, 0x00000001, 0x00000015, 0x00000011, 0x00000015, 0x00000016, 0x00000011, 0x00000013, 0x00000015, 0x00000017, 0x00000012, 0x00000018, 0x00000017, 0x00000013, 0x00000012, 0x00000002, 0x00000019, 0x00000005, 0x00000002, 0x00000013, 0x00000019, 0x00000017, 0x00000019, 0x00000013, 0x00000017, 0x0000001a, 0x00000019, 0x0000001b, 0x00000016, 0x0000001c, 0x00000008, 0x0000000a, 0x0000001c, 0x0000000a, 0x0000001b, 0x0000001c, 0x00000016, 0x00000015, 0x0000001c, 0x00000001, 0x00000008, 0x0000001c, 0x00000015, 0x00000001, 0x0000001c, 0x00000006, 0x0000001d, 0x0000001e, 0x00000006, 0x00000005, 0x0000001d, 0x0000001f, 0x00000020, 0x0000001e, 0x0000001f, 0x00000021, 0x00000020, 0x00000006, 0x00000020, 0x00000003, 0x00000006, 0x0000001e, 0x00000020, 0x00000022, 0x0000001d, 0x00000023, 0x00000022, 0x0000001e, 0x0000001d, 0x0000001f, 0x00000024, 0x00000025, 0x0000001f, 0x0000001e, 0x00000024, 0x00000022, 0x00000024, 0x0000001e, 0x00000022, 0x00000026, 0x00000024, 0x0000000f, 0x00000003, 0x00000027, 0x00000028, 0x00000010, 0x00000027, 0x00000010, 0x0000000f, 0x00000027, 0x00000003, 0x00000020, 0x00000027, 0x00000021, 0x00000028, 0x00000027, 0x00000020, 0x00000021, 0x00000027, 0x00000019, 0x0000001a, 0x00000029, 0x0000001d, 0x00000005, 0x00000029, 0x00000005, 0x00000019, 0x00000029, 0x0000001a, 0x0000002a, 0x00000029, 0x00000023, 0x0000001d, 0x00000029, 0x0000002a, 0x00000023, 0x00000029]),
          "type" => "Uint32Array",
          "itemSize" => 3)),
      "type" => "BufferGeometry")],
    "object" => Dict{String,Any}(
      "material" => "1c892160-b1fe-11ea-1947-6f8e0ec55cac",
      "uuid" => "1c892160-b1fe-11ea-25b0-51dfca2aef35",
      "geometry" => "1c892160-b1fe-11ea-3542-9fc0f2d5a2ea",
      "matrix" => [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
      "type" => "Mesh"),
    "materials" => Any[Dict{String,Any}(
      "vertexColors" => 0,
      "transparent" => false,
      "opacity" => 1.0f0,
      "depthTest" => true,
      "linewidth" => 1.0,
      "depthFunc" => 3,
      "side" => 2,
      "uuid" => "1c892160-b1fe-11ea-1947-6f8e0ec55cac",
      "color" => "0xFFFFFF",
      "type" => "MeshLambertMaterial",
      "depthWrite" => true)],
    "metadata" => Dict{String,Any}(
      "version" => 4.5,"type" => "Object")))
=#

#=
create_meshcat_spline(path, vs) =
    let geom_id = string(uuid1()),
        material_id = string(uuid1()),
        obj_id = string(uuid1())
      Dict("type="set_object",
           "path=path,
           "object=Dict(
              "metadata"=> Dict{String, Any}(
                  "type="Object",
                  "version=4.5),
              "geometries=[Dict{String, Any}(
                  "type="SplineCurve",
                  "uuid=geom_id,
                  "points=[(-10,0), (-5,5), (0,0), (5,-5), (10,0)])],
              "materials=[Dict{String, Any}(
                  "vertexColors=0,
                  "transparent=false,
                  "opacity=1.0,
                  "depthTest=true,
                  "linewidth=1.0,
                  "depthFunc=3,
                  "side=2,
                  "color="0xFFFFFF",
                  "reflectivity=0.5,
                  "type="LineBasicMaterial",
                  "uuid=material_id,
                  "depthWrite=true)],
              "object"=> Dict("geometry=geom_id,
                              "material=material_id,
                              "matrix=[1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                              "type="Mesh",
                              "uuid=obj_id)))
    end

send_meshcat(v, create_meshcat_spline("/meshcat/box1", []))
=#
#=

  Drawing lines
  Let's say you want to draw a line or a circle, not a wireframe Mesh. First we need to set up the renderer, scene and camera (see the Creating a scene page).

  Here is the code that we will use:

  var renderer = new THREE.WebGLRenderer();
  renderer.setSize( window.innerWidth, window.innerHeight );
  document.body.appendChild( renderer.domElement );

  var camera = new THREE.PerspectiveCamera( 45, window.innerWidth / window.innerHeight, 1, 500 );
  camera.position.set( 0, 0, 100 );
  camera.lookAt( 0, 0, 0 );

  var scene = new THREE.Scene();
  Next thing we will do is define a material. For lines we have to use LineBasicMaterial or LineDashedMaterial.

  //create a blue LineBasicMaterial
  var material = new THREE.LineBasicMaterial( { color: 0x0000ff } );
  After material we will need a geometry with some vertices:

  var points = [];
  points.push( new THREE.Vector3( - 10, 0, 0 ) );
  points.push( new THREE.Vector3( 0, 10, 0 ) );
  points.push( new THREE.Vector3( 10, 0, 0 ) );

  var geometry = new THREE.BufferGeometry().setFromPoints( points );
  Note that lines are drawn between each consecutive pair of vertices, but not between the first and last (the line is not closed.)

  Now that we have points for two lines and a material, we can put them together to form a line.

  var line = new THREE.Line( geometry, material );
  All that's left is to add it to the scene and call render.

  scene.add( line );
  renderer.render( scene, camera );
=#

#=

{
    "type": "set_object",
    "path": "/meshcat/box",
    "object": {
        "metadata": {"type": "Object", "version": 4.5},
        "geometries": [{"depth": 0.5,
                        "height": 0.5,
                        "type": "BoxGeometry",
                        "uuid": "fbafc3d6-18f8-11e8-b16e-f8b156fe4628",
                        "width": 0.5}],
        "materials": [{"color": 16777215,
                       "reflectivity": 0.5,
                       "type": "MeshPhongMaterial",
                       "uuid": "e3c21698-18f8-11e8-b16e-f8b156fe4628"}],
        "object": {"geometry": "fbafc3d6-18f8-11e8-b16e-f8b156fe4628",
                   "material": "e3c21698-18f8-11e8-b16e-f8b156fe4628",
                   "matrix": [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
                   "type": "Mesh",
                   "uuid": "fbafc3d7-18f8-11e8-b16e-f8b156fe4628"}},
}
Translating that box by the vector [2, 3, 4]:

{
    "type": "set_transform",
    "path": "/meshcat/box",
    "matrix": [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 2.0, 3.0, 4.0, 1.0]
}

=#


####################################################
abstract type MCATKey end
const MCATId = String #MCATlyJS.GenericTrace
const MCATRef = GenericRef{MCATKey, MCATId}
const MCATNativeRef = NativeRef{MCATKey, MCATId}
const MCATUnionRef = UnionRef{MCATKey, MCATId}
const MCATSubtractionRef = SubtractionRef{MCATKey, MCATId}

abstract type MCATMaterial end

mutable struct MCATBackend{K,T} <: Backend{K,T}
  connection::LazyParameter{Visualizer}
  count::Int64
  layer::String
  material::Union{Nothing, MCATMaterial}
end

const MCAT = MCATBackend{MCATKey, MCATId}

create_MCAT_connection() =
  let (width, height) = render_size(),
      vis = Visualizer()
    display(render(vis))
    MeshCat.wait(vis)
    vis
  end

meshcat = MCAT(LazyParameter(Visualizer, create_MCAT_connection),
               0,
               "default",
               nothing)

connection(b::MCATBackend{K,T}) where {K,T} = b.connection()

const meshcat_root_path = "/Khepri"

next_id(b::MCATBackend{K,T}) where {K,T} =
  begin
      b.count += 1
      string(meshcat_root_path, "/", b.layer, "/", b.count - 1)
  end
#=
add_object!(b::MCAT, obj) =
  let vis = connection(b),
      id = next_id(b)
    isnothing(b.material) ?
      send_meshcat!(vis, joinpath(b.layer, id), obj) :
      error("Go handle material!") #setobject!(vis[b.layer][id], obj, b.material)
  end


function setobject!(vis::Visualizer, obj::AbstractObject)
    send(vis.core, SetObject(obj, vis.path))
    vis
end

function setobject!(vis::ArrowVisualizer, material::AbstractMaterial=defaultmaterial();
        shaft_material::AbstractMaterial=material,
        head_material::AbstractMaterial=material)
    settransform!(vis, zero(Point{3, Float64}), zero(Vec{3, Float64}))
    shaft = Cylinder(zero(Point{3, Float64}), Point(0.0, 0.0, 1.0), 1.0)
    setobject!(vis.vis[:shaft], shaft, shaft_material)
    head = Cone(zero(Point{3, Float64}), Point(0.0, 0.0, 1.0), 1.0)
    setobject!(vis.vis[:head], head, head_material)
    vis
end
=#
#=
For each object, we need to pack two messages, one to encode the object properties,
the other to encode its transform.
=#

has_boolean_ops(::Type{MCAT}) = HasBooleanOps{false}()
void_ref(b::MCAT) = MCATNativeRef("")

reset_backend(b::MCAT) =
  display(render(connection(b)))

new_backend(b::MCAT) =
  begin
    reset(b.connection)
    display(render(connection(b)))
  end
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

struct MCATColoredMaterial <: MCATMaterial
  color
end
=#

set_view(camera::Loc, target::Loc, lens::Real, aperture::Real, b::MCAT) =
  begin
    b.camera = camera
    b.target = target
    b.lens = lens
  end

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
  add_object!(
    b,
    GeometryTypes.Circle(convert(GeometryTypes.Point, path.center),
                         Float64(path.radius)))

    #=
backend_stroke(b::MCAT, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        @remote(b, ClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)]))
    end
backend_stroke(b::MCAT, path::ArcPath) =
    backend_stroke_arc(b, path.center, path.radius, path.start_angle, path.amplitude)

backend_stroke(b::MCAT, path::OpenPolygonalPath) =
  	@remote(b, PolyLine(path.vertices))
backend_stroke(b::MCAT, path::ClosedPolygonalPath) =
    @remote(b, ClosedPolyLine(path.vertices))
backend_fill(b::MCAT, path::ClosedPolygonalPath) =
    @remote(b, SurfaceClosedPolyLine(path.vertices))
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



realize(b::MCAT, s::EmptyShape) =
  MCATEmptyRef()
realize(b::MCAT, s::UniversalShape) =
  MCATUniversalRef()
realize(b::MCAT, s::Point) =
  @remote(b, Point(s.position))

realize(b::MCAT, s::Line) =
  let pts = map(in_world, s.vertices),
      r = MCATlyJS.scatter3d(
         x=map(cx, pts),
         y=map(cy, pts),
         z=map(cz, pts),
         line_shape="linear",
         marker_size=2,
         autocolorscale=false,
         showscale=false,
         hoverinfo="skip")
    MCATlyJS.addtraces!(connection(b), r)
    MCATNativeRef(r)
  end

realize(b::MCAT, s::Spline) = # This should be merged with opensplinepath
  if (s.v0 == false) && (s.v1 == false)
    let pts = map_division(in_world, spline_path(s.points), length(s.points)*4),
        r = MCATlyJS.scatter3d(
           x=map(cx, pts),
           y=map(cy, pts),
           z=map(cz, pts),
           line_shape="spline",
           marker_size=2,
           autocolorscale=false,
           showscale=false,
           hoverinfo="skip")
      MCATlyJS.addtraces!(connection(b), r)
      MCATNativeRef(r)
    end
  elseif (s.v0 != false) && (s.v1 != false)
    @remote(b, InterpSpline(s.points, s.v0, s.v1))
  else
    @remote(b, InterpSpline(
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1))
  end
  =#

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
  send_setobject(connection(b), next_id(b), meshcat_sphere(s.center, s.radius))
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
realize(b::MCAT, s::RightCuboid) =
  MCATCenteredBox(connection(b), s.cb, s.width, s.height, s.h)
=#

realize(b::MCAT, s::Box) =
  let buf = buffer(b),
      bot = in_world(s.c),
      top = in_world(s.c + vxyz(s.dx, s.dy, s.dz, s.c.cs)),
      mat = get_material(b, s)
    write_MCAT_object(buf, "box", mat, bot, top)
    void_ref(b)
  end
=#

realize(b::MCAT, s::Cone) =
  send_setobject(connection(b), next_id(b), meshcat_cone(s.cb, s.r, s.h))

#=
realize(b::MCAT, s::ConeFrustum) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_MCAT_object(buf, "cone", mat, bot, s.rb, top, s.rt)
    void_ref(b)
  end

=#

realize(b::MCAT, s::Cylinder) =
  send_setobject(connection(b), next_id(b), meshcat_cylinder(s.cb, s.r, s.h))

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
    send_setobject(connection(b), next_id(b), meshcat_mesh(reshape(permutedims(pts),:), idxs)) #reshape(permutedims(pts), :), idxs))
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

# BIM

realize_box(b::MCAT, mat, p, dx, dy, dz) =
  let buf = buffer(b),
      bot = in_world(p),
      top = in_world(add_xyz(p, dx, dy, dz))
    write_MCAT_object(buf, "box", mat, bot, top)
    void_ref(b)
  end

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

# Layers
current_layer(b::MCAT) =
  default_MCAT_material()

current_layer(layer, b::MCAT) =
  default_MCAT_material(layer)

backend_create_layer(b::MCAT, name::String, active::Bool, color::RGB) =
  begin
    @assert active
    MCAT_material(name, red=red(color), green=green(color), blue=blue(color))
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

used_materials(b::MCAT) =
  unique(map(f -> realize(s.family, b), b.shapes))

####################################################
=#





#=




c = Cone(Point(1,2,3), Point(4,5,6), 2)
vis = Visualizer()
open(vis)

shape_counter = 0



add_object!(vis, Cone(Point(3,5,3), Point(4,5,6), 2))
add_object!(vis, GeometryTypes.Cylinder(GeometryTypes.Point(-1.,-1.,-1.), GeometryTypes.Point(2., 3., 4.), 0.1))


backend(meshcat)

test_surface(surface_grid(map_division((u,v)->xyz(u,v,sin(u+v)), 0, 6pi, 20, 0, 6pi, 20)))





using Khepri
for i in 1:1000
    setobject!(vis["cone"][string(i)],
        Cone(Point(random_range(-100, 100),random_range(-100, 100),random_range(-100, 100)),
             Point(random_range(-100, 100),random_range(-100, 100),random_range(-100, 100)),
             random_range(1, 5)))
         end

for i in 1:2000
 setobject!(vis2["cone"][string(i)],
     HyperEllipsoid(Point(random_range(-100, 100),random_range(-100, 100),random_range(-100, 100)),
                    let r = random_range(1, 5); Vec(r, r, r) end))
      end


render(vis)


vis2 = Visualizer()
render(vis2)
p4 = polyhedron(v, CDDLib.Library())

# Project that polyhedron down to 3 dimensions for visualization
v1 = [1, -1,  0,  0]
v2 = [1,  1, -2,  0]
v3 = [1,  1,  1, -3]
p3 = project(p4, [v1 v2 v3])

# Show the result
setobject!(vis, Polyhedra.Mesh(p3))



p1 = Point(3, 1)
p2 = Point(1, 3)
p3 = Point(4, 4)

l1 = Line(p1, p2)
l2 = Line(p2, p3)
ls1 = LineString([l1, l2])
ls2 = LineString([p1, p2, p3])

pl = Polygon(Point{2, Int}[(3, 1), (4, 4), (2, 4), (1, 2), (3, 1)])

rect = Rect(Vec(0.0, 0.0), Vec(1.0, 1.0))

rect_faces = decompose(TriangleFace{Int}, rect)
rect_vertices = decompose(Point{2, Float64}, rect)

mesh = Mesh(rect_vertices, rect_faces)

vis = Visualizer()
open(vis)

setobject!(vis, pl)
=#

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
