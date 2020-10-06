export meshcat,
       meshcat_material,
       meshcat_glass_material,
       meshcat_metal_material,
       meshcat_line_material


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
    #vis.core.tree[path].object = data
    write(vis.core, data)
    path
  end

send_settransform(vis, path, transform) =
  let msg = (type="set_transform", path=path, matrix=convert(Vector{Float32}, transform)),
      data = pack(msg)
    #vis.core.tree[path].transform = data
    write(vis.core, data)
    nothing
  end

send_setproperty(vis, path, property, value) =
  let msg = (type="set_property", path=path, property=property, value=value),
      data = pack(msg)
    #vis.core.tree[path].properties[property] = data
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

meshcat_metadata(type="Object") =
  (type=type, version=4.5)

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
meshcat_glass_material(opacity=0.5, color=RGB(0.9,0.9,1.0)) =
  (uuid=string(uuid1()),
   type="MeshPhysicalMaterial",
    #"vertexColors" => 0,
   transparent=true,
   opacity=opacity,
    #"depthTest" => true,
    #"linewidth" => 1.0,
    #"depthFunc" => 3,
   side=2,
   color="0x$(hex(color))",
   reflectivity=0.1,
    #depthWrite=true
  )
meshcat_metal_material(roughness=0.5, color=RGB(1.0,1.0,1.0)) =
  (uuid=string(uuid1()),
   type="MeshStandardMaterial",
   metalness=1,
   roughness=roughness,
   side=2,
   color="0x$(hex(color))",
  )

meshcat_line_material(color) =
  (uuid=string(uuid1()),
   type="LineBasicMaterial",
   color="0x$(hex(color))",
   #linewidth=2, Due to limitations of the OpenGL Core Profile with the WebGL renderer on most platforms linewidth will always be 1 regardless of the set value.
   #depthFunc=3,
   #depthTest=true,
   depthWrite=true,
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

meshcat_object_2D(type, geom, shapes, material) =
 (metadata=meshcat_metadata(),
  shapes=shapes,
  geometries=[geom],
  materials=[material],
  object=(uuid=string(uuid1()),
          type=type,
          geometry=geom.uuid,
          material=material.uuid,
          matrix=(1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0,-1, 0,
                  0, 0, 0, 1)))

meshcat_object_shapes(type, geom, shapes, material, p=u0(world_cs)) =
 (metadata=meshcat_metadata(),
  shapes=shapes,
  geometries=[geom],
  materials=[material],
  object=(uuid=string(uuid1()),
          type=type,
          geometry=geom.uuid,
          material=material.uuid,
          matrix=(translated_cs(p.cs, p.x, p.y, p.z).transform*[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1])[:])) # meshcat_transform(p))) #(1, 0, 0, 0, 0, 1, 0, 0, 0, 0,-1, 0, 0, 0, 0, 1)))


meshcat_buffer_geometry_attributes_position(vertices) =
  (itemSize=3,
   type="Float32Array",
   array=convert(Vector{Float32}, reduce(vcat, meshcat_point.(vertices))))

meshcat_line(vertices, material) =
  let geom = (uuid=string(uuid1()),
              type="BufferGeometry",
              data=(
                attributes=(
                  position=meshcat_buffer_geometry_attributes_position(vertices),),))
    meshcat_object("Line", geom, material)
  end

#=
Three.js uses 2D locations and 3D locations
=#

abstract type Meshcat2D end

convert(::Type{Meshcat2D}, p::Loc) = (cx(p),cy(p))
meshcat_2d(p::Loc) = let z =  @assert(abs(cz(p)) < 1e-10); (cx(p),cy(p)) end
meshcat_3d(p::Loc) = (cx(p),cy(p),cz(p))
meshcat_line_curve_2d(v1::Loc, v2::Loc) =
  (type="LineCurve", v1=meshcat_2d(v1), v2=meshcat_2d(v2))
meshcat_line_curve_3d(v1::Loc, v2::Loc) =
  (type="LineCurve3", v1=meshcat_3d(v1), v2=meshcat_3d(v2))

#=
Three.js provides a hierarchy of curves.
Curve - Abstract
2D curves:
  ArcCurve
  CubicBezierCurve
  EllipseCurve
  LineCurve
  QuadraticBezierCurve
  SplineCurve
3D curves:
  CatmullRomCurve3
  CubicBezierCurve3
  LineCurve3
  QuadraticBezierCurve3
Sequences:
  CurvePath - Abstract
    Path
      Shape
=#
abstract type MeshcatCurve end
abstract type MeshcatCurve2D <: MeshcatCurve end
abstract type MeshcatCurve3D <: MeshcatCurve end
abstract type MeshcatCurvePath <: MeshcatCurve end
abstract type MeshcatPath <: MeshcatCurvePath end
abstract type MeshcatShape <: MeshcatPath end

abstract type MeshcatCurves end

meshcat_curve(path) = convert(MeshcatCurve, path)
convert(::Type{MeshcatCurve}, p::CircularPath) =
  (type="EllipseCurve",
   aX=cx(p.center), aY=cy(p.center),
   xRadius=p.radius, yRadius=p.radius,
   aStartAngle=0, aEndAngle=2π,
   aClockwise=false,
   aRotation=0)
convert(::Type{MeshcatCurve}, p::OpenPolygonalPath) =
  let ps = path_vertices(p)
    length(ps) == 2 ?
      meshcat_line_curve_2d(ps[1], ps[2]) :
      convert(MeshcatPath, p)
  end

meshcat_path(path) = convert(MeshcatPath, path)
convert(::Type{MeshcatPath}, path) =
  (type="Path",
   curves=meshcat_curves(path),
   autoclose=false,
   currentPoint=(0,0))

meshcat_curves(path) = convert(MeshcatCurves, path)
convert(::Type{MeshcatCurves}, path) =
  [meshcat_curve(path)]
convert(::Type{MeshcatCurves}, vs::Locs) =
  [meshcat_line_curve_2d(v1, v2)
   for (v1,v2) in zip(vs, circshift(vs, -1))]
convert(::Type{MeshcatCurves}, p::Union{RectangularPath, ClosedPolygonalPath}) =
  convert(MeshcatCurves, path_vertices(p))

meshcat_shape(path) = convert(MeshcatShape, path)
convert(::Type{MeshcatShape}, p::PathSet) =
  (uuid=string(uuid1()),
   type="Shape",
   curves=meshcat_curves(p.paths[1]),
   autoclose=false,
   currentPoint=(0,0),
   holes=meshcat_path.(p.paths[2:end]))
convert(::Type{MeshcatShape}, p::Path) =
  (uuid=string(uuid1()),
   type="Shape",
   curves=meshcat_curves(p),
   autoclose=false,
   currentPoint=(0,0),
   holes=[])

meshcat_surface_2d(path, material, p=u0(world_cs)) =
  let shape = meshcat_shape(path),
      geom = (uuid=string(uuid1()),
              type="ShapeBufferGeometry",
              shapes=[shape.uuid],
              curveSegments=64)
    meshcat_object_shapes("Mesh", geom, [shape], material, p)
  end

#=
backend(meshcat)
delete_all_shapes()
add_object(
  meshcat,
  meshcat_surface_polygon_2d(
    path_vertices(closed_polygonal_path(regular_polygon_vertices(5, u0(), 3))),
    [path_vertices(closed_polygonal_path(regular_polygon_vertices(4, u0(), 1)))],
  material(meshcat)))

dump(meshcat_surface_polygon_2d(
  path_vertices(closed_polygonal_path(regular_polygon_vertices(5, u0(), 3))),
  [path_vertices(closed_polygonal_path(regular_polygon_vertices(4, u0(), 1)))],
  material(meshcat)))
dump(meshcat_surface_2d(
    path_set(closed_polygonal_path(regular_polygon_vertices(5, u0(), 3)),
             closed_polygonal_path(regular_polygon_vertices(4, u0(), 1))),
    material(meshcat)))

add_object(meshcat, meshcat_surface_polygon(regular_polygon_vertices(5), material(meshcat)))
add_object(
  meshcat,
  meshcat_surface_2d(
    path_set(
      circular_path(u0(), 2),
      circular_path(u0(), 1)),
    material(meshcat)))

add_object(
 meshcat,
  meshcat_surface_2d(circular_path(u0(), 3), material(meshcat)))

add_object(
  meshcat,
  meshcat_surface_2d(
    closed_polygonal_path(regular_polygon_vertices(5, u0(), 5)),
    material(meshcat)))
add_object(
  meshcat,
  meshcat_surface_2d(
    closed_polygonal_path(regular_polygon_vertices(4, u0(), 1)),
    material(meshcat)))



add_object(
  meshcat,
  meshcat_surface_2d(
    path_set(
      closed_polygonal_path(regular_polygon_vertices(5, u0(), 3)),
      closed_polygonal_path(regular_polygon_vertices(4, u0(), 1))),
    material(meshcat)))


dump(meshcat_surface_polygon(regular_polygon_vertices(5), material(meshcat)))


dump(meshcat_surface_2d(
  path_set(
    circular_path(u0(), 5),
    circular_path(u0(), 1)),
  material(meshcat)))
=#

meshcat_surface_polygon(vertices, material) =
  let n = length(vertices)
    n <= 0 ?
      meshcat_mesh(vertices, n < 4 ? [(0,1,2)] : [(0,1,2),(2,3,0)], material) :
      let ps = in_world.(vertices),
          n = vertices_normal(ps),
          cs = cs_from_o_vz(ps[1], n),
          vs = [in_cs(v, cs) for v in vertices]
        meshcat_surface_2d(closed_polygonal_path(vs), material, u0(cs))
      end
  end

#=
reset_backend()
delete_all_shapes()
add_object(meshcat, meshcat_surface_polygon([x(1), x(3), xz(3,4), xz(1,4)], material(meshcat)))
=#

meshcat_circle_mesh(center, radius, start_angle, amplitude, material) =
  let geom = (uuid=string(uuid1()),
              type="CircleBufferGeometry",
              radius=radius,
              segments=64,
              thetaStart=start_angle,
              thetaLength=amplitude),
      cs = cs_from_o_vz(center, vy(1, center.cs))
    meshcat_object("Mesh", geom, material, u0(cs))
  end

meshcat_centered_box(p, dx, dy, dz, material) =
  let geom = (uuid=string(uuid1()),
              type="BoxBufferGeometry",
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
              type="SphereBufferGeometry",
              radius=r,
              widthSegments=64,
              heightSegments=64)
    meshcat_object("Mesh", geom, material, p)
  end

meshcat_torus(p, re, ri, material) =
  let geom = (uuid=string(uuid1()),
              type="TorusBufferGeometry",
              radius=re,
              tube=ri,
              radialSegments=64,
              tubularSegments=16)
    meshcat_object_shapes("Mesh", geom, [], material, p)
  end

meshcat_centered_cone(p, r, h, material) =
  let geom = (uuid=string(uuid1()),
              type="ConeBufferGeometry",
              radius=r,
              height=h,
              radialSegments=64)
    meshcat_object("Mesh", geom, material, p)
  end

meshcat_cone(p, r, h, material) =
  meshcat_centered_cone(p+vz(h/2, p.cs), r, h, material)

meshcat_centered_cone_frustum(p, rb, rt, h, material) =
  let geom = (uuid=string(uuid1()),
              type="CylinderBufferGeometry",
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

meshcat_extrusion_z(profile, h, material, p=u0(world_cs)) =
  let shape = meshcat_shape(profile),
      geom = (uuid=string(uuid1()),
              type="ExtrudeBufferGeometry",
              shapes=[shape.uuid],
              options=(
                #steps=2,
                depth=-h,
                bevelEnabled=false,
                #bevelThickness=1,
                #bevelSize=1,
                #bevelOffset=0,
                #bevelSegments=1,
                #extrudePath=meshcat_path(path),
                curveSegments=64
                ))
    meshcat_object_shapes("Mesh", geom, [shape], material, p)
  end
#

#=
delete_all_shapes()
add_object(
  meshcat,
  meshcat_extrusion_z(
    path_set(
      circular_path(u0(), 2),
      circular_path(u0(), 1)),
    10,
    material(meshcat),
    y(10)))

delete_all_shapes()
add_object(
  meshcat,
  meshcat_extrusion(open_polygonal_path([x(0), z(10)]),
                    closed_polygonal_path(regular_polygon_vertices(5, u0(), 3)),
                    material(meshcat)))
dump(meshcat_extrusion(open_polygonal_path([x(0), z(10)]),
                  closed_polygonal_path(regular_polygon_vertices(5, u0(), 3)),
                  material(meshcat)))
=#

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
#=
meshcat_extruded_surface_polygon(outer_vertices, holes_vertices, material) =
  let ps = in_world.(outer_vertices),
      n = vertices_normal(ps),
      cs = cs_from_o_vz(ps[1], n),
      vs = [in_cs(v, cs) for v in outer_vertices],
      hsvs = [[in_cs(v, cs) for v in hvs] for hvs in holes_vertices],
      curves = [(#metadata=meshcat_metadata("Curve"),
                 type="LineCurve",
                 v1=(cx(v1),cy(v1)), v2=(cx(v2),cy(v2)))
                for (v1,v2) in zip(vs, circshift(vs, -1))],
      shape = (uuid=string(uuid1()),
               type="Shape",
               curves=curves,
               autoclose=false,
               currentPoint=(cx(vs[1]), cy(vs[1])),
               holes=[],
               ),
      geom = (uuid=string(uuid1()),
              type="ShapeBufferGeometry",
              shapes=[shape.uuid],
              #curveSegments=12,
              )
    #meshcat_object_2D("Mesh", geom, [shape], material, )
    meshcat_object_shapes("Mesh", geom, [shape], material, u0(cs))
  end
=#

letter_glyph = Dict{Char, NamedTuple}(
  ' '=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[]),
  '!'=>(bb=[xy(0, 0), xy(0, 1)], vss=[[xy(0, 0), xy(0, 1/6)], [xy(0, 1/3), xy(0, 1)]]),
  '"'=>(bb=[xy(0, 2/3), xy(1/3, 1)], vss=[[xy(0, 2/3), xy(1/6, 1)], [xy(1/3, 1), xy(1/6, 2/3)]]),
  '#'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/3), xy(2/3, 1/3)], [xy(2/3, 2/3), xy(0, 2/3)], [xy(1/6, 1), xy(1/6, 0)], [xy(1/2, 0), xy(1/2, 1)]]),
  '$'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/6), xy(1/2, 1/6), xy(2/3, 1/3), xy(1/2, 1/2), xy(1/6, 1/2), xy(0, 2/3), xy(1/6, 5/6), xy(2/3, 5/6)], [xy(1/3, 1), xy(1/3, 0)]]),
  '%'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(0, 5/6), xy(1/6, 5/6), xy(1/6, 1), xy(0, 1)], [xy(2/3, 1), xy(0, 0)], [xy(2/3, 0), xy(1/2, 0), xy(1/2, 1/6), xy(2/3, 1/6), xy(2/3, 0)]]),
  '&'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(2/3, 1/3), xy(1/3, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 1/3), xy(1/3, 2/3), xy(1/3, 5/6), xy(1/6, 1), xy(0, 5/6), xy(0, 2/3), xy(2/3, 0)]]),
  '\''=>(bb=[xy(0, 2/3), xy(1/6, 1)], vss=[[xy(0, 2/3), xy(1/6, 1)]]),
  '('=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(1/3, 1), xy(0, 2/3), xy(0, 1/3), xy(1/3, 0)]]),
  ')'=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(0, 1), xy(1/3, 2/3), xy(1/3, 1/3), xy(0, 0)]]),
  '*'=>(bb=[xy(0, 1/6), xy(2/3, 5/6)], vss=[[xy(1/3, 1/6), xy(1/3, 5/6)], [xy(2/3, 1/2), xy(0, 1/2)], [xy(2/3, 5/6), xy(0, 1/6)], [xy(0, 5/6), xy(2/3, 1/6)]]),
  '+'=>(bb=[xy(0, 1/6), xy(2/3, 5/6)], vss=[[xy(1/3, 1/6), xy(1/3, 5/6)], [xy(2/3, 1/2), xy(0, 1/2)]]),
  ','=>(bb=[xy(0, -1/6), xy(1/6, 1/6)], vss=[[xy(1/6, 1/6), xy(1/6, 0), xy(0, -1/6)]]),
  '-'=>(bb=[xy(0, 1/2), xy(2/3, 1/2)], vss=[[xy(0, 1/2), xy(2/3, 1/2)]]),
  '.'=>(bb=[xy(0, 0), xy(0, 1/6)], vss=[[xy(0, 0), xy(0, 1/6)]]),
  '/'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(2/3, 1)]]),
  '0'=>(bb=[xy(0, 0), xy(1/2, 1)], vss=[[xy(1/6, 0), xy(0, 1/6), xy(0, 5/6), xy(1/6, 1), xy(1/3, 1), xy(1/2, 5/6), xy(1/2, 1/6), xy(1/3, 0), xy(1/6, 0)]]),
  '1'=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(0, 5/6), xy(1/6, 1), xy(1/6, 0)], [xy(0, 0), xy(1/3, 0)]]),
  '2'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 5/6), xy(1/6, 1), xy(1/2, 1), xy(2/3, 5/6), xy(2/3, 2/3), xy(1/2, 1/2), xy(1/6, 1/2), xy(0, 1/3), xy(0, 0), xy(2/3, 0)]]),
  '3'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 5/6), xy(1/6, 1), xy(1/2, 1), xy(2/3, 5/6), xy(2/3, 2/3), xy(1/2, 1/2), xy(1/3, 1/2)], [xy(1/2, 1/2), xy(2/3, 1/3), xy(2/3, 1/6), xy(1/2, 0), xy(1/6, 0), xy(0, 1/6)]]),
  '4'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(2/3, 1/3), xy(0, 1/3), xy(1/2, 1), xy(1/2, 0)]]),
  '5'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/6), xy(1/6, 0), xy(1/2, 0), xy(2/3, 1/6), xy(2/3, 1/2), xy(1/2, 2/3), xy(0, 2/3), xy(0, 1), xy(2/3, 1)]]),
  '6'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/2), xy(1/2, 1/2), xy(2/3, 1/3), xy(2/3, 1/6), xy(1/2, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 2/3), xy(1/3, 1), xy(1/2, 1)]]),
  '7'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(2/3, 1), xy(1/6, 0)]]),
  '8'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(1/6, 0), xy(0, 1/6), xy(0, 1/3), xy(1/6, 1/2), xy(1/2, 1/2), xy(2/3, 2/3), xy(2/3, 5/6), xy(1/2, 1), xy(1/6, 1), xy(0, 5/6), xy(0, 2/3), xy(1/6, 1/2)], [xy(1/2, 1/2), xy(2/3, 1/3), xy(2/3, 1/6), xy(1/2, 0), xy(1/6, 0)]]),
  '9'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(1/6, 0), xy(1/3, 0), xy(2/3, 1/3), xy(2/3, 5/6), xy(1/2, 1), xy(1/6, 1), xy(0, 5/6), xy(0, 2/3), xy(1/6, 1/2), xy(2/3, 1/2)]]),
  ':'=>(bb=[xy(0, 1/6), xy(0, 2/3)], vss=[[xy(0, 2/3), xy(0, 1/2)], [xy(0, 1/3), xy(0, 1/6)]]),
  ';'=>(bb=[xy(0, -1/6), xy(1/6, 2/3)], vss=[[xy(1/6, 2/3), xy(1/6, 1/2)], [xy(1/6, 1/3), xy(1/6, 0), xy(0, -1/6)]]),
  '<'=>(bb=[xy(0, 0), xy(1/2, 1)], vss=[[xy(1/2, 1), xy(0, 1/2), xy(1/2, 0)]]),
  '='=>(bb=[xy(0, 1/3), xy(2/3, 2/3)], vss=[[xy(0, 2/3), xy(2/3, 2/3)], [xy(2/3, 1/3), xy(0, 1/3)]]),
  '>'=>(bb=[xy(0, 0), xy(1/2, 1)], vss=[[xy(0, 1), xy(1/2, 1/2), xy(0, 0)]]),
  '?'=>(bb=[xy(0, 0), xy(1/2, 1)], vss=[[xy(0, 5/6), xy(1/6, 1), xy(1/3, 1), xy(1/2, 5/6), xy(1/2, 2/3), xy(1/3, 1/2), xy(1/3, 1/3)], [xy(1/3, 1/6), xy(1/3, 0)]]),
  '@'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(1/2, 1/2), xy(1/3, 1/3), xy(1/6, 1/3), xy(1/6, 1/2), xy(1/3, 2/3), xy(1/2, 2/3), xy(1/2, 1/3), xy(2/3, 1/2), xy(2/3, 5/6), xy(1/2, 1), xy(1/6, 1), xy(0, 5/6), xy(0, 1/6), xy(1/6, 0), xy(2/3, 0)]]),
  'A'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1/3), xy(1/3, 1), xy(2/3, 1/3), xy(2/3, 0)], [xy(0, 1/3), xy(2/3, 1/3)]]),
  'B'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(1/2, 0), xy(2/3, 1/6), xy(2/3, 1/3), xy(1/2, 1/2), xy(1/6, 1/2)], [xy(1/2, 1/2), xy(2/3, 2/3), xy(2/3, 5/6), xy(1/2, 1), xy(0, 1)], [xy(1/6, 1), xy(1/6, 0)]]),
  'C'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(2/3, 1/6), xy(1/2, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 5/6), xy(1/6, 1), xy(1/2, 1), xy(2/3, 5/6)]]),
  'D'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(1/2, 0), xy(2/3, 1/6), xy(2/3, 5/6), xy(1/2, 1), xy(0, 1)], [xy(1/6, 1), xy(1/6, 0)]]),
  'E'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(2/3, 1)], [xy(0, 1/2), xy(1/3, 1/2)], [xy(0, 0), xy(2/3, 0)]]),
  'F'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(2/3, 1)], [xy(0, 1/2), xy(1/3, 1/2)]]),
  'G'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(1/2, 1/2), xy(2/3, 1/2), xy(2/3, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 5/6), xy(1/6, 1), xy(2/3, 1)]]),
  'H'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1)], [xy(0, 1/2), xy(2/3, 1/2)], [xy(2/3, 1), xy(2/3, 0)]]),
  'I'=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(0, 1), xy(1/3, 1)], [xy(1/6, 1), xy(1/6, 0)], [xy(0, 0), xy(1/3, 0)]]),
  'J'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/6), xy(1/6, 0), xy(1/2, 0), xy(2/3, 1/6), xy(2/3, 1)]]),
  'K'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1)], [xy(2/3, 1), xy(1/6, 1/2), xy(0, 1/2)], [xy(1/6, 1/2), xy(2/3, 0)]]),
  'L'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(0, 0), xy(2/3, 0)]]),
  'M'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(1/3, 1/3), xy(2/3, 1), xy(2/3, 0)]]),
  'N'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(2/3, 0), xy(2/3, 1)]]),
  'O'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(2/3, 1), xy(2/3, 0), xy(0, 0)]]),
  'P'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(1/2, 1), xy(2/3, 5/6), xy(2/3, 2/3), xy(1/2, 1/2), xy(0, 1/2)]]),
  'Q'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(1/3, 1/3), xy(1/2, 1/6), xy(1/3, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 5/6), xy(1/6, 1), xy(1/2, 1), xy(2/3, 5/6), xy(2/3, 1/3), xy(1/2, 1/6), xy(2/3, 0)]]),
  'R'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(1/2, 1), xy(2/3, 5/6), xy(2/3, 2/3), xy(1/2, 1/2), xy(0, 1/2)], [xy(1/6, 1/2), xy(2/3, 0)]]),
  'S'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/6), xy(1/6, 0), xy(1/2, 0), xy(2/3, 1/6), xy(0, 5/6), xy(1/6, 1), xy(1/2, 1), xy(2/3, 5/6)]]),
  'T'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(2/3, 1)], [xy(1/3, 1), xy(1/3, 0)]]),
  'U'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(0, 1/6), xy(1/6, 0), xy(1/2, 0), xy(2/3, 1/6), xy(2/3, 1)]]),
  'V'=>(bb=[xy(0, 0), xy(1, 1)], vss=[[xy(0, 1), xy(1/2, 0), xy(1, 1)]]),
  'W'=>(bb=[xy(0, 0), xy(1, 1)], vss=[[xy(0, 1), xy(1/3, 0), xy(1/2, 1/2), xy(2/3, 0), xy(1, 1)]]),
  'X'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(2/3, 1)], [xy(0, 1), xy(2/3, 0)]]),
  'Y'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(1/3, 1/2), xy(1/3, 0)], [xy(1/3, 1/2), xy(2/3, 1)]]),
  'Z'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(2/3, 1), xy(0, 0), xy(2/3, 0)]]),
  '['=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(0, 0), xy(0, 1), xy(1/3, 1)], [xy(1/3, 0), xy(0, 0)]]),
  '\\'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1), xy(2/3, 0)]]),
  ']'=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(0, 1), xy(1/3, 1), xy(1/3, 0), xy(0, 0)]]),
  '^'=>(bb=[xy(0, 2/3), xy(2/3, 1)], vss=[[xy(0, 2/3), xy(1/3, 1), xy(2/3, 2/3)]]),
  '_'=>(bb=[xy(0, -1/6), xy(2/3, -1/6)], vss=[[xy(0, -1/6), xy(2/3, -1/6)]]),
  '`'=>(bb=[xy(0, 2/3), xy(1/6, 1)], vss=[[xy(0, 1), xy(1/6, 2/3)]]),
  'a'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(1/3, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 1/2), xy(1/6, 2/3), xy(1/3, 2/3), xy(1/2, 1/2), xy(1/2, 1/6), xy(1/3, 0)], [xy(1/2, 1/6), xy(2/3, 0)]]),
  'b'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1)], [xy(0, 1/3), xy(1/3, 2/3), xy(1/2, 2/3), xy(2/3, 1/2), xy(2/3, 1/6), xy(1/2, 0), xy(1/3, 0), xy(0, 1/3)]]),
  'c'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(2/3, 2/3), xy(1/6, 2/3), xy(0, 1/2), xy(0, 1/6), xy(1/6, 0), xy(2/3, 0)]]),
  'd'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(2/3, 1/3), xy(1/3, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 1/2), xy(1/6, 2/3), xy(1/3, 2/3), xy(2/3, 1/3)], [xy(2/3, 1), xy(2/3, 0)]]),
  'e'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 1/3), xy(1/2, 1/3), xy(2/3, 1/2), xy(1/2, 2/3), xy(1/6, 2/3), xy(0, 1/2), xy(0, 1/6), xy(1/6, 0), xy(1/2, 0)]]),
  'f'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 1/2), xy(1/2, 1/2)], [xy(2/3, 5/6), xy(1/2, 1), xy(1/3, 1), xy(1/6, 5/6), xy(1/6, 0)]]),
  'g'=>(bb=[xy(0, -1/3), xy(2/3, 2/3)], vss=[[xy(0, -1/6), xy(1/6, -1/3), xy(1/2, -1/3), xy(2/3, -1/6), xy(2/3, 1/2), xy(1/2, 2/3), xy(1/6, 2/3), xy(0, 1/2), xy(0, 1/6), xy(1/6, 0), xy(2/3, 0)]]),
  'h'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1)], [xy(0, 1/3), xy(1/3, 2/3), xy(1/2, 2/3), xy(2/3, 1/2), xy(2/3, 0)]]),
  'i'=>(bb=[xy(0, 0), xy(0, 1)], vss=[[xy(0, 0), xy(0, 2/3)], [xy(0, 5/6), xy(0, 1)]]),
  'j'=>(bb=[xy(0, -1/3), xy(1/2, 1)], vss=[[xy(0, -1/6), xy(1/6, -1/3), xy(1/3, -1/3), xy(1/2, -1/6), xy(1/2, 2/3)], [xy(1/2, 5/6), xy(1/2, 1)]]),
  'k'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 0), xy(0, 1)], [xy(0, 1/3), xy(1/3, 1/3), xy(2/3, 2/3)], [xy(1/3, 1/3), xy(2/3, 0)]]),
  'l'=>(bb=[xy(0, 0), xy(1/6, 1)], vss=[[xy(0, 1), xy(0, 1/6), xy(1/6, 0)]]),
  'm'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 0), xy(0, 2/3)], [xy(0, 1/2), xy(1/6, 2/3), xy(1/3, 1/2), xy(1/3, 1/3)], [xy(1/3, 1/2), xy(1/2, 2/3), xy(2/3, 1/2), xy(2/3, 0)]]),
  'n'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 0), xy(0, 2/3)], [xy(0, 1/3), xy(1/3, 2/3), xy(1/2, 2/3), xy(2/3, 1/2), xy(2/3, 0)]]),
  'o'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(1/2, 0), xy(1/6, 0), xy(0, 1/6), xy(0, 1/2), xy(1/6, 2/3), xy(1/2, 2/3), xy(2/3, 1/2), xy(2/3, 1/6), xy(1/2, 0)]]),
  'p'=>(bb=[xy(0, -1/3), xy(2/3, 2/3)], vss=[[xy(0, -1/3), xy(0, 2/3)], [xy(0, 1/2), xy(1/6, 2/3), xy(1/2, 2/3), xy(2/3, 1/2), xy(2/3, 1/6), xy(1/2, 0), xy(0, 0)]]),
  'q'=>(bb=[xy(0, -1/3), xy(2/3, 2/3)], vss=[[xy(2/3, -1/3), xy(2/3, 2/3)], [xy(2/3, 1/2), xy(1/2, 2/3), xy(1/6, 2/3), xy(0, 1/2), xy(0, 1/6), xy(1/6, 0), xy(2/3, 0)]]),
  'r'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 0), xy(0, 2/3)], [xy(0, 1/3), xy(1/3, 2/3), xy(1/2, 2/3), xy(2/3, 1/2)]]),
  's'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 0), xy(1/2, 0), xy(2/3, 1/6), xy(1/2, 1/3), xy(1/6, 1/3), xy(0, 1/2), xy(1/6, 2/3), xy(2/3, 2/3)]]),
  't'=>(bb=[xy(0, 0), xy(2/3, 1)], vss=[[xy(0, 2/3), xy(2/3, 2/3)], [xy(1/3, 1), xy(1/3, 1/6), xy(1/2, 0), xy(2/3, 1/6)]]),
  'u'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 2/3), xy(0, 1/6), xy(1/6, 0), xy(1/3, 0), xy(2/3, 1/3)], [xy(2/3, 2/3), xy(2/3, 0)]]),
  'v'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 2/3), xy(1/3, 0), xy(2/3, 2/3)]]),
  'w'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 2/3), xy(1/6, 0), xy(1/3, 2/3), xy(1/2, 0), xy(2/3, 2/3)]]),
  'x'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 0), xy(2/3, 2/3)], [xy(0, 2/3), xy(2/3, 0)]]),
  'y'=>(bb=[xy(0, -1/3), xy(2/3, 2/3)], vss=[[xy(0, 2/3), xy(1/3, 0)], [xy(2/3, 2/3), xy(1/6, -1/3), xy(0, -1/3)]]),
  'z'=>(bb=[xy(0, 0), xy(2/3, 2/3)], vss=[[xy(0, 2/3), xy(2/3, 2/3), xy(0, 0), xy(2/3, 0)]]),
  '{'=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(1/3, 1), xy(1/6, 5/6), xy(1/6, 2/3), xy(0, 1/2), xy(1/6, 1/3), xy(1/6, 1/6), xy(1/3, 0)]]),
  '|'=>(bb=[xy(0, 0), xy(0, 1)], vss=[[xy(0, 0), xy(0, 1)]]),
  '}'=>(bb=[xy(0, 0), xy(1/3, 1)], vss=[[xy(0, 0), xy(1/6, 1/6), xy(1/6, 1/3), xy(1/3, 1/2), xy(1/6, 2/3), xy(1/6, 5/6), xy(0, 1)]]),
  '~'=>(bb=[xy(0, 1/2), xy(2/3, 2/3)], vss=[[xy(0, 1/2), xy(1/6, 2/3), xy(1/2, 1/2), xy(2/3, 2/3)]])
)

inter_letter_spacing_factor = 1/3

meshcat_text(str::String, p::Loc, size::Real) =
  for c in str
    let glyph = letter_glyph[c]
      for vs in glyph.vss
        line([add_xy(p, q.x*size, q.y*size) for q in vs])
      end
      p = add_x(p, (glyph.bb[2].x-glyph.bb[1].x + inter_letter_spacing_factor)*size)
    end
  end

send_setview(v, camera::Loc, target::Loc, lens::Real, aperture::Real) =
  let (x1,y1,z1) = raw_point(camera),
      (x2,y2,z2) = raw_point(target)
    send_settransform(v, "/Cameras/default", [
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      x2, y2, z2, 1])
    send_setproperty(v, "/Cameras/default/rotated/<object>", "zoom", 1)
    send_setproperty(v, "/Cameras/default/rotated/<object>", "fov", view_angles(lens)[2])
    send_setproperty(v, "/Cameras/default/rotated/<object>", "position", [x1-x2,z1-z2,y2-y1])
  end

#send_setproperty(connection(meshcat), "/Orbit/<object>", "target", [2,0,0])
####################################################
abstract type MCATKey end
const MCATId = String #MCATlyJS.GenericTrace
const MCATRef = GenericRef{MCATKey, MCATId}
const MCATNativeRef = NativeRef{MCATKey, MCATId}
const MCATUnionRef = UnionRef{MCATKey, MCATId}
const MCATSubtractionRef = SubtractionRef{MCATKey, MCATId}

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
  camera::Loc
  target::Loc
  lens::Real
  sun_altitude::Real
  sun_azimuth::Real

end

const MCAT = MCATBackend{MCATKey, MCATId}

connection(b::MCAT) = b.connection()
material(b::MCAT) = b.layer.material
line_material(b::MCAT) = b.layer.line_material

create_MCAT_connection() =
  let (width, height) = render_size(),
      vis = Visualizer()
    display(MCATViewer(vis))
    send_setproperty(vis, "/Cameras/default/rotated/<object>", "far", 50000)
    vis
  end

meshcat = MCAT(LazyParameter(Visualizer, create_MCAT_connection),
               0,
               mcat_layer("default", RGB(1,1,1)),
               xyz(10,10,10),
               xyz(0,0,0),
               35,
               90,
               0)

#=
To visualize, we piggyback on Julia's display mechanisms
=#

display_meshcat(io, vis, (w, h) = render_size()) =
  let frame = vis.core
    print(io, """
    <div style="height: $(h)px; width: $(w)px; overflow-x: auto; overflow-y: hidden; resize: both">
    <iframe src="$(MeshCat.url(frame))" style="width: 100%; height: 100%; border: none"></iframe>
    </div>
    """)
    MeshCat.wait_for_server(frame)
  end

struct MCATViewer
  visualizer
end

Base.show(
  io::IO,
  ::Union{
    MIME"text/html",
    MIME"juliavscode/html",
    MIME"application/prs.juno.plotpane+html"},
  v::MCATViewer) = display_meshcat(io, v.visualizer)

export display_view
display_view(b::Backend=current_backend()) = error("Needs specialization")
display_view(b::MCAT) = MCATViewer(connection(b))

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
  begin
    display(render(connection(b)))
    set_view(get_view()...)
  end

new_backend(b::MCAT) =
  begin
    reset(b.connection)
    display_view(b)
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
=#

set_view(camera::Loc, target::Loc, lens::Real, aperture::Real, b::MCAT) =
  let v = connection(b)
    send_setview(v, camera, target, lens, aperture)
    b.camera = camera
    b.target = target
    b.lens = lens
  end

get_view(b::MCAT) =
  b.camera, b.target, b.lens

###################################

render_view(path::String, b::MCAT) = false
  #=
  let v = connection(b)
      @remote(b, Render(
               render_width(), render_height(),
               path,
               convert_render_quality(b, render_quality()),
               convert_render_exposure(b, render_exposure())))
=#
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

realize(b::MCAT, s::Text) =
  (meshcat_text(s.str, s.corner, s.height); MCATNativeRef(""))

backend_stroke(b::MCAT, path::CircularPath) =
  add_object(b, meshcat_line(path_frames(path), line_material(b)))

backend_stroke(b::MCAT, path::RectangularPath) =
  let c = path.corner,
      dx = path.dx,
      dy = path.dy
    add_object(
      b,
      meshcat_line([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy), c],
                   line_material(b)))
  end
backend_polygon(b::MCAT, vs::Locs) =
  backend_stroke_line(b, [vs..., vs[1]])
backend_surface_polygon(b::MCAT, vs::Locs) =
  add_object(b, meshcat_surface_polygon(vs, material(b)))
backend_stroke(b::MCAT, path::ArcPath) =
  add_object(
    b,
    meshcat_line(
      path_frames(arc_path(path.center, path.radius, path.start_angle, path.amplitude)),
      line_material(b)))

backend_fill(b::MCAT, path::ClosedPolygonalPath) =
  backend_surface_polygon(b, path.vertices)

backend_fill(b::MCAT, path::RectangularPath) =
  let c = path.corner,
      dx = path.dx,
      dy = path.dy
    add_object(b, meshcat_surface_polygon([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)], material(b)))
  end

backend_stroke(b::MCAT, path::OpenSplinePath) =
  if (path.v0 == false) && (path.v1 == false)
    add_object(b, meshcat_line(path_frames(path), line_material(b)))
  elseif (path.v0 != false) && (path.v1 != false)
    @remote(b, InterpSpline(path.vertices, path.v0, path.v1))
  else
    @remote(b, InterpSpline(
                     path.vertices,
                     path.v0 == false ? path.vertices[2]-path.vertices[1] : path.v0,
                     path.v1 == false ? path.vertices[end-1]-path.vertices[end] : path.v1))
  end
backend_stroke(b::MCAT, path::ClosedSplinePath) =
    add_object(b, meshcat_line(path_frames(path), line_material(b)))
backend_fill(b::MCAT, path::ClosedSplinePath) =
    add_object(b, meshcat_surface_polygon(path_frames(path), material(b)))

#=
backend_fill_curves(b::MCAT, refs::MCATIds) = @remote(b, SurfaceFromCurves(refs))
backend_fill_curves(b::MCAT, ref::MCATId) = @remote(b, SurfaceFromCurves([ref]))
=#
backend_stroke_line(b::MCAT, vs) =
	add_object(b, meshcat_line(vs, line_material(b)))
#=
realize(b::MCAT, s::Point) =
  @remote(b, Point(s.position))
  =#
realize(b::MCAT, s::Line) =
  backend_stroke_line(b, s.vertices)

realize(b::MCAT, s::Spline) =
 # This should be merged with opensplinepath
  if (s.v0 == false) && (s.v1 == false)
    backend_stroke(b, open_spline_path(s.points))
  elseif (s.v0 != false) && (s.v1 != false)
    error("Finish this")
  else
    error("Finish this")
  end

realize(b::MCAT, s::ClosedSpline) =
  backend_stroke(b, closed_spline_path(s.points))

realize(b::MCAT, s::Circle) =
  backend_stroke(b, circular_path(s.center, s.radius))
realize(b::MCAT, s::Arc) =
  backend_stroke(b, arc_path(s.center, s.radius, s.start_angle, s.amplitude))
#=
realize(b::MCAT, s::Ellipse) =
  if s.radius_x > s.radius_y
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x))
  else
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y))
  end
realize(b::MCAT, s::EllipticArc) =
  error("Finish this")
=#
realize(b::MCAT, s::SurfaceCircle) =
  add_object(b, meshcat_circle_mesh(s.center, s.radius, 0, 2pi, material(b)))

realize(b::MCAT, s::Surface) =
  let #ids = map(r->@remote(b, NurbSurfaceFrom(r)), @remote(b, SurfaceFromCurves(collect_ref(s.frontier))))
      ids = @remote(b, SurfaceFromCurves(collect_ref(s.frontier)))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::MCAT, s::Shape2D) =
    map(c -> shape_from_ref(c, b), @remote(b, CurvesFromSurface(ref(s).value)))

realize(b::MCAT, s::Sphere) =
  add_object(b, meshcat_sphere(s.center, s.radius, material(b)))

realize(b::MCAT, s::Torus) =
  add_object(b, meshcat_torus(s.center, s.re, s.ri, material(b)))

backend_right_cuboid(b::MCAT, cb, width, height, h, material) =
  add_object(b, meshcat_box(add_xy(cb, -width/2, -height/2), width, height, h, material))

realize(b::MCAT, s::Box) =
  add_object(b, meshcat_box(s.c, s.dx, s.dy, s.dz, material(b)))

realize(b::MCAT, s::Cone) =
  add_object(b, meshcat_cone(s.cb, s.r, s.h, material(b)))

realize(b::MCAT, s::ConeFrustum) =
  add_object(b, meshcat_cone_frustum(s.cb, s.rb, s.rt, s.h, material(b)))

backend_cylinder(b::MCAT, cb::Loc, r::Real, h::Real, material=material(b)) =
  send_setobject(connection(b), next_id(b), meshcat_cylinder(cb, r, h, material))


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
  backend_surface_grid(b, s.points, s.closed_u, s.closed_v, s.smooth_u, s.smooth_v)

backend_surface_grid(b::MCAT, pts, closed_u, closed_v, smooth_u, smooth_v) =
  let si = size(pts, 1),
      sj = size(pts, 2),
      idxs = meshcat_faces(si, sj, closed_u, closed_v)
    add_object(b, meshcat_mesh(reshape(permutedims(pts),:), idxs, material(b)))
  end

smooth_pts(pts) = in_world.(path_frames(open_spline_path(pts)))

backend_surface_mesh(b::MCAT, vertices, faces) =
  add_object(b, meshcat_mesh(vertices, faces, material(b)))

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

realize_prism(b::MCAT, top, bot, side, path::PathSet, h::Real) =
  let p = path_start(path.paths[1]),
      n = planar_path_normal(path.paths[1]),
      cs = cs_from_o_vz(p, n)
    add_object(b, meshcat_extrusion_z(in_cs(path, cs), h, top, u0(cs)))
  end
#=
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
=#

realize_pyramid_frustum(b::MCAT, top, bot, side, bot_vs::Locs, top_vs::Locs, closed=true) =
  let
    if closed
      add_object(b, meshcat_surface_polygon(reverse(bot_vs), bot))
      add_object(b, meshcat_surface_polygon(top_vs, top))
    end
    add_object(b,
      meshcat_mesh([top_vs...,bot_vs...],
                   meshcat_faces(2, length(bot_vs), false, false),
                   material(b)))
#    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
#      add_object(b, meshcat_surface_polygon(vs, side))
#    end
    void_ref(b)
  end

#meshcat_faces(2, 3, false, false)

#regular_pyramid_frustum(5)


backend_surface_polygon(b::MCAT, mat, path::PathSet, acw=true) =
  backend_surface_polygon(b, mat, path_vertices(path.paths[1])) # WARNING: Ignore the remaining paths that should be subtracted

backend_surface_polygon(b::MCAT, mat, vs::Locs, acw=true) =
  add_object(b, meshcat_surface_polygon(vs, mat))
#=
# Polygons with holes need a PathSets in MCAT

subtract_paths(b::MCAT, c_r_w_path::PathSet, c_l_w_path::PathSet, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path.paths..., c_r_op_path),
  path_set(c_l_w_path.paths..., c_l_op_path)

subtract_paths(b::MCAT, c_r_w_path::Path, c_l_w_path::Path, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path, c_r_op_path),
  path_set(c_l_w_path, c_l_op_path)
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
MeshcatMaterial = NamedTuple

meshcat_material_family(mat::MeshcatMaterial) =
  BackendMaterialFamily(mat)

meshcat_slab_family(top::MeshcatMaterial, bot::MeshcatMaterial=top, side::MeshcatMaterial=bot) =
  BackendSlabFamily(top, bot, side)

meshcat_roof_family(top::MeshcatMaterial, bot::MeshcatMaterial=top, side::MeshcatMaterial=bot) =
  BackendRoofFamily(top, bot, side)

meshcat_wall_family(right::MeshcatMaterial, left::MeshcatMaterial=right) =
  BackendWallFamily(right, left)

export meshcat_material_family,
       meshcat_slab_family,
       meshcat_roof_family,
       meshcat_wall_family,
       default_meshcat_material

meshcat_generic_floor_20 = meshcat_material(RGB(0.2,0.2,0.2))
meshcat_generic_ceiling_80 = meshcat_material(RGB(0.8,0.8,0.8))
meshcat_generic_interior_wall_70 = meshcat_material(RGB(0.7,0.7,0.7))
meshcat_outside_facade_30 = meshcat_material(RGB(0.3,0.3,0.3))
meshcat_generic_metal = meshcat_metal_material()
meshcat_generic_furniture_50 = meshcat_material(RGB(0.5,0.5,0.5))
meshcat_glass = meshcat_glass_material()

set_backend_family(default_wall_family(), meshcat, meshcat_wall_family(meshcat_generic_interior_wall_70))
#set_backend_family(default_slab_family(), meshcat, meshcat_slab_family(meshcat_generic_floor_20, meshcat_generic_ceiling_80))
set_backend_family(default_slab_family(), meshcat, meshcat_material_family(meshcat_generic_ceiling_80))
set_backend_family(default_roof_family(), meshcat, meshcat_roof_family(meshcat_generic_floor_20, meshcat_outside_facade_30))
set_backend_family(default_beam_family(), meshcat, meshcat_material_family(meshcat_generic_metal))
set_backend_family(default_column_family(), meshcat, meshcat_material_family(meshcat_generic_metal))
set_backend_family(default_door_family(), meshcat, meshcat_material_family(meshcat_generic_furniture_50))
set_backend_family(default_panel_family(), meshcat, meshcat_material_family(meshcat_glass))

set_backend_family(default_truss_node_family(), meshcat, meshcat_material_family(meshcat_generic_metal))
set_backend_family(default_truss_bar_family(), meshcat, meshcat_material_family(meshcat_generic_metal))

#set_backend_family(default_table_family(), meshcat, meshcat_resource_family("Default/Prefabs/Table"))
#set_backend_family(default_chair_family(), meshcat, meshcat_resource_family("Default/Prefabs/Chair"))
#set_backend_family(default_table_chair_family(), meshcat, meshcat_resource_family("Default/Prefabs/TableChair"))

set_backend_family(default_curtain_wall_family().panel, meshcat, meshcat_material_family(meshcat_glass))
set_backend_family(default_curtain_wall_family().boundary_frame, meshcat, meshcat_material_family(meshcat_generic_metal))
set_backend_family(default_curtain_wall_family().transom_frame, meshcat, meshcat_material_family(meshcat_generic_metal))
set_backend_family(default_curtain_wall_family().mullion_frame, meshcat, meshcat_material_family(meshcat_generic_metal))

realize_slab(b::MCAT, contour::ClosedPath, holes::Vector{<:ClosedPath}, level::Level, family::Family) =
  let base = vz(level.height + slab_family_elevation(b, family)),
      thickness = slab_family_thickness(b, family),
      (mattop, matbot, matside) = slab_materials(b, family_ref(b, family))
    add_object(b, meshcat_extrusion_z(path_set(contour, holes...), thickness, mattop, base))
  end
