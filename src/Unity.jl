export unity, fast_unity,
       unity_material_family

# We need some additional Encoders
encode_GameObject = encode_int
decode_GameObject = decode_int_or_error
encode_GameObject_array = encode_int_array
decode_GameObject_array = decode_int_array
encode_Material = encode_int
decode_Material = decode_int_or_error

# Must convert from local to world coordinates
encode_Vector3(c::IO, v::Union{XYZ, VXYZ}) = begin
  v = in_world(v)
#  encode_float3(c, v.x, v.y, v.z)
  encode_float3(c, v.x, v.z, v.y)
end
decode_Vector3(c::IO) =
  let x = decode_float(c)
      z = decode_float(c)
      y = decode_float(c)
    xyz(x, y, z, world_cs)
  end
encode_Vector3_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_Vector3(c, e) end
end
decode_Vector3_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{XYZ}(undef, len)
  for i in 1:len
    r[i] = decode_Vector3(c)
  end
  r
end

encode_Vector3_array_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_Vector3_array(c, e) end
end
decode_Vector3_array_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{Vector{XYZ}}(undef, len)
  for i in 1:len
    r[i] = decode_Vector3_array(c)
  end
  r
end

encode_Color(c::IO, v::RGB) = begin
  encode_float3(c, v.r, v.g, v.b)
end
decode_Color(c::IO) =
  let r = decode_float(c)
      g = decode_float(c)
      b = decode_float(c)
    rgb(r, g, b)
  end

#=
encode_Quaternion(c::IO, pv::Union{XYZ, VXYZ}) =
  let tr = pv.cs.transform
    for i in 1:3
      for j in 1:3
        encode_float(c, tr[i,j])
      end
    end
  end
=#

unity_api = @remote_functions :CS """
public void SetApplyMaterials(bool apply)
public void SetApplyColliders(bool apply)
public GameObject SurfacePolygon(Vector3[] ps)
public GameObject Text(string txt, Vector3 position, Vector3 vx, Vector3 vy, string fontName, int fontSize)
public GameObject Sphere(Vector3 center, float radius)
public GameObject Pyramid(Vector3[] ps, Vector3 q)
public GameObject PyramidFrustum(Vector3[] ps, Vector3[] qs)
public GameObject PyramidFrustumWithMaterial(Vector3[] ps, Vector3[] qs, Material material)
public GameObject RightCuboid(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz, float angle)
public GameObject RightCuboidWithMaterial(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz, float angle, Material material)
public GameObject Box(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz)
public GameObject Cylinder(Vector3 bottom, float radius, Vector3 top)
public GameObject Unite(GameObject s0, GameObject s1)
public GameObject Intersect(GameObject s0, GameObject s1)
public GameObject Subtract(GameObject s0, GameObject s1)
public void SubtractFrom(GameObject s0, GameObject s1)
public GameObject Canonicalize(GameObject s)
public void Move(GameObject s, Vector3 v)
public void Scale(GameObject s, Vector3 p, float scale)
public void Rotate(GameObject s, Vector3 p, Vector3 n, float a)
public GameObject SurfaceFromGrid(int m, int n, Vector3[] pts, bool closedM, bool closedN, int level)
public GameObject LoadResource(String name)
public Material LoadMaterial(String name)
public Material CurrentMaterial()
public void SetCurrentMaterial(Material material)
public GameObject InstantiateResource(GameObject family, Vector3 pos, Vector3 vx, Vector3 vy, float scale)
public GameObject InstantiateBIMElement(GameObject family, Vector3 pos, float angle)
public GameObject Window(Vector3 position, Quaternion rotation, float dx, float dy, float dz)
public GameObject Shelf(Vector3 position, int rowLength, int lineLength, float cellWidth, float cellHeight, float cellDepth)
public GameObject Slab(Vector3[] contour, Vector3[][] holes, float h, Material material)
public GameObject BeamRectSection(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz, float angle, Material material)
public GameObject BeamCircSection(Vector3 bot, float radius, Vector3 top, Material material)
public GameObject Panel(Vector3[] pts, Vector3 n, Material material)
public void SetView(Vector3 position, Vector3 target, float lens)
public Vector3 ViewCamera()
public Vector3 ViewTarget()
public float ViewLens()
public void DeleteAll()
public void DeleteMany(GameObject[] objs)
public GameObject CreateParent(String name)
public GameObject CurrentParent()
public GameObject SetCurrentParent(GameObject newParent)
public void SetActive(GameObject obj, bool state)
public void DeleteAllInParent(GameObject parent)
public void SwitchToParent(GameObject newParent)
public int SetMaxNonInteractiveRequests(int n)
public void SetNonInteractiveRequests()
public void SetInteractiveRequests()
public GameObject CreateBlockInstance(GameObject block, Vector3 position, Vector3 vx, Vector3 vy, float scale)
public GameObject CreateBlockFromShapes(String name, GameObject[] objs)
public GameObject PointLight(Vector3 position, Color color, float range, float intensity)
public Point3d[] GetPosition(string prompt)
public ObjectId[] GetPoint(string prompt)
public ObjectId[] GetCurve(string prompt)
public ObjectId[] GetSurface(string prompt)
public ObjectId[] GetSolid(string prompt)
public ObjectId[] GetShape(string prompt)
public long GetHandleFromShape(Entity e)
public ObjectId GetShapeFromHandle(long h)
public void RegisterForChanges(ObjectId id)
public void UnregisterForChanges(ObjectId id)
public ObjectId[] ChangedShape()
public void DetectCancel()
public void UndetectCancel()
public bool WasCanceled()
public ObjectId[] GetAllShapes()
public ObjectId[] GetAllShapesInLayer(ObjectId layerId)
public void SetResolution(int width, int height)
public void ScreenShot(String path)
public void SelectGameObjects(GameObject[] objs)
public void StartSelectingGameObject()
public int SelectedGameObjectId(bool existing)
public void SetSun(float x, float y, float z)
public Vector3 GetSunRotation()
public string GetRenderResolution()
public float GetCurrentFPS()
public int GetViewTriangleCount()
public int GetViewVertexCount()
"""

abstract type UnityKey end
const UnityId = Int
const UnityIds = Vector{UnityId}
const UnityRef = GenericRef{UnityKey, UnityId}
const UnityRefs = Vector{UnityRef}
const UnityEmptyRef = EmptyRef{UnityKey, UnityId}
const UnityUniversalRef = UniversalRef{UnityKey, UnityId}
const UnityNativeRef = NativeRef{UnityKey, UnityId}
const UnityUnionRef = UnionRef{UnityKey, UnityId}
const UnitySubtractionRef = SubtractionRef{UnityKey, UnityId}
const Unity = SocketBackend{UnityKey, UnityId}

void_ref(b::Unity) = UnityNativeRef(-1)

create_Unity_connection() =
    begin
        #check_plugin()
        create_backend_connection("Unity", 11002)
    end

const unity = Unity(LazyParameter(TCPSocket, create_Unity_connection), unity_api)

backend_name(b::Unity) = "Unity"

realize(b::Unity, s::EmptyShape) =
  UnityEmptyRef()
realize(b::Unity, s::UniversalShape) =
  UnityUniversalRef()

fast_unity() =
  begin
    @remote(b, SetApplyMaterials(false))
    @remote(b, SetApplyColliders(false))
  end

#=
realize(b::Unity, s::Point) =
  @remote(b, Point(s.position))
realize(b::Unity, s::Line) =
  @remote(b, PolyLine(s.vertices))
realize(b::Unity, s::Spline) =
  if (s.v0 == false) && (s.v1 == false)
    #@remote(b, Spline(s.points))
    @remote(b, InterpSpline(connection(b),
                     s.points,
                     s.points[2]-s.points[1],
                     s.points[end]-s.points[end-1]))
  elseif (s.v0 != false) && (s.v1 != false)
    @remote(b, InterpSpline(s.points, s.v0, s.v1))
  else
    @remote(b, InterpSpline(connection(b),
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1))
  end
realize(b::Unity, s::ClosedSpline) =
  @remote(b, InterpClosedSpline(s.points))
realize(b::Unity, s::Circle) =
  @remote(b, Circle(s.center, vz(1, s.center.cs), s.radius))
realize(b::Unity, s::Arc) =
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

realize(b::Unity, s::Ellipse) =
  if s.radius_x > s.radius_y
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x))
  else
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y))
  end
realize(b::Unity, s::EllipticArc) =
  error("Finish this")

realize(b::Unity, s::Polygon) =
  @remote(b, ClosedPolyLine(s.vertices))
realize(b::Unity, s::RegularPolygon) =
  @remote(b, ClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))
realize(b::Unity, s::Rectangle) =
  @remote(b, ClosedPolyLine(
    connection(b),
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)]))
realize(b::Unity, s::SurfaceCircle) =
  @remote(b, SurfaceCircle(s.center, vz(1, s.center.cs), s.radius))
realize(b::Unity, s::SurfaceArc) =
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
            @remote(b, SurfaceFromCurves(connection(b),
                [@remote(b, Arc(s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle)),
                 @remote(b, PolyLine([add_pol(s.center, s.radius, end_angle),
                                              add_pol(s.center, s.radius, s.start_angle)]))])
        else
            @remote(b, SurfaceFromCurves(connection(b),
                [@remote(b, Arc(s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle),
                 @remote(b, PolyLine([add_pol(s.center, s.radius, s.start_angle)),
                                              add_pol(s.center, s.radius, end_angle)]))])
        end
    end

#realize(b::Unity, s::SurfaceElliptic_Arc) = @remote(b, Circle(connection(b),
#realize(b::Unity, s::SurfaceEllipse) = @remote(b, Circle(connection(b),
=#

realize(b::Unity, s::SurfacePolygon) =
  @remote(b, SurfacePolygon(reverse(s.vertices)))

backend_fill(b::Unity, path::ClosedPolygonalPath) =
  @remote(b, SurfacePolygon(path.vertices))

#=
realize(b::Unity, s::SurfaceRegularPolygon) =
  @remote(b, SurfaceClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))
realize(b::Unity, s::SurfaceRectangle) =
  @remote(b, SurfaceClosedPolyLine(
    connection(b),
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)])
realize(b::Unity, s::Surface) =
  let #ids = map(r->@remote(b, NurbSurfaceFrom(connection(b),r), @remote(b, SurfaceFromCurves(collect_ref(s.frontier))))
      ids = @remote(b, SurfaceFromCurves(collect_ref(s.frontier)))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::Unity, s::Shape2D) =
    map(shape_from_ref, @remote(b, CurvesFromSurface(ref(s).value)))
=#
backend_fill(b::Unity, path::ClosedPathSequence) =
  backend_fill(b, convert(ClosedPolygonalPath, path))

#=
# Iterating over curves and surfaces

Unity"public double[] CurveDomain(Entity ent)"
Unity"public double CurveLength(Entity ent)"
Unity"public Frame3d CurveFrameAt(Entity ent, double t)"
Unity"public Frame3d CurveFrameAtLength(Entity ent, double l)"
=#

backend_map_division(b::Unity, f::Function, s::Shape1D, n::Int) =
  let (t1, t2) = curve_domain(s)
    map_division(t1, t2, n) do t
      f(frame_at(s, t))
    end
  end
#=
Unity"public Vector3d RegionNormal(Entity ent)"
Unity"public Point3d RegionCentroid(Entity ent)"
Unity"public double[] SurfaceDomain(Entity ent)"
Unity"public Frame3d SurfaceFrameAt(Entity ent, double u, double v)"

backend_surface_domain(b::Unity, s::Shape2D) =
    tuple(@remote(b, SurfaceDomain(ref(s).value)...))

backend_map_division(b::Unity, f::Function, s::Shape2D, nu::Int, nv::Int) =
    let conn = connection(b)
        r = ref(s).value
        (u1, u2, v1, v2) = @remote(b, SurfaceDomain(r))
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(@remote(b, SurfaceFrameAt(r, u, v)))
            end
        end
    end

# The previous method cannot be applied to meshes in AutoCAD, which are created by surface_grid

backend_map_division(b::Unity, f::Function, s::SurfaceGrid, nu::Int, nv::Int) =
    let (u1, u2, v1, v2) = @remote(b, SurfaceDomain(r))
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(@remote(b, SurfaceFrameAt(r, u, v)))
            end
        end
    end
=#



realize(b::Unity, s::Text) =
  @remote(b, Text(s.str, s.corner, vz(-1, s.corner.cs), vy(1, s.corner.cs), "Fonts/Inconsolata-Regular", s.height))

realize(b::Unity, s::Sphere) =
  @remote(b, Sphere(s.center, s.radius))

#=
realize(b::Unity, s::Torus) =
  @remote(b, Torus(s.center, vz(1, s.center.cs), s.re, s.ri))
=#

realize(b::Unity, s::Cuboid) =
  @remote(b, PyramidFrustum([s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3]))

realize(b::Unity, s::RegularPyramidFrustum) =
    @remote(b, PyramidFrustum(
                    regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                    regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed)))

realize(b::Unity, s::RegularPyramid) =
  @remote(b, Pyramid(
               regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
               add_z(s.cb, s.h)))

realize(b::Unity, s::IrregularPyramid) =
  @remote(b, Pyramid(s.bs, s.t))

realize(b::Unity, s::RegularPrism) =
  let bs = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
    @remote(b, PyramidFrustum(bs, map(p -> add_z(p, s.h), bs)))
  end

realize(b::Unity, s::IrregularPyramidFrustum) =
    @remote(b, PyramidFrustum(s.bs, s.ts))

realize(b::Unity, s::IrregularPrism) =
  @remote(b, PyramidFrustum(s.bs, map(p -> (p + s.v), s.bs)))

realize(b::Unity, s::RightCuboid) =
  @remote(b, RightCuboid(s.cb, vz(1, s.cb.cs), vx(1, s.cb.cs), s.height, s.width, s.h, s.angle))


#unity"public GameObject Box2(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz)"

realize(b::Unity, s::Box) =
  @remote(b, Box(s.c, vz(1, s.c.cs), vx(1, s.c.cs), s.dy, s.dx, s.dz))

realize(b::Unity, s::Cone) =
  @remote(b, Pyramid(regular_polygon_vertices(64, s.cb, s.r), add_z(s.cb, s.h)))

realize(b::Unity, s::ConeFrustum) =
  @remote(b, PyramidFrustum(
    regular_polygon_vertices(64, s.cb, s.rb),
    regular_polygon_vertices(64, s.cb + vz(s.h, s.cb.cs), s.rt)))

realize(b::Unity, s::Cylinder) =
  @remote(b, Cylinder(s.cb, s.r, s.cb + vz(s.h, s.cb.cs)))

#=
backend_extrusion(b::Unity, s::Shape, v::Vec) =
    and_mark_deleted(
        map_ref(s) do r
            @remote(b, Extrude(r, v))
        end,
        s)

backend_sweep(b::Unity, path::Shape, profile::Shape, rotation::Real, scale::Real) =
  map_ref(profile) do profile_r
    map_ref(path) do path_r
      @remote(b, Sweep(path_r, profile_r, rotation, scale))
    end
  end

realize(b::Unity, s::Revolve) =
  and_delete_shape(
    map_ref(s.profile) do r
      @remote(b, Revolve(r, s.p, s.n, s.start_angle, s.amplitude))
    end,
    s.profile)

backend_loft_curves(b::Unity, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  and_delete_shapes(UnityLoft(connection(b),
                             collect_ref(profiles),
                             collect_ref(rails),
                             ruled, closed),
                    vcat(profiles, rails))

            MAYBE USE THIS
            ruled_surface(s1, s2) =
                let pts1 = map_division(in_world, s1, 20),
                    pts2 = map_division(in_world, s2, 20)
                  iterate_quads((p0, p1, p2, p3)->(surface_polygon([p0,p1,p3]), surface_polygon([p1,p2,p3])),
                                [pts1, pts2])
                end

            ruled_surface(s1, s2)


backend_loft_surfaces(b::Unity, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
    backend_loft_curves(b, profiles, rails, ruled, closed)

backend_loft_curve_point(b::Unity, profile::Shape, point::Shape) =
    and_delete_shapes(UnityLoft(connection(b),
                               vcat(collect_ref(profile), collect_ref(point)),
                               [],
                               true, false),
                      [profile, point])

backend_loft_surface_point(b::Unity, profile::Shape, point::Shape) =
    backend_loft_curve_point(b, profile, point)

=#






###
unite_ref(b::Unity, r0::UnityNativeRef, r1::UnityNativeRef) =
    ensure_ref(b, @remote(b, Unite(r0.value, r1.value)))

intersect_ref(b::Unity, r0::UnityNativeRef, r1::UnityNativeRef) =
    ensure_ref(b, @remote(b, Intersect(r0.value, r1.value)))

subtract_ref(b::Unity, r0::UnityNativeRef, r1::UnityNativeRef) =
    let r = @remote(b, Subtract(r0.value, r1.value))
      @remote(b, DeleteMany([r0.value, r1.value]))
      r
    end

#=
subtract_ref(b::Unity, r0::UnityNativeRef, r1::UnityNativeRef) =
    begin
      @remote(b, SubtractFrom(r0.value, r1.value))
      r0.value
    end
=#

#=
slice_ref(b::Unity, r::UnityNativeRef, p::Loc, v::Vec) =
    (@remote(b, Slice(r.value, p, v); r))

slice_ref(b::Unity, r::UnityUnionRef, p::Loc, v::Vec) =
    map(r->slice_ref(b, r, p, v), r.values)

=#
unite_refs(b::Unity, refs::Vector{<:UnityRef}) =
    UnityUnionRef(tuple(refs...))

#
realize(b::Unity, s::UnionShape) =
  let r = foldl((r0,r1)->unite_ref(b,r0,r1), map(ref, s.shapes),
                init=UnityEmptyRef())
    delete_shapes(s.shapes)
    #@remote(b, Canonicalize(r.value))
    r
  end

realize(b::Unity, s::IntersectionShape) =
  let r = foldl((r0,r1)->intersect_ref(b,r0,r1), map(ref, s.shapes),
                init=UnityUniversalRef())
    delete_shapes(s.shapes)
    r
  end

realize(b::Unity, s::Slice) =
  slice_ref(b, ref(s.shape), s.p, s.n)





realize(b::Unity, s::Move) =
  let r = map_ref(s.shape) do r
            @remote(b, Move(r, s.v))
            r
          end
    mark_deleted(s.shape)
    r
  end
#=
realize(b::Unity, s::Transform) =
  let r = map_ref(s.shape) do r
            @remote(b, Transform(r, s.xform))
            r
          end
    mark_deleted(s.shape)
    r
  end
=#
realize(b::Unity, s::Scale) =
  let r = map_ref(s.shape) do r
            @remote(b, Scale(r, s.p, s.s))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::Unity, s::Rotate) =
  let r = map_ref(s.shape) do r
            @remote(b, Rotate(r, s.p, s.v, s.angle))
            r
          end
    mark_deleted(s.shape)
    r
  end

#=
realize(b::Unity, s::Mirror) =
  and_delete_shape(map_ref(s.shape) do r
                    @remote(b, Mirror(r, s.p, s.n, false))
                   end,
                   s.shape)

realize(b::Unity, s::UnionMirror) =
  let r0 = ref(s.shape),
      r1 = map_ref(s.shape) do r
            @remote(b, Mirror(r, s.p, s.n, true))
          end
    UnionRef((r0,r1))
  end
=#

realize(b::Unity, s::SurfaceGrid) =
    @remote(b, SurfaceFromGrid(
        size(s.points,1),
        size(s.points,2),
        reshape(s.points,:),
        s.closed_u,
        s.closed_v,
        2))
#=
realize(b::Unity, s::Thicken) =
  and_delete_shape(
    map_ref(s.shape) do r
      @remote(b, Thicken(r, s.thickness))
    end,
    s.shape)

# backend_frame_at
backend_frame_at(b::Unity, s::Circle, t::Real) = add_pol(s.center, s.radius, t)

backend_frame_at(b::Unity, c::Shape1D, t::Real) = @remote(b, CurveFrameAt(ref(c).value, t))

#backend_frame_at(b::Unity, s::Surface, u::Real, v::Real) =
    #What should we do with v?
#    backend_frame_at(b, s.frontier[1], u)

#backend_frame_at(b::Unity, s::SurfacePolygon, u::Real, v::Real) =

backend_frame_at(b::Unity, s::Shape2D, u::Real, v::Real) = @remote(b, SurfaceFrameAt(ref(s).value, u, v))

=#

# BIM








# Families

abstract type UnityFamily <: Family end

struct UnityMaterialFamily <: UnityFamily
  name::String
  parameter_map::Dict{Symbol,String}
  ref::Parameter{Any}
end

unity_material_family(name, pairs...) = UnityMaterialFamily(name, Dict(pairs...), Parameter{Any}(nothing))
backend_get_family_ref(b::Unity, f::Family, uf::UnityMaterialFamily) = @remote(b, LoadMaterial(uf.name))

struct UnityResourceFamily <: UnityFamily
  name::String
  parameter_map::Dict{Symbol,String}
  ref::Parameter{Any}
end

unity_resource_family(name, pairs...) = UnityResourceFamily(name, Dict(pairs...), Parameter{Any}(nothing))
backend_get_family_ref(b::Unity, f::Family, uf::UnityResourceFamily) = @remote(b, LoadResource(uf.name))

set_backend_family(default_wall_family(), unity, unity_material_family("Default/Materials/Plaster"))
set_backend_family(default_slab_family(), unity, unity_material_family("Default/Materials/Concrete"))
set_backend_family(default_roof_family(), unity, unity_material_family("Default/Materials/Concrete"))
set_backend_family(default_beam_family(), unity, unity_material_family("Default/Materials/Aluminum"))
set_backend_family(default_column_family(), unity, unity_material_family("Default/Materials/Concrete"))
set_backend_family(default_door_family(), unity, unity_material_family("Default/Materials/Wood"))
set_backend_family(default_panel_family(), unity, unity_material_family("Default/Materials/Glass"))

set_backend_family(default_table_family(), unity, unity_resource_family("Default/Prefabs/Table"))
set_backend_family(default_chair_family(), unity, unity_resource_family("Default/Prefabs/Chair"))
set_backend_family(default_table_chair_family(), unity, unity_resource_family("Default/Prefabs/TableChair"))

set_backend_family(default_curtain_wall_family().panel, unity, unity_material_family("Default/Materials/Glass"))
set_backend_family(default_curtain_wall_family().boundary_frame, unity, unity_material_family("Default/Materials/Steel"))
set_backend_family(default_curtain_wall_family().transom_frame, unity, unity_material_family("Default/Materials/Steel"))
set_backend_family(default_curtain_wall_family().mullion_frame, unity, unity_material_family("Default/Materials/Steel"))

backend_rectangular_table(b::Unity, c, angle, family) =
    @remote(b, InstantiateBIMElement(realize(b, family), c, -angle))

backend_chair(b::Unity, c, angle, family) =
    @remote(b, InstantiateBIMElement(realize(b, family), c, -angle))

backend_rectangular_table_and_chairs(b::Khepri.Unity, c, angle, family) =
    @remote(b, InstantiateBIMElement(realize(b, family), c, -angle))


backend_slab(b::Unity, profile, holes, thickness, family) =
  let bot_vs = path_vertices(profile)
    @remote(b, Slab(bot_vs, map(path_vertices, holes), thickness, realize(b, family)))
  end

# A poor's man approach to deal with Z-fighting
const support_z_fighting_factor = 0.999
const wall_z_fighting_factor = 0.998

realize_beam_profile(b::Unity, s::Union{Beam,FreeColumn,Column}, profile::CircularPath, cb::Loc, length::Real) =
  @remote(b, BeamCircSection(
    cb,
    profile.radius,
    add_z(cb, length*support_z_fighting_factor), #We reduce height just a bit to avoid Z-fighting
    realize(b, s.family)))

realize_beam_profile(b::Unity, s::Union{Beam,FreeColumn,Column}, profile::RectangularPath, cb::Loc, length::Real) =
  let profile_u0 = profile.corner,
      c = add_xy(cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
    @remote(b, BeamRectSection(
      c, vz(1, c.cs), vx(1, c.cs), profile.dy, profile.dx,
      length*support_z_fighting_factor,
      -s.angle,
      realize(b, s.family)))
  end

realize(b::Unity, s::Panel) =
  let #p1 = s.vertices[1],
      #p2 = s.vertices[2],
      #p3 = s.vertices[3],
      #n = vz(s.family.thickness, cs_from_o_vx_vy(p1, p2-p1, p3-p1))
      verts = in_world.(s.vertices)
      n = vertices_normal(verts)*(s.family.thickness/2)
    @remote(b, Panel(
      map(p -> p - n, verts),
      n*2,
      realize(b, s.family)))
  end

#=
sweep_fractions(b, verts, height, l_thickness, r_thickness) =
  let p = add_z(verts[1], height/2),
      q = add_z(verts[2], height/2)
    if distance(p, q) > 1e-14
      let (c, h) = position_and_height(p, q),
          thickness = r_thickness + l_thickness,
          s = UnityNativeRef(@remote(b, RightCuboid(c, vz(1, c.cs), vx(1, c.cs), height, thickness, h, 0)))
        if length(verts) > 2
          (s, sweep_fractions(b, verts[2:end], height, l_thickness, r_thickness)...)
        else
          (s, )
        end
      end
    else
      sweep_fractions(b, verts[2:end], height, l_thickness, r_thickness)
    end
  end

backend_wall(b::Unity, path, height, l_thickness, r_thickness, family) =
  path_length(path) < path_tolerance() ?
    UnityEmptyRef() :
    begin
      @remote(b, SetCurrentMaterial(realize(b, family)))
      backend_wall_path(
          b,
          path,
          height*0.999, #We reduce height just a bit to avoid Z-fighting
          l_thickness, r_thickness)
    end

backend_wall_path(b::Unity, path::OpenPolygonalPath, height, l_thickness, r_thickness) =
    UnityUnionRef(sweep_fractions(b, path.vertices, height, l_thickness, r_thickness))

backend_wall_path(b::Unity, path::Path, height, l_thickness, r_thickness) =
    backend_wall_path(b, convert(OpenPolygonalPath, path), height, l_thickness, r_thickness)
=#

backend_wall(b::Unity, w_path, w_height, l_thickness, r_thickness, family) =
  path_length(w_path) < path_tolerance() ?
    UnityEmptyRef() :
    let w_paths = subpaths(w_path),
        r_w_paths = subpaths(offset(w_path, -r_thickness)),
        l_w_paths = subpaths(offset(w_path, l_thickness)),
        w_height = w_height*wall_z_fighting_factor,
        prevlength = 0,
        material = realize(b, family),
        refs = UnityNativeRef[]
      for (w_seg_path, r_w_path, l_w_path) in zip(w_paths, r_w_paths, l_w_paths)
        let currlength = prevlength + path_length(w_seg_path),
            c_r_w_path = closed_path_for_height(r_w_path, w_height),
            c_l_w_path = closed_path_for_height(l_w_path, w_height)
          push!(refs, realize_pyramid_fustrum(b, c_l_w_path, c_r_w_path, material))
          prevlength = currlength
        end
      end
      refs
    end

realize(b::Unity, w::Wall) =
  let w_base_height = w.bottom_level.height,
      w_height = (w.top_level.height - w_base_height)*wall_z_fighting_factor,
      r_thickness = r_thickness(w),
      l_thickness = l_thickness(w),
      w_path = translate(w.path, vz(w_base_height)),
      w_paths = subpaths(w_path),
      r_w_paths = subpaths(offset(w_path, -r_thickness)),
      l_w_paths = subpaths(offset(w_path, l_thickness)),
      openings = [w.doors..., w.windows...],
      prevlength = 0,
      material = realize(b, w.family),
      refs = UnityNativeRef[]
    for (w_seg_path, r_w_path, l_w_path) in zip(w_paths, r_w_paths, l_w_paths)
      let currlength = prevlength + path_length(w_seg_path),
          c_r_w_path = closed_path_for_height(r_w_path, w_height),
          c_l_w_path = closed_path_for_height(l_w_path, w_height)
        push!(refs, realize_pyramid_fustrum(b, c_l_w_path, c_r_w_path, material))
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
              push!(refs, realize_pyramid_fustrum(b, c_r_op_path, c_l_op_path, material))
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
        #realize_polygon(b, modleft, modright, modright, c_l_w_path, false)
        #realize_polygon(b, modleft, modright, modright, c_r_w_path, true)
      end
    end
    refs
  end

realize_pyramid_fustrum(b::Unity, bot_path::Path, top_path::Path, material) =
  realize_pyramid_fustrum(b, path_vertices(bot_path), path_vertices(top_path), material)

realize_pyramid_fustrum(b::Unity, bot_vs, top_vs, material) =
    UnityNativeRef(@remote(b, PyramidFrustumWithMaterial(bot_vs, top_vs, material)))

############################################
#=
backend_bounding_box(b::Unity, shapes::Shapes) =
  @remote(b, BoundingBox(collect_ref(shapes)))
=#

set_view(camera::Loc, target::Loc, lens::Real, b::Unity) =
  let c = connection(b)
    @remote(b, SetView(camera, target, lens))
    interrupt_processing(c)
  end

get_view(b::Unity) =
    (@remote(b, ViewCamera()), @remote(b, ViewTarget()), @remote(b, ViewLens()))

zoom_extents(b::Unity) = @remote(b, ZoomExtents())

view_top(b::Unity) = @remote(b, ViewTop())


delete_all_shapes(b::Unity) = @remote(b, DeleteAll())

backend_delete_shapes(b::Unity, shapes::Shapes) =
  @remote(b, DeleteMany(collect_ref(shapes)))

set_length_unit(unit::String, b::Unity) = nothing # Unused, for now

#=
# Dimensions

const UnityDimensionStyles = Dict(:architectural => "_ARCHTICK", :mechanical => "")

dimension(p0::Loc, p1::Loc, p::Loc, scale::Real, style::Symbol, b::Unity=current_backend()) =
    @remote(b, CreateAlignedDimension(p0, p1, p,
        scale,
        UnityDimensionStyles[style]))

dimension(p0::Loc, p1::Loc, sep::Real, scale::Real, style::Symbol, b::Unity=current_backend()) =
    let v = p1 - p0
        angle = pol_phi(v)
        dimension(p0, p1, add_pol(p0, sep, angle + pi/2), scale, style, b)
    end

=#

# Layers
# Experiment for multiple, simultaneous, alternative layers
# Layers







UnityLayer = Int

current_layer(b::Unity)::UnityLayer =
  @remote(b, CurrentParent())

current_layer(layer::UnityLayer, b::Unity) =
  @remote(b, SetCurrentParent(layer))

create_layer(name::String, b::Unity) =
  @remote(b, CreateParent(name))

create_layer(name::String, color::RGB, b::Unity) =
  let layer = @remote(b, CreateParent(name))
    @warn "Ignoring color in create_layer for Unity"
    #@remote(b, SetLayerColor(layer, color.r, color.g, color.b))
    layer
  end


set_layer_active(layer::UnityLayer, status, b::Unity) =
  let c = connection(b)
    @remote(b, SetActive(layer, status))
    interrupt_processing(c)
  end

delete_all_shapes_in_layer(layer::UnityLayer, b::Unity) =
  @remote(b, DeleteAllInParent(layer))

switch_to_layer(layer::UnityLayer, b::Unity) =
  @remote(b, SwitchToParent(layer))

# To preserve interactiveness during background


preserving_interactiveness(f, b::Unity=current_backend()) =
  let prev = @remote(b, SetMaxNonInteractiveRequests(0))
    f()
    @remote(b, SetMaxNonInteractiveRequests(prev))
  end

# Experiment to speed up things

canonicalize_layer(layer::UnityLayer, b::Unity) =
  @remote(b, Canonicalize(layer))

# Materials

UnityMaterial = Int

current_material(b::Unity)::UnityMaterial =
  @remote(b, CurrentMaterial())

current_material(material::UnityMaterial, b::Unity) =
  @remote(b, SetCurrentMaterial(material))

get_material(name::String, b::Unity) =
  @remote(b, LoadMaterial(name))


# Blocks




realize(b::Unity, s::Block) =
  s.shapes == [] ?
    @remote(b, LoadResource(s.name)) :
    @remote(b, Canonicalize(@remote(b, CreateBlockFromShapes(s.name, collect_ref(s.shapes)))))

realize(b::Unity, s::BlockInstance) =
    @remote(b, CreateBlockInstance(
        ref(s.block).value,
        s.loc, vy(1, s.loc.cs), vz(1, s.loc.cs), s.scale))
#=

# Manual process
@time for i in 1:1000 for r in 1:10 circle(x(i*10), r) end end

# Create block...
Khepri.create_block("Foo", [circle(radius=r) for r in 1:10])

# ...and instantiate it
@time for i in 1:1000 Khepri.instantiate_block("Foo", x(i*10), 0) end

=#

# Lights


backend_pointlight(b::Unity, loc::Loc, color::RGB, range::Real, intensity::Real) =
    @remote(b, PointLight(loc, color, range, intensity))
#=
backend_spotlight(b::Unity, loc::Loc, dir::Vec, hotspot::Real, falloff::Real) =
    @remote(b, SpotLight(loc, hotspot, falloff, loc + dir))

backend_ieslight(b::Unity, file::String, loc::Loc, dir::Vec, alpha::Real, beta::Real, gamma::Real) =
    @remote(b, IESLight(file, loc, loc + dir, vxyz(alpha, beta, gamma)))

# User Selection
=#

shape_from_ref(r, b::Unity) =
  let idx = findfirst(s -> r in collect_ref(s), collected_shapes())
    if isnothing(idx)
          unknown(r, backend=b, ref=LazyRef(b, UnityNativeRef(r), 0, 0))
          #code = @remote(b, ShapeCode(r)),
          #ref = LazyRef(b, UnityNativeRef(r))
          #error("Unknown shape with code $(code)")
    else
      collected_shapes()[idx]
    end
  end
#
#=


select_position(prompt::String, b::Unity) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = @remote(b, GetPosition(prompt))
      length(ans) > 0 && ans[1]
    end
  end

select_with_prompt(prompt::String, b::Backend, f::Function) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = f(connection(b), prompt)
      length(ans) > 0 && shape_from_ref(ans[1], b)
    end
  end



# HACK: The next operations should receive a set of shapes to avoid re-creating already existing shapes

select_point(prompt::String, b::Unity) =
  select_with_prompt(prompt, b, UnityGetPoint)



select_curve(prompt::String, b::Unity) =
  select_with_prompt(prompt, b, UnityGetCurve)



select_surface(prompt::String, b::Unity) =
  select_with_prompt(prompt, b, UnityGetSurface)



select_solid(prompt::String, b::Unity) =
  select_with_prompt(prompt, b, UnityGetSolid)



select_shape(prompt::String, b::Unity) =
  select_with_prompt(prompt, b, UnityGetShape)




captured_shape(b::Unity, handle) =
  shape_from_ref(@remote(b, GetShapeFromHandle(handle)),
                 b)

generate_captured_shape(s::Shape, b::Unity) =
    println("captured_shape(autocad, $(@remote(b, GetHandleFromShape(ref(s).value)))"))

# Register for notification








register_for_changes(s::Shape, b::Unity) =
    negin
        @remote(b, RegisterForChanges(ref(s).value))
        @remote(b, DetectCancel())
        s
    end

unregister_for_changes(s::Shape, b::Unity) =
    begin
        @remote(b, UnregisterForChanges(ref(s).value))
        @remote(b, UndetectCancel())
        s
    end

waiting_for_changes(s::Shape, b::Unity) =
    ! @remote(b, WasCanceled())

changed_shape(ss::Shapes, b::Unity) =
    let changed = []
        while length(changed) == 0 && ! @remote(b, WasCanceled())
            changed =  @remote(b, ChangedShape())
            sleep(0.1)
        end
        if length(changed) > 0
            shape_from_ref(changed[1], b)
        else
            nothing
        end
    end




# HACK: This should be filtered on the plugin, not here.
all_shapes(b::Unity) =
  Shape[shape_from_ref(r, b)
        for r in filter(r -> @remote(b, ShapeCode(r) != 0, @remote(b, GetAllShapes())))]

all_shapes_in_layer(layer, b::Unity) =
  Shape[shape_from_ref(r, b) for r in @remote(b, GetAllShapesInLayer(layer))]

disable_update(b::Unity) =
  @remote(b, DisableUpdate())

enable_update(b::Unity) =
  @remote(b, EnableUpdate())

# Render

=#



#render exposure: [-3, +3] -> [-6, 21]
convert_render_exposure(b::Unity, v::Real) = -4.05*v + 8.8
#render quality: [-1, +1] -> [+1, +50]
convert_render_quality(b::Unity, v::Real) = round(Int, 25.5 + 24.5*v)

render_view(path::String, b::Unity) =
    let c = connection(b)
      @remote(b, SetResolution(render_width(), render_height()))
      interrupt_processing(c)
      @remote(b, ScreenShot(path))
    end



highlight_shape(s::Shape, b::Unity) =
    @remote(b, SelectGameObjects(collect_ref(s)))

highlight_shapes(ss::Shapes, b::Unity) =
    @remote(b, SelectGameObjects(collect_ref(ss)))


select_position(prompt::String, b::Unity) =
  begin
    @info "$(prompt) on the $(b) backend."
    @remote(b, StartSelectingPosition())
    let s = u0() # Means not found
      while s == u0()
        sleep(0.1)
        s = @remote(b, SelectedPosition())
      end
      s
    end
  end

select_shape(prompt::String, b::Unity) =
  select_one_with_prompt(prompt, b, (c, prompt) ->
    let s = -2 # Means not found
      @remote(b, StartSelectingGameObject())
      while s == -2
        sleep(0.1)
        s = @remote(b, SelectedGameObjectId(true))
      end
      [s]
    end)
