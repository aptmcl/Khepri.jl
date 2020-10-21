export rhino

macro rhino_str(str)
    rpc("RH", str)
end

# We need some additional Encoders
encode_Brep = encode_Guid
decode_Brep = decode_Guid
decode_Brep_array = decode_Guid_array
encode_RhinoObject = encode_Guid
encode_RhinoObject_array = encode_Guid_array
decode_RhinoObject = decode_Guid

encode_Plane(c::IO, v::XYZ) =
  let o = in_world(v)
      vx = in_world(vx(1, c.cs))
      vy = in_world(vy(1, c.cs))
    encode_Point3d(c, o)
    encode_Vector3d(c, vx)
    encode_Vector3d(c, vy)
  end
decode_Plane(c::IO) =
    loc_from_o_vx_vy(decode_Point3d(c), decode_Vector3d(c), decode_Vector3d(c))

encode_TableFamily = encode_int
decode_TableFamily = decode_int
encode_ChairFamily = encode_int
decode_ChairFamily = decode_int
encode_TableChairFamily = encode_int
decode_TableChairFamily = decode_int

rhino_api = @remote_functions :CS """
public void RunScript(string script)
public void SetView(Point3d position, Point3d target, double lens, bool perspective, string style)
public void View(Point3d position, Point3d target, double lens)
public void ViewTop()
public Point3d ViewCamera()
public Point3d ViewTarget()
public double ViewLens()
public Guid Point(Point3d p)
public Point3d PointPosition(Guid ent)
public Guid PolyLine(Point3d[] pts)
public Point3d[] LineVertices(RhinoObject id)
public Guid Spline(Point3d[] pts)
public Guid ClosedPolyLine(Point3d[] pts)
public Guid ClosedSpline(Point3d[] pts)
public Guid Circle(Point3d c, Vector3d n, double r)
public Point3d CircleCenter(RhinoObject obj)
public Vector3d CircleNormal(RhinoObject obj)
public double CircleRadius(RhinoObject obj)
public Guid Ellipse(Point3d c, Vector3d n, double radiusX, double radiusY)
public Guid Arc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)
public Guid JoinCurves(Guid[] objs)
public Guid SurfaceCircle(Point3d c, Vector3d n, double r)
public Guid SurfaceEllipse(Point3d c, Vector3d n, double radiusX, double radiusY)
public Guid SurfaceArc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)
public Guid SurfaceClosedPolyLine(Point3d[] pts)
public Guid SurfaceFrom(Guid[] objs)
public Guid Sphere(Point3d c, double r)
public Guid Torus(Point3d c, Vector3d vz, double majorRadius, double minorRadius)
public Brep Cylinder(Point3d bottom, double radius, Point3d top)
public Brep Cone(Point3d bottom, double radius, Point3d top)
public Brep ConeFrustum(Point3d bottom, double bottom_radius, Point3d top, double top_radius)
public Brep Box(Point3d corner, Vector3d vx, Vector3d vy, double dx, double dy, double dz)
public Brep XYCenteredBox(Point3d corner, Vector3d vx, Vector3d vy, double dx, double dy, double dz)
public Brep IrregularPyramid(Point3d[] pts, Point3d apex)
public Brep IrregularPyramidFrustum(Point3d[] bpts, Point3d[] tpts)
public Guid SurfaceFromGrid(int nU, int nV, Point3d[] pts, bool closedU, bool closedV, int degreeU, int degreeV)
public Brep[] Thicken(RhinoObject obj, double thickness)
public double[] CurveDomain(RhinoObject obj)
public double CurveLength(RhinoObject obj)
public Plane CurveFrameAt(RhinoObject obj, double t)
public Plane CurveFrameAtLength(RhinoObject obj, double l)
public Brep Extrusion(RhinoObject obj, Vector3d dir)
public Brep[] SweepPathProfile(RhinoObject path, RhinoObject profile, double rotation, double scale)
public Guid[] Intersect(RhinoObject obj0, RhinoObject obj1)
public Guid[] Subtract(RhinoObject obj0, RhinoObject obj1)
public Guid Move(Guid id, Vector3d v)
public Guid Scale(Guid id, Point3d p, double s)
public Guid Rotate(Guid id, Point3d p, Vector3d n, double a)
public Guid Mirror(Guid id, Point3d p, Vector3d n, bool copy)
public void Render(int width, int height, string path)
public void SaveView(int width, int height, string path)
public bool IsPoint(RhinoObject e)
public bool IsCircle(RhinoObject e)
public bool IsPolyLine(RhinoObject e)
public bool IsSpline(RhinoObject e)
public bool IsInterpSpline(RhinoObject e)
public bool IsClosedPolyLine(RhinoObject e)
public bool IsClosedSpline(RhinoObject e)
public bool IsInterpClosedSpline(RhinoObject e)
public bool IsEllipse(RhinoObject e)
public bool IsArc(RhinoObject e)
public bool IsText(RhinoObject e)
public byte ShapeCode(RhinoObject id)
public Guid ClosedPathCurveArray(Point3d[] pts, double[] angles)
public Brep[] PathWall(RhinoObject obj, double lThickness, double rThickness, double height)
public Brep RectangularTable(Point3d c, double angle, double length, double width, double height, double top_thickness, double leg_thickness)
public Brep[] BaseRectangularTable(double length, double width, double height, double top_thickness, double leg_thickness)
public int CreateRectangularTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)
public Guid Table(Point3d c, double angle, int family)
public Brep[] BaseChair(double length, double width, double height, double seat_height, double thickness)
public int CreateChairFamily(double length, double width, double height, double seat_height, double thickness)
public Guid Chair(Point3d c, double angle, int family)
public GeometryBase[] BaseRectangularTableAndChairs(int tableFamily, int chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)
public int CreateRectangularTableAndChairsFamily(int tableFamily, int chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)
public Guid TableAndChairs(Point3d c, double angle, int family)
public String CreateLayer(String name, bool active, byte r, byte g, byte b)
public void SetLayerColor(ObjectId id, byte r, byte g, byte b)
public void SetShapeColor(ObjectId id, byte r, byte g, byte b)
public String CurrentLayer()
public void SetCurrentLayer(String id)
public void SetLayerActive(String name, bool active)
public String ShapeLayer(RhinoObject objId)
public void SetShapeLayer(RhinoObject objId, String layerId)
public Point3d[] CurvePointsAt(RhinoObject obj, double[] ts)
public Vector3d[] CurveTangentsAt(RhinoObject obj, double[] ts)
public Vector3d[] CurveNormalsAt(RhinoObject obj, double[] ts)
public int DeleteAll()
public int DeleteAllInLayer(String name)
public void Delete(Guid id)
public void DeleteMany(Guid[] ids)
public Guid Clone(RhinoObject e)
public Point3d[] GetPosition(string prompt)
public Guid[] GetPoint(string prompt)
public Guid[] GetPoints(string prompt)
public Guid[] GetCurve(string prompt)
public Guid[] GetCurves(string prompt)
public Guid[] GetSurface(string prompt)
public Guid[] GetSurfaces(string prompt)
public Guid[] GetSolid(string prompt)
public Guid[] GetSolids(string prompt)
public Guid[] GetShape(string prompt)
public Guid[] GetShapes(string prompt)
public Guid[] GetAllShapes()
public Guid[] GetAllShapesInLayer(String name)
"""
#=
public byte Sync()
public byte Disconnect()
public Entity InterpSpline(Point3d[] pts, Vector3d tan0, Vector3d tan1)
public Entity InterpClosedSpline(Point3d[] pts)
public Entity Text(string str, Point3d corner, Vector3d vx, Vector3d vy, double height)
public Guid SurfaceFromCurve(Entity curve)
public Guid IrregularPyramidMesh(Point3d[] pts, Point3d apex)
public Entity SolidFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level, double thickness)
public Guid NurbSurfaceFrom(ObjectId id)
public double[] SurfaceDomain(Entity ent)
public Frame3d SurfaceFrameAt(Entity ent, double u, double v)
public Guid Sweep(ObjectId pathId, ObjectId profileId, double rotation, double scale)
public Guid Loft(ObjectId[] profilesIds, ObjectId[] guidesIds, bool ruled, bool closed)
public void Unite(ObjectId objId0, ObjectId objId1)
public Guid Revolve(ObjectId profileId, Point3d p, Vector3d n, double startAngle, double amplitude)
public Point3d[] GetPoint(string prompt)
public Point3d[] BoundingBox(ObjectId[] ids)
public void ZoomExtents()
public void SetSystemVariableInt(string name, int value)
public int Command(string cmd)
public void DisableUpdate()
public void EnableUpdate()
public BIMLevel FindOrCreateLevelAtElevation(double elevation)
public BIMLevel UpperLevel(BIMLevel currentLevel, double addedElevation)
public double GetLevelElevation(BIMLevel level)
public Entity MeshFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN)

public FloorFamily FloorFamilyInstance(double totalThickness, double coatingThickness)
public Entity LightweightPolyLine(Point2d[] pts, double[] angles, double elevation)
public Entity SurfaceLightweightPolyLine(Point2d[] pts, double[] angles, double elevation)
public ObjectId CreatePathFloor(Point2d[] pts, double[] angles, BIMLevel level, FloorFamily family)
#public TableFamily FindOrCreateTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)
#public ChairFamily FindOrCreateChairFamily(double length, double width, double height, double seat_height, double thickness)
#public TableChairFamily FindOrCreateTableChairFamily(TableFamily tableFamily, ChairFamily chairFamily, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)
#public Brep RectangularTable(Point3d c, double angle, TableFamily family)
#public Brep Chair(Point3d c, double angle, ChairFamily family)
#public Guid InstanceChair(Point3d c, double angle, ChairFamily family)
#public Brep[] RowOfChairs(Point3d c, double angle, int n, double spacing, ChairFamily family)
#public Brep[] CenteredRowOfChairs(Point3d c, double angle, int n, double spacing, ChairFamily family)
#public Brep[] RectangularTableAndChairs(Point3d c, double angle, TableChairFamily f)
=#


abstract type RHKey end
const RHId = Guid
const RHRef = GenericRef{RHKey, RHId}
const RHRefs = Vector{RHRef}
const RHEmptyRef = EmptyRef{RHKey, RHId}
const RHNativeRef = NativeRef{RHKey, RHId}
const RHUnionRef = UnionRef{RHKey, RHId}
const RHSubtractionRef = SubtractionRef{RHKey, RHId}
const RH = SocketBackend{RHKey, RHId}

void_ref(b::RH) = RHNativeRef(zeros(UInt8, 16))

create_RH_connection() = create_backend_connection("Rhinoceros", 12000)

const rhino = RH(LazyParameter(TCPSocket, create_RH_connection), rhino_api)

# This should not be done automatically
backend_name(b::RH) = "Rhino"

realize(b::RH, s::EmptyShape) =
  RHEmptyRef()
realize(b::RH, s::UniversalShape) =
  RHUniversalRef()

backend_stroke(b::RH, path::CircularPath) =
    @remote(b, Circle(path.center, vz(1, path.center.cs), path.radius))
backend_stroke(b::RH, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        @remote(b, ClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)]))
    end
backend_stroke(b::RH, path::ArcPath) =
  backend_stroke_arc(b, path.center, path.radius, path.start_angle, path.amplitude)

backend_stroke(b::RH, path::OpenPolygonalPath) =
  	@remote(b, PolyLine(path.vertices))
backend_stroke(b::RH, path::ClosedPolygonalPath) =
    @remote(b, ClosedPolyLine(path.vertices))

backend_stroke(b::RH, path::OpenSplinePath) =
  if (path.v0 == false) && (path.v1 == false)
    #ACADSpline(connection(b), path.vertices)
    @remote(b, Spline(path.vertices))
  elseif (path.v0 != false) && (path.v1 != false)
    @remote(b, InterpSpline(path.vertices, path.v0, path.v1))
  else
    @remote(b, InterpSpline(connection(b),
                     path.vertices,
                     path.v0 == false ? path.vertices[2]-path.vertices[1] : path.v0,
                     path.v1 == false ? path.vertices[end-1]-path.vertices[end] : path.v1))
  end
backend_stroke(b::RH, path::ClosedSplinePath) =
    @remote(b, InterpClosedSpline(path.vertices))
backend_fill(b::RH, path::ClosedSplinePath) =
    backend_fill_curves(b, @remote(b, InterpClosedSpline(path.vertices)))

backend_stroke_unite(b::RH, refs) = @remote(b, JoinCurves(refs))

#=backend_fill(b::RH, path::ClosedPolygonalPath) =
    @remote(b, SurfaceClosedPolyLine(path.vertices))
    backend_fill(b::RH, path::RectangularPath) =
        let c = path.corner,
            dx = path.dx,
            dy = path.dy
            @remote(b, SurfaceClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)]))
        end
backend_fill(b::RH, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        SurfaceClosedPolyLine(connection(b), [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)])
    end
=#

backend_fill_curves(b::RH, gs::Guids) = @remote(b, SurfaceFrom(gs))
backend_fill_curves(b::RH, g::Guid) = @remote(b, SurfaceFrom([g]))

backend_stroke_line(b::RH, vs) = @remote(b, PolyLine(vs))

backend_stroke_arc(b::RH, center::Loc, radius::Real, start_angle::Real, amplitude::Real) =
    let end_angle = start_angle + amplitude
        if end_angle > start_angle
            @remote(b, Arc(center, vz(1, center.cs), radius, start_angle, end_angle))
        else
            @remote(b, Arc(center, vz(1, center.cs), radius, end_angle, start_angle))
        end
    end

realize(b::RH, s::Point) =
  @remote(b, Point(s.position))

realize(b::RH, s::Line) =
  @remote(b, PolyLine(s.vertices))

realize(b::RH, s::Spline) =
  if (s.v0 == false) && (s.v1 == false)
    @remote(b, Spline(s.points))
  elseif (s.v0 != false) && (s.v1 != false)
    @remote(b, InterpSpline(s.points, s.v0, s.v1))
  else
    @remote(b, InterpSpline(connection(b),
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1))
  end

realize(b::RH, s::ClosedSpline) =
  @remote(b, InterpClosedSpline(s.points))

realize(b::RH, s::Circle) =
  @remote(b, Circle(s.center, vz(1, s.center.cs), s.radius))

realize(b::RH, s::Arc) =
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

realize(b::RH, s::Ellipse) =
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), s.radius_x, s.radius_y))

realize(b::RH, s::EllipticArc) =
  error("Finish this")

realize(b::RH, s::Polygon) =
  @remote(b, ClosedPolyLine(s.vertices))

realize(b::RH, s::RegularPolygon) =
  @remote(b, ClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))

realize(b::RH, s::Rectangle) =
  @remote(b, ClosedPolyLine([s.corner, add_x(s.corner, s.dx), add_xy(s.corner, s.dx, s.dy), add_y(s.corner, s.dy)]))

backend_map_division(b::RH, f::Function, s::Shape1D, n::Int) =
  let conn = connection(b),
      r = ref(s).value,
      (t1, t2) = @remote(b, CurveDomain(r)),
      ti = division(t1, t2, n),
      ps = @remote(b, CurvePointsAt(r, ti)),
      ts = @remote(b, CurveTangentsAt(r, ti)),
      #ns = ACADCurveNormalsAt(conn, r, ti),
      frames = rotation_minimizing_frames(@remote(b, CurveFrameAt(r, t1)), ps, ts)
    map(f, frames)
  end

realize(b::RH, s::SurfaceCircle) =
  @remote(b, SurfaceCircle(s.center, vz(1, s.center.cs), s.radius))

realize(b::RH, s::SurfaceArc) =
  @remote(b, SurfaceArc(s.center, vz(1, s.center.cs), s.radius, s.start_angle, s.start_angle + s.amplitude))

realize(b::RH, s::SurfaceEllipticArc) =
  error("Finish this")

realize(b::RH, s::SurfaceEllipse) =
  @remote(b, SurfaceEllipse(s.center, vz(1, s.center.cs), s.radius_x, s.radius_y))

backend_surface_polygon(b::RH, vs::Locs) =
  @remote(b, SurfaceClosedPolyLine(vs))

realize(b::RH, s::Surface) =
  let ids = @remote(b, SurfaceFrom(collect_ref(s.frontier)))
    foreach(mark_deleted, s.frontier)
    ids
  end
#=
realize(b::RH, s::Text) =
  @remote(b, Text(s.str, s.c, vx(1, s.c.cs), vy(1, s.c.cs), s.h))
=#

backend_surface_domain(b::RH, s::Shape2D) =
    tuple(@remote(b, SurfaceDomain(ref(s).value)...))

backend_map_division(b::RH, f::Function, s::Shape2D, nu::Int, nv::Int) =
    let conn = connection(b)
        r = ref(s).value
        (u1, u2, v1, v2) = @remote(b, SurfaceDomain(r))
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(@remote(b, SurfaceFrameAt(r, u, v)))
            end
        end
    end

realize(b::RH, s::Sphere) =
    @remote(b, Sphere(s.center, s.radius))
realize(b::RH, s::Torus) =
    @remote(b, Torus(s.center, vz(1, s.center.cs), s.re, s.ri))
realize(b::RH, s::Cuboid) =
    @remote(b, IrregularPyramidFrustum([s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3]))
realize(b::RH, s::RegularPyramidFrustum) =
    @remote(b, IrregularPyramidFrustum(regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                                       regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed)))
backend_pyramid(b::RH, bs::Locs, t::Loc) =
  @remote(b, IrregularPyramid(bs, t))
backend_pyramid_frustum(b::RH, bs::Locs, ts::Locs) =
  @remote(b, IrregularPyramidFrustum(bs, ts))
backend_right_cuboid(b::RH, cb, width, height, h, material) =
  @remote(b, XYCenteredBox(cb, vx(1, cb.cs), vy(1, cb.cs), width, height, h))

realize(b::RH, s::Box) =
    @remote(b, Box(s.c, vx(1,s.c.cs), vy(1,s.c.cs), s.dx, s.dy, s.dz))

realize(b::RH, s::Cone) =
  @remote(b, Cone(add_z(s.cb, s.h), s.r, s.cb))

realize(b::RH, s::ConeFrustum) =
  @remote(b, ConeFrustum(s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt))

realize(b::RH, s::Cylinder) =
  @remote(b, Cylinder(s.cb, s.r, s.cb + vz(s.h, s.cb.cs)))

backend_extrusion(b::RH, s::Shape, v::Vec) =
    and_mark_deleted(
        map_ref(s) do r
            @remote(b, Extrusion(r, v))
        end,
        s)

backend_sweep(b::RH, path::Shape, profile::Shape, rotation::Real, scale::Real) =
  map_ref(profile) do profile_r
    map_ref(path) do path_r
      @remote(b, SweepPathProfile(path_r, profile_r, rotation, scale))
    end
  end

#=
realize(b::RH, s::Revolve) =
  and_delete_shape(
    map_ref(s.profile) do r
      @remote(b, Revolve(r, s.p, s.n, s.start_angle, s.amplitude))
    end,
    s.profile)

backend_loft_points(b::Backend, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  let f = (ruled ? (closed ? polygon : line) : (closed ? closed_spline : spline))
    and_delete_shapes(ref(f(map(point_position, profiles), backend=b)),
                      vcat(profiles, rails))
  end

backend_loft_curves(b::RH, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  and_delete_shapes(RHLoft(connection(b),
                             collect_ref(profiles),
                             collect_ref(rails),
                             ruled, closed),
                    vcat(profiles, rails))

backend_loft_surfaces(b::RH, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  backend_loft_curves(b, profiles, rails, ruled, closed)

backend_loft_curve_point(b::RH, profile::Shape, point::Shape) =
  and_delete_shapes(RHLoft(connection(b),
                             vcat(collect_ref(profile), collect_ref(point)),
                             [],
                             true, false),
                    [profile, point])

backend_loft_surface_point(b::RH, profile::Shape, point::Shape) =
  backend_loft_curve_point(b, profile, point)
=#
unite_refs(b::RH, refs::Vector{<:RHRef}) =
    RHUnionRef(tuple(refs...))

unite_refs(b::RH, r::RHUnionRef) =
  r

unite_ref(b::RH, r0::RHNativeRef, r1::RHNativeRef) =
  (@remote(b, Unite(r0.value, r1.value)); r0)

intersect_ref(b::RH, r0::RHNativeRef, r1::RHNativeRef) =
    let refs = @remote(b, Intersect(r0.value, r1.value))
        n = length(refs)
        if n == 0
            RHEmptyRef()
        elseif n == 1
            RHNativeRef(refs[1])
        else
            RHUnionRef(map(RHNativeRef, tuple(refs...)))
        end
    end

subtract_ref(b::RH, r0::RHNativeRef, r1::RHNativeRef) =
  let refs = try @remote(b, Subtract(r0.value, r1.value)); catch e; [nothing] end,
      n = length(refs)
    if n == 0
        RHEmptyRef()
    elseif n == 1
        if isnothing(refs[1]) # failed
            RHSubtractionRef(r0, tuple(r1))
        else
            RHNativeRef(refs[1])
        end
    else
        RHUnionRef(map(RHNativeRef, tuple(refs...)))
    end
  end

subtract_ref(b::RH, r0::SubtractionRef, r1::RHNativeRef) =
  subtract_ref(b,
    subtract_ref(b, r0.value, r1),
    length(r0.values) == 1 ? r0.values[1] : RHUnionRef(r0.values))

subtract_ref(b::RH, r0::RHRef, r1::RHUnionRef) =
  foldl((r0,r1)->subtract_ref(b,r0,r1), r1.values, init=r0)

#=
slice_ref(b::RH, r::RHNativeRef, p::Loc, v::Vec) =
  (@remote(b, Slice(r.value, p, v); r))

slice_ref(b::RH, r::RHUnionRef, p::Loc, v::Vec) =
  map(r->slice_ref(b, r, p, v), r.values)

=#

#=
realize(b::RH, s::IntersectionShape) =
  foldl((r0,r1)->intersect_ref(b,r0,r1), UniversalRef{RHId}(), map(ref, s.shapes))
  =#

realize(b::RH, s::Slice) =
  slice_ref(b, ref(s.shape), s.p, s.n)

realize(b::RH, s::Move) =
  let r = map_ref(s.shape) do r
            @remote(b, Move(r, s.v))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::RH, s::Scale) =
  let r = map_ref(s.shape) do r
            @remote(b, Scale(r, s.p, s.s))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::RH, s::Rotate) =
  let r = map_ref(s.shape) do r
            @remote(b, Rotate(r, s.p, s.v, s.angle))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::RH, s::Mirror) =
  and_delete_shape(map_ref(s.shape) do r
                    @remote(b, Mirror(r, s.p, s.n, false))
                   end,
                   s.shape)

realize(b::RH, s::UnionMirror) =
  let r0 = ref(s.shape),
      r1 = map_ref(r0) do r
            @remote(b, Mirror(r, s.p, s.n, true))
          end
    UnionRef((r0,r1))
  end

backend_surface_grid(b::RH, points, closed_u, closed_v, smooth_u, smooth_v) =
  let (nu, nv) = size(points),
      order(n) = min(max(2*floor(Int,n/30) + 1, 2), 5)
    @remote(b, SurfaceFromGrid(nu, nv,
                               reshape(permutedims(points), :),
                               closed_u, closed_v,
                               smooth_u ? order(nu) : 1,
                               smooth_v ? order(nv) : 1))
  end

realize(b::RH, s::Thicken) =
  and_delete_shape(
    map_ref(s.shape) do r
      @remote(b, Thicken(r, s.thickness))
    end,
    s.shape)

# BIM

# Families

realize(b::RH, f::TableFamily) =
    @remote(b, CreateRectangularTableFamily(f.length, f.width, f.height, f.top_thickness, f.leg_thickness))
realize(b::RH, f::ChairFamily) =
    @remote(b, CreateChairFamily(f.length, f.width, f.height, f.seat_height, f.thickness))
realize(b::RH, f::TableChairFamily) =
    @remote(b, CreateRectangularTableAndChairsFamily(
        realize(b, f.table_family), realize(b, f.chair_family),
        f.table_family.length, f.table_family.width,
        f.chairs_top, f.chairs_bottom, f.chairs_right, f.chairs_left,
        f.spacing))

backend_rectangular_table(b::RH, c, angle, family) =
    @remote(b, Table(c, angle, realize(b, family)))

backend_chair(b::RH, c, angle, family) =
    @remote(b, Chair(c, angle, realize(b, family)))

backend_rectangular_table_and_chairs(b::RH, c, angle, family) =
    @remote(b, TableAndChairs(c, angle, realize(b, family)))

backend_slab(b::RH, profile, holes, thickness, family) =
    map_ref(b,
            r->@remote(b, Extrusion(r, vz(thickness))),
            ensure_ref(b, backend_fill(b, profile)))

#=
#Beams are aligned along the top axis.
realize(b::RH, s::Beam) =
  let profile = s.family.profile
      profile_u0 = profile.corner
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
      # need to test whether it is rotation on center or on axis
      o = add_y(loc_from_o_phi(s.cb, s.angle), -profile.dy/2)
    @remote(b, XYCenteredBox(o, vx(1, o.cs), vy(1, o.cs), profile.dx, profile.dy, s.h))
  end

#Columns are aligned along the center axis.
realize(b::RH, s::FreeColumn) =
  let profile = s.family.profile
      profile_u0 = profile.corner
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
      # need to test whether it is rotation on center or on axis
      o = loc_from_o_phi(c, s.angle)
    @remote(b, XYCenteredBox(o, vx(1, o.cs), vy(1, o.cs), profile.dx, profile.dy, s.h))
  end

realize(b::RH, s::Column) =
  let profile = s.family.profile,
      profile_u0 = profile.corner,
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2),
      base_height = s.bottom_level.height,
      height = s.top_level.height - base_height,
      # need to test whether it is rotation on center or on axis
      o = loc_from_o_phi(s.cb + vz(base_height), s.angle)
    @remote(b, XYCenteredBox(o, vx(1, o.cs), vy(1, o.cs), profile.dx, profile.dy, height))
  end
=#

backend_wall(b::RH, path, height, l_thickness, r_thickness, family) =
    @remote(b, PathWall(backend_stroke(b, offset(path, l_thickness - r_thickness)), #offset(path, (l_thickness - r_thickness)/2)),
                        l_thickness, r_thickness,
                        height))

############################################

backend_bounding_box(b::RH, shapes::Shapes) =
  @remote(b, BoundingBox(collect_ref(shapes)))

import Base.view
set_view(camera::XYZ, target::XYZ, lens::Real, aperture::Real, b::RH) =
  @remote(b, View(camera, target, lens))

get_view(b::RH) =
  @remote(b, ViewCamera()), @remote(b, ViewTarget()), @remote(b, ViewLens())

zoom_extents(b::RH) =
  @remote(b, ZoomExtents())

view_top(b::RH) =
  @remote(b, ViewTop())

backend_delete_shapes(b::RH, shapes::Shapes) =
  @remote(b, DeleteMany(collect_ref(shapes)))

delete_all_shapes(b::RH) =
  @remote(b, DeleteAll())

# Layers
RHLayer = String

current_layer(b::RH)::RHLayer =
  @remote(b, CurrentLayer())

current_layer(layer::RHLayer, b::RH) =
  @remote(b, SetCurrentLayer(layer))

backend_create_layer(b::RH, name::String, active::Bool, color::RGB) =
  let to255(x) = round(UInt8, x*255)
    @remote(b, CreateLayer(name, true, to255(red(color)), to255(green(color)), to255(blue(color))))
  end

delete_all_shapes_in_layer(layer::RHLayer, b::RH) =
  @remote(b, DeleteAllInLayer(layer))

shape_from_ref(r, b::RH) =
  let c = connection(b)
    let code = @remote(b, ShapeCode(r)),
        ref = LazyRef(b, RHNativeRef(r))
      if code == 1
        point(@remote(b, PointPosition(r)),
              backend=b, ref=ref)
      elseif code == 2
        circle(maybe_loc_from_o_vz(@remote(b, CircleCenter(r)), @remote(b, CircleNormal(r))),
               @remote(b, CircleRadius(r)),
               backend=b, ref=ref)
      elseif 3 <= code <= 6
        line(@remote(b, LineVertices(r)),
             backend=b, ref=ref)
      elseif code == 7
        spline([xy(0,0)], false, false, #HACK obtain interpolation points
               backend=b, ref=ref)
      elseif 103 <= code <= 106
        polygon(@remote(b, LineVertices(r)),
                backend=b, ref=ref)
      elseif code == 40
        # FIXME: frontier is missing
        surface(frontier=[], backend=b, ref=ref)
      elseif code == 41
        surface(frontier=[], backend=b, ref=ref)
      elseif code == 81
        sphere(backend=b, ref=ref)
      elseif code == 82
        cylinder(backend=b, ref=ref)
      elseif code == 83
        cone(backend=b, ref=ref)
      elseif code == 84
        torus(backend=b, ref=ref)
      else
        error("Unknown shape")
      end
    end
  end

#

select_position(prompt::String, b::RH) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = @remote(b, GetPosition(prompt))
      length(ans) > 0 ? ans[1] : nothing
    end
  end

select_positions(prompt::String, b::RH) =
  let ps = Loc[]
      p = nothing
    @info "$(prompt) on the $(b) backend."
    while ((p = @remote(b, GetPosition(prompt)) |> length)) > 0
        push!(ps, p...)
    end
    ps
  end

# HACK: The next operations should receive a set of shapes to avoid re-creating already existing shapes

select_point(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, @get_remote b GetPoint)

select_points(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, @get_remote b GetPoints)

select_curve(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, @get_remote b GetCurve)

select_curves(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, @get_remote b GetCurves)


select_surface(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, @get_remote b GetSurface)

select_surfaces(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, @get_remote b GetSurfaces)


select_solid(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, @get_remote b GetSolid)

select_solids(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, @get_remote b GetSolids)

select_shape(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, @get_remote b GetShape)

select_shapes(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, @get_remote b GetShapes)

captured_shape(b::RH, handle) =
  shape_from_ref(handle, b)
#
captured_shapes(b::RH, handles) =
  map(handles) do handle
      shape_from_ref(handle, b)
  end

generate_captured_shape(s::Shape, b::RH) =
    println("captured_shape(rhino, $(ref(s).value))")

generate_captured_shapes(ss::Shapes, b::RH) =
  begin
    print("captured_shapes(rhino, [")
    for s in ss
      print(ref(s).value)
      print(", ")
    end
    println("])")
  end

all_shapes(b::RH) =
  Shape[shape_from_ref(r, b)
        for r in filter(r -> @remote(b, ShapeCode(r)) != 0, @remote(b, GetAllShapes()))]

all_shapes_in_layer(layer, b::RH) =
  Shape[shape_from_ref(r, b) for r in @remote(b, GetAllShapesInLayer(layer))]

render_view(path::String, b::RH) =
  @remote(b, Render(render_width(), render_height(), path))

save_view(path::String, b::RH) =
  @remote(b, SaveView(render_width(), render_height(), path))


#=
BIM families for Rhino
=#

abstract type RhinoFamily <: Family end

struct RhinoLayerFamily <: RhinoFamily
  name::String
  color::RGB
  ref::Parameter{Any}
end

rhino_layer_family(name, color::RGB=rgb(1,1,1)) =
  RhinoLayerFamily(name, color, Parameter{Any}(nothing))

backend_get_family_ref(b::RH, f::Family, af::RhinoLayerFamily) =
  backend_create_layer(b, af.name, true, af.color)

set_backend_family(default_wall_family(), rhino, rhino_layer_family("Wall"))
set_backend_family(default_slab_family(), rhino, rhino_layer_family("Slab"))
set_backend_family(default_roof_family(), rhino, rhino_layer_family("Roof"))
set_backend_family(default_beam_family(), rhino, rhino_layer_family("Beam"))
set_backend_family(default_column_family(), rhino, rhino_layer_family("Column"))
set_backend_family(default_door_family(), rhino, rhino_layer_family("Door"))
set_backend_family(default_panel_family(), rhino, rhino_layer_family("Panel"))

set_backend_family(default_table_family(), rhino, rhino_layer_family("Table"))
set_backend_family(default_chair_family(), rhino, rhino_layer_family("Chair"))
set_backend_family(default_table_chair_family(), rhino, rhino_layer_family("TableChairs"))

set_backend_family(default_curtain_wall_family(), rhino, rhino_layer_family("CurtainWall"))
set_backend_family(default_curtain_wall_family().panel, rhino, rhino_layer_family("CurtainWall-Panel"))
set_backend_family(default_curtain_wall_family().boundary_frame, rhino, rhino_layer_family("CurtainWall-Boundary"))
set_backend_family(default_curtain_wall_family().transom_frame, rhino, rhino_layer_family("CurtainWall-Transom"))
set_backend_family(default_curtain_wall_family().mullion_frame, rhino, rhino_layer_family("CurtainWall-Mullion"))
#current_backend(rhino)

set_backend_family(default_truss_node_family(), rhino, rhino_layer_family("TrussNodes"))
set_backend_family(default_truss_bar_family(), rhino, rhino_layer_family("TrussBars"))
