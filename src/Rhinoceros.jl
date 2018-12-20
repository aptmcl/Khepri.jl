export rhino,
       reset_rhino_connection


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

#=
rhino"public void SetView(Point3d position, Point3d target, double lens, bool perspective, string style)"
rhino"public void View(Point3d position, Point3d target, double lens)"
rhino"public void ViewTop()"
rhino"public Point3d ViewCamera()"
rhino"public Point3d ViewTarget()"
rhino"public double ViewLens()"
rhino"public byte Sync()"
rhino"public byte Disconnect()"
rhino"public void Delete(ObjectId id)"
rhino"public void DeleteMany(ObjectId[] ids)"
rhino"public Guid Copy(ObjectId id)"
=#
rhino"public Guid Point(Point3d p)"
rhino"public Point3d PointPosition(Guid ent)"

rhino"public Guid PolyLine(Point3d[] pts)"
rhino"public Point3d[] LineVertices(RhinoObject id)"
rhino"public Guid Spline(Point3d[] pts)"
#=
rhino"public Entity InterpSpline(Point3d[] pts, Vector3d tan0, Vector3d tan1)"
=#
rhino"public Guid ClosedPolyLine(Point3d[] pts)"
rhino"public Guid ClosedSpline(Point3d[] pts)"
#=
rhino"public Entity InterpClosedSpline(Point3d[] pts)"
=#
rhino"public Guid Circle(Point3d c, Vector3d n, double r)"
rhino"public Point3d CircleCenter(RhinoObject obj)"
rhino"public Vector3d CircleNormal(RhinoObject obj)"
rhino"public double CircleRadius(RhinoObject obj)"
rhino"public Guid Ellipse(Point3d c, Vector3d n, double radiusX, double radiusY)"
rhino"public Guid Arc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)"
rhino"public Guid JoinCurves(Guid[] objs)"
#=

rhino"public Entity Text(string str, Point3d corner, Vector3d vx, Vector3d vy, double height)"

rhino"public Guid SurfaceFromCurve(Entity curve)"
=#
rhino"public Guid SurfaceCircle(Point3d c, Vector3d n, double r)"
rhino"public Guid SurfaceEllipse(Point3d c, Vector3d n, double radiusX, double radiusY)"
rhino"public Guid SurfaceArc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)"
rhino"public Guid SurfaceClosedPolyLine(Point3d[] pts)"
#rhino"public Guid[] SurfaceFromCurves(RhinoObject[] objs)"
rhino"public Guid Sphere(Point3d c, double r)"
rhino"public Guid Torus(Point3d c, Vector3d vz, double majorRadius, double minorRadius)"
rhino"public Brep Cylinder(Point3d bottom, double radius, Point3d top)"
rhino"public Brep Cone(Point3d bottom, double radius, Point3d top)"
rhino"public Brep ConeFrustum(Point3d bottom, double bottom_radius, Point3d top, double top_radius)"
rhino"public Brep Box(Point3d corner, Vector3d vx, Vector3d vy, double dx, double dy, double dz)"
rhino"public Brep XYCenteredBox(Point3d corner, Vector3d vx, Vector3d vy, double dx, double dy, double dz)"
#=
rhino"public Guid IrregularPyramidMesh(Point3d[] pts, Point3d apex)"
=#
rhino"public Brep IrregularPyramid(Point3d[] pts, Point3d apex)"
rhino"public Brep IrregularPyramidFrustum(Point3d[] bpts, Point3d[] tpts)"
#=rhino"public Entity MeshFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN)"
rhino"public Entity SurfaceFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level)"
rhino"public Entity SolidFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level, double thickness)"
=#
rhino"public Brep[] Thicken(RhinoObject obj, double thickness)"
#=
rhino"public double[] CurveDomain(Entity ent)"
rhino"public double CurveLength(Entity ent)"
rhino"public Frame3d CurveFrameAt(Entity ent, double t)"
rhino"public Frame3d CurveFrameAtLength(Entity ent, double l)"
rhino"public Guid NurbSurfaceFrom(ObjectId id)"
rhino"public double[] SurfaceDomain(Entity ent)"
rhino"public Frame3d SurfaceFrameAt(Entity ent, double u, double v)"
=#
rhino"public Brep Extrusion(RhinoObject obj, Vector3d dir)"
#=
rhino"public Guid Sweep(ObjectId pathId, ObjectId profileId, double rotation, double scale)"
rhino"public Guid Loft(ObjectId[] profilesIds, ObjectId[] guidesIds, bool ruled, bool closed)"
rhino"public void Unite(ObjectId objId0, ObjectId objId1)"
=#
rhino"public Guid[] Intersect(RhinoObject obj0, RhinoObject obj1)"
rhino"public Guid[] Subtract(RhinoObject obj0, RhinoObject obj1)"
#=
rhino"public Guid Revolve(ObjectId profileId, Point3d p, Vector3d n, double startAngle, double amplitude)"
rhino"public void Move(ObjectId id, Vector3d v)"
rhino"public void Scale(ObjectId id, Point3d p, double s)"
rhino"public void Rotate(ObjectId id, Point3d p, Vector3d n, double a)"
rhino"public Guid Mirror(ObjectId id, Point3d p, Vector3d n, bool copy)"
rhino"public Point3d[] GetPoint(string prompt)"
rhino"public Point3d[] BoundingBox(ObjectId[] ids)"
rhino"public void ZoomExtents()"
rhino"public void SetSystemVariableInt(string name, int value)"
rhino"public int Render(int width, int height, string path)"
rhino"public int Command(string cmd)"
rhino"public void DisableUpdate()"
rhino"public void EnableUpdate()"

=#
rhino"public bool IsPoint(RhinoObject e)"
rhino"public bool IsCircle(RhinoObject e)"
rhino"public bool IsPolyLine(RhinoObject e)"
rhino"public bool IsSpline(RhinoObject e)"
rhino"public bool IsInterpSpline(RhinoObject e)"
rhino"public bool IsClosedPolyLine(RhinoObject e)"
rhino"public bool IsClosedSpline(RhinoObject e)"
rhino"public bool IsInterpClosedSpline(RhinoObject e)"
rhino"public bool IsEllipse(RhinoObject e)"
rhino"public bool IsArc(RhinoObject e)"
rhino"public bool IsText(RhinoObject e)"
rhino"public byte ShapeCode(RhinoObject id)"
#=
rhino"public BIMLevel FindOrCreateLevelAtElevation(double elevation)"
rhino"public BIMLevel UpperLevel(BIMLevel currentLevel, double addedElevation)"
rhino"public double GetLevelElevation(BIMLevel level)"

rhino"public FloorFamily FloorFamilyInstance(double totalThickness, double coatingThickness)"
rhino"public Entity LightweightPolyLine(Point2d[] pts, double[] angles, double elevation)"
rhino"public Entity SurfaceLightweightPolyLine(Point2d[] pts, double[] angles, double elevation)"
rhino"public ObjectId CreatePathFloor(Point2d[] pts, double[] angles, BIMLevel level, FloorFamily family)"
=#
rhino"public Guid ClosedPathCurveArray(Point3d[] pts, double[] angles)"
rhino"public Brep[] PathWall(RhinoObject obj, double thickness, double height)"
rhino"public Brep RectangularTable(Point3d c, double angle, double length, double width, double height, double top_thickness, double leg_thickness)"


#rhino"public TableFamily FindOrCreateTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)"
#rhino"public ChairFamily FindOrCreateChairFamily(double length, double width, double height, double seat_height, double thickness)"
#rhino"public TableChairFamily FindOrCreateTableChairFamily(TableFamily tableFamily, ChairFamily chairFamily, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)"
#rhino"public Brep RectangularTable(Point3d c, double angle, TableFamily family)"
#rhino"public Brep Chair(Point3d c, double angle, ChairFamily family)"
#rhino"public Guid InstanceChair(Point3d c, double angle, ChairFamily family)"
#rhino"public Brep[] RowOfChairs(Point3d c, double angle, int n, double spacing, ChairFamily family)"
#rhino"public Brep[] CenteredRowOfChairs(Point3d c, double angle, int n, double spacing, ChairFamily family)"
#rhino"public Brep[] RectangularTableAndChairs(Point3d c, double angle, TableChairFamily f)"

rhino"public Brep[] BaseRectangularTable(double length, double width, double height, double top_thickness, double leg_thickness)"
rhino"public int CreateRectangularTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)"
rhino"public Guid Table(Point3d c, double angle, int family)"
rhino"public Brep[] BaseChair(double length, double width, double height, double seat_height, double thickness)"
rhino"public int CreateChairFamily(double length, double width, double height, double seat_height, double thickness)"
rhino"public Guid Chair(Point3d c, double angle, int family)"
rhino"public GeometryBase[] BaseRectangularTableAndChairs(int tableFamily, int chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)"
rhino"public int CreateRectangularTableAndChairsFamily(int tableFamily, int chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)"
rhino"public Guid TableAndChairs(Point3d c, double angle, int family)"

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

const rhino = RH(LazyParameter(TCPSocket, create_RH_connection))
reset_rhino_connection() = reset(rhino.connection)

# This should not be done automatically
backend_name(b::RH) = "Rhino"

realize(b::RH, s::EmptyShape) =
  RHEmptyRef()
realize(b::RH, s::UniversalShape) =
  RHUniversalRef()

backend_stroke(b::RH, path::CircularPath) =
    RHCircle(connection(b), path.center, vz(1, path.center.cs), path.radius)
backend_stroke(b::RH, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        RHClosedPolyLine(connection(b), [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)])
    end
backend_stroke(b::RH, path::OpenPolygonalPath) =
  	RHPolyLine(connection(b), path.vertices)
backend_stroke(b::RH, path::ClosedPolygonalPath) =
    RHClosedPolyLine(connection(b), path.vertices)

backend_stroke_unite(b::RH, refs) = RHJoinCurves(connection(b), refs)

#=backend_fill(b::RH, path::ClosedPolygonalPath) =
    RHSurfaceClosedPolyLine(connection(b), path.vertices)
    backend_fill(b::RH, path::RectangularPath) =
        let c = path.corner,
            dx = path.dx,
            dy = path.dy
            RHSurfaceClosedPolyLine(connection(b), [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)])
        end
backend_fill(b::RH, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        SurfaceClosedPolyLine(connection(b), [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)])
    end
=#

backend_fill_curves(b::RH, gs::Guids) = RHSurfaceFromCurves(connection(b), gs)
backend_fill_curves(b::RH, g::Guid) = RHSurfaceFromCurves(connection(b), [g])

backend_stroke_line(b::RH, vs) = RHPolyLine(connection(b), vs)

backend_stroke_arc(b::RH, center::Loc, radius::Real, start_angle::Real, amplitude::Real) =
    let end_angle = start_angle + amplitude
        if end_angle > start_angle
            RHArc(connection(b), center, vz(1, center.cs), radius, start_angle, end_angle)
        else
            RHArc(connection(b), center, vz(1, center.cs), radius, end_angle, start_angle)
        end
    end

realize(b::RH, s::Point) =
  RHPoint(connection(b), s.position)

realize(b::RH, s::Line) =
  RHPolyLine(connection(b), s.vertices)

#=
realize(b::RH, s::Spline) =
  if (s.v0 == false) && (s.v1 == false)
    RHSpline(connection(b), s.points)
  elseif (s.v0 != false) && (s.v1 != false)
    RHInterpSpline(connection(b), s.points, s.v0, s.v1)
  else
    RHInterpSpline(connection(b),
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1)
  end

realize(b::RH, s::ClosedSpline) =
  RHInterpClosedSpline(connection(b), s.points)
=#

realize(b::RH, s::Circle) =
  RHCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)

realize(b::RH, s::Arc) =
  if s.radius == 0
    RHPoint(connection(b), s.center)
  elseif s.amplitude == 0
    RHPoint(connection(b), s.center + vpol(s.radius, s.start_angle, s.center.cs))
  elseif abs(s.amplitude) >= 2*pi
    RHCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
  else
    end_angle = s.start_angle + s.amplitude
    if end_angle > s.start_angle
      RHArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle)
    else
      RHArc(connection(b), s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle)
    end
  end

realize(b::RH, s::Ellipse) =
    RHEllipse(connection(b), s.center, vz(1, s.center.cs), s.radius_x, s.radius_y)

realize(b::RH, s::EllipticArc) =
  error("Finish this")

realize(b::RH, s::Polygon) =
  RHClosedPolyLine(connection(b), s.vertices)

realize(b::RH, s::RegularPolygon) =
  RHClosedPolyLine(connection(b), regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))

realize(b::RH, s::Rectangle) =
  RHClosedPolyLine(connection(b), [s.c, add_x(s.c, s.dx), add_xy(s.c, s.dx, s.dy), add_y(s.c, s.dy)])

rhino"public Guid SurfaceCircle(Point3d c, Vector3d n, double r)"

realize(b::RH, s::SurfaceCircle) =
  RHSurfaceCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)

realize(b::RH, s::SurfaceArc) =
  RHSurfaceArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, s.start_angle + s.amplitude)

realize(b::RH, s::SurfaceElliptic_Arc) =
  error("Finish this")

realize(b::RH, s::SurfaceEllipse) =
    RHSurfaceEllipse(connection(b), s.center, vz(1, s.center.cs), s.radius_x, s.radius_y)

rhino"public Guid SurfaceClosedPolyLine(Point3d[] pts)"

realize(b::RH, s::SurfacePolygon) =
  RHSurfaceClosedPolyLine(connection(b), s.vertices)
realize(b::RH, s::SurfaceRegularPolygon) =
  RHSurfaceClosedPolyLine(connection(b), regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
realize(b::RH, s::SurfaceRectangle) =
  let c = s.corner
    RHSurfaceClosedPolyLine(connection(b), [c, add_x(c, s.dx), add_xy(c, s.dx, s.dy), add_y(c, s.dy)])
  end
#=
realize(b::RH, s::Surface) =
  let ids = RHSurfaceFromCurves(connection(b), collect_ref(s.frontier))
    foreach(mark_deleted, s.frontier)
    ids
  end
realize(b::RH, s::Text) =
  RHText(connection(b), s.str, s.c, vx(1, s.c.cs), vy(1, s.c.cs), s.h)
=#

rhino"public double[] SurfaceDomain(RhinoObject obj)"
rhino"public Plane SurfaceFrameAt(RhinoObject obj, double u, double v)"

backend_surface_domain(b::RH, s::Shape2D) =
    tuple(RHSurfaceDomain(connection(b), ref(s).value)...)

backend_map_division(b::RH, f::Function, s::Shape2D, nu::Int, nv::Int) =
    let conn = connection(b)
        r = ref(s).value
        (u1, u2, v1, v2) = RHSurfaceDomain(conn, r)
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(RHSurfaceFrameAt(conn, r, u, v))
            end
        end
    end

realize(b::RH, s::Sphere) =
    RHSphere(connection(b), s.center, s.radius)
realize(b::RH, s::Torus) =
    RHTorus(connection(b), s.center, vz(1, s.center.cs), s.re, s.ri)
realize(b::RH, s::Cuboid) =
    RHIrregularPyramidFrustum(connection(b), [s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3])
realize(b::RH, s::RegularPyramidFrustum) =
    RHIrregularPyramidFrustum(connection(b),
                              regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                              regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed))
realize(b::RH, s::RegularPyramid) =
    RHIrregularPyramid(connection(b),
                       regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                       add_z(s.cb, s.h))
realize(b::RH, s::IrregularPyramid) =
    RHIrregularPyramid(connection(b), s.cbs, s.ct)
realize(b::RH, s::RegularPrism) =
    let cbs = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
        RHIrregularPyramidFrustum(connection(b),
                                  cbs,
                                  map(p -> add_z(p, s.h), cbs))
    end
realize(b::RH, s::IrregularPrism) =
    RHIrregularPyramidFrustum(connection(b),
                              s.cbs,
                              map(p -> (p + s.v), s.cbs))

realize(b::RH, s::RightCuboid) =
    RHXYCenteredBox(connection(b), s.cb, vx(1, s.cb.cs), vy(1, s.cb.cs), s.width, s.height, s.h)

realize(b::RH, s::Box) =
    RHBox(connection(b), s.c, vx(1,s.c.cs), vy(1,s.c.cs), s.dx, s.dy, s.dz)
#=
realize(b::RH, s::Cone) =
  RHCone(connection(b), add_z(s.cb, s.h), s.r, s.cb)
realize(b::RH, s::ConeFrustum) =
  RHConeFrustum(connection(b), s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt)
=#

realize(b::RH, s::Cylinder) =
  RHCylinder(connection(b), s.cb, s.r, s.cb + vz(s.h, s.cb.cs))

#realize(b::RH, s::Circle) = RHCircle(connection(b),

backend_extrusion(b::RH, s::Shape, v::Vec) =
    and_mark_deleted(
        map_ref(s) do r
            RHExtrusion(connection(b), r, v)
        end,
        s)
#=
realize(b::Backend, s::Sweep) =
  backend_sweep(b, s.path, s.profile, s.rotation, s.scale)

backend_sweep(b::RH, path::Shape, profile::Shape, rotation::Real, scale::Real) =
  map_ref(profile) do profile_r
    map_ref(path) do path_r
      RHSweep(connection(b), path_r, profile_r, rotation, scale)
    end
  end

realize(b::RH, s::Revolve) =
  and_delete_shape(
    map_ref(s.profile) do r
      RHRevolve(connection(b), r, s.p, s.n, s.start_angle, s.amplitude)
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
  (RHUnite(connection(b), r0.value, r1.value); r0)

intersect_ref(b::RH, r0::RHNativeRef, r1::RHNativeRef) =
    let refs = RHIntersect(connection(b), r0.value, r1.value)
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
    let refs = RHSubtract(connection(b), r0.value, r1.value)
        n = length(refs)
        if n == 0
            RHEmptyRef()
        elseif n == 1
            if is_empty_guid(refs[1]) # failed
                RHSubtractionRef(r0, tuple(r1))
            else
                RHNativeRef(refs[1])
            end
        else
            RHUnionRef(map(RHNativeRef, tuple(refs...)))
        end
    end

subtract_ref(b::RH, r0::RHRef, r1::RHUnionRef) =
  foldl((r0,r1)->subtract_ref(b,r0,r1), r1.values, init=r0)

#=
slice_ref(b::RH, r::RHNativeRef, p::Loc, v::Vec) =
  (RHSlice(connection(b), r.value, p, v); r)

slice_ref(b::RH, r::RHUnionRef, p::Loc, v::Vec) =
  map(r->slice_ref(b, r, p, v), r.values)

=#

#=
realize(b::RH, s::IntersectionShape) =
  foldl((r0,r1)->intersect_ref(b,r0,r1), UniversalRef{RHId}(), map(ref, s.shapes))

realize(b::RH, s::Slice) =
  slice_ref(b, ref(s.shape), s.p, s.n)

realize(b::RH, s::Move) =
  let r = map_ref(s.shape) do r
            RHMove(connection(b), r, s.v)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::RH, s::Scale) =
  let r = map_ref(s.shape) do r
            RHScale(connection(b), r, s.p, s.s)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::RH, s::Rotate) =
  let r = map_ref(s.shape) do r
            RHRotate(connection(b), r, s.p, s.v, s.angle)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::RH, s::Mirror) =
  and_delete_shape(map_ref(s.shape) do r
                    RHMirror(connection(b), r, s.p, s.n, false)
                   end,
                   s.shape)

realize(b::RH, s::UnionMirror) =
  let r0 = ref(s.shape),
      r1 = map_ref(r0) do r
            RHMirror(connection(b), r, s.p, s.n, true)
          end
    UnionRef((r0,r1))
  end

realize(b::RH, s::SurfaceGrid) =
  RHSurfaceFromGrid(connection(b), length(s.ptss), length(s.ptss[1]), vcat(s.ptss...), s.closed_u, s.closed_v, 2)

realize(b::RH, s::Thicken) =
  and_delete_shape(
    map_ref(s.shape) do r
      RHThicken(connection(b), r, s.thickness)
    end,
    s.shape)

=#
# BIM

# Families

backend_get_family(b::RH, f::TableFamily) =
    RHCreateRectangularTableFamily(connection(b), f.length, f.width, f.height, f.top_thickness, f.leg_thickness)
backend_get_family(b::RH, f::ChairFamily) =
    RHCreateChairFamily(connection(b), f.length, f.width, f.height, f.seat_height, f.thickness)
backend_get_family(b::RH, f::TableChairFamily) =
    RHCreateRectangularTableAndChairsFamily(connection(b),
        ref(f.table_family), ref(f.chair_family),
        f.table_family.length, f.table_family.width,
        f.chairs_top, f.chairs_bottom, f.chairs_right, f.chairs_left,
        f.spacing)

backend_rectangular_table(b::RH, c, angle, family) =
    RHTable(connection(b), c, angle, ref(family))

backend_chair(b::RH, c, angle, family) =
    RHChair(connection(b), c, angle, ref(family))

backend_rectangular_table_and_chairs(b::RH, c, angle, family) =
    RHTableAndChairs(connection(b), c, angle, ref(family))

backend_slab(b::RH, profile, thickness) =
    map_ref(b,
            r->RHExtrusion(connection(b), r, thickness),
            ensure_ref(b, backend_fill(b, profile)))

#Beams are aligned along the top axis.
realize(b::RH, s::Beam) =
    let o = loc_from_o_phi(s.cb, s.angle)
        RHXYCenteredBox(connection(b), add_y(o, -s.family.height/2), vx(1, o.cs), vy(1, o.cs), s.family.width, s.family.height, s.h)
    end

#Columns are aligned along the center axis.
realize(b::RH, s::Column) =
    let o = loc_from_o_phi(s.cb, s.angle)
        RHXYCenteredBox(connection(b), o, vx(1, o.cs), vy(1, o.cs), s.family.width, s.family.height, s.h)
    end

backend_wall(b::RH, path, height, thickness) =
    RHPathWall(connection(b), backend_stroke(b, path), thickness, height)

############################################

backend_bounding_box(b::RH, shapes::Shapes) =
  RHBoundingBox(connection(b), collect_ref(shapes))

import Base.view
set_view(camera::XYZ, target::XYZ, lens::Real, b::RH) =
  RHView(connection(b), camera, target, lens)

get_view(b::RH) =
  let c = connection(b)
    RHViewCamera(c), RHViewTarget(c), RHViewLens(c)
  end

zoom_extents(b::RH) =
  RHZoomExtents(connection(b))

view_top(b::RH) =
    RHViewTop(connection(b))

#
rhino"public int DeleteAll()"
rhino"public int DeleteAllInLayer(String name)"
rhino"public void Delete(Guid id)"
rhino"public void DeleteMany(Guid[] ids)"

backend_delete_shapes(b::RH, shapes::Shapes) =
  RHDeleteMany(connection(b), collect_ref(shapes))

delete_all_shapes(b::RH) =
  RHDeleteAll(connection(b))

# Layers
rhino"public String CreateLayer(string name)"
rhino"public void SetLayerColor(ObjectId id, byte r, byte g, byte b)"
rhino"public void SetShapeColor(ObjectId id, byte r, byte g, byte b)"
rhino"public String CurrentLayer()"
rhino"public void SetCurrentLayer(String id)"
rhino"public String ShapeLayer(RhinoObject objId)"
rhino"public void SetShapeLayer(RhinoObject objId, String layerId)"

RHLayer = String

current_layer(b::RH)::RHLayer =
  RHCurrentLayer(connection(b))

current_layer(layer::RHLayer, b::RH) =
  RHSetCurrentLayer(connection(b), layer)

create_layer(name::String, b::RH) =
  RHCreateLayer(connection(b), name)

create_layer(name::String, color::RGB, b::RH) =
  let layer = RHCreateLayer(connection(b), name)
    RHSetLayerColor(connection(b), layer, color.r, color.g, color.b)
    layer
  end

delete_all_shapes_in_layer(layer::RHLayer, b::RH) =
  RHDeleteAllInLayer(connection(b), layer)

shape_from_ref(r, b::RH) =
  let c = connection(b)
    let code = RHShapeCode(c, r)
        ref = LazyRef(b, RHNativeRef(RHClone(connection(b), r))) # HACK CLONING!!!!!
      if code == 1
        point(RHPointPosition(c, r),
              backend=b, ref=ref)
      elseif code == 2
        circle(maybe_loc_from_o_vz(RHCircleCenter(c, r), RHCircleNormal(c, r)),
               RHCircleRadius(c, r),
               backend=b, ref=ref)
      elseif 3 <= code <= 6
        line(RHLineVertices(c, r),
             backend=b, ref=ref)
      elseif code == 7
        spline([xy(0,0)], false, false, #HACK obtain interpolation points
               backend=b, ref=ref)
      elseif 103 <= code <= 106
        polygon(RHLineVertices(c, r),
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

"""            obj = _geometry_from_id(r)
            ref = native_ref(r)
            if isinstance(obj, geo.Point):
                return point.new_ref(ref, fromPt(obj.Location))
            elif isinstance(obj, geo.Curve):
                if rh.IsLine(r) or rh.IsPolyline(r):
                    if rh.IsCurveClosed(r):
                        return polygon.new_ref(ref, [fromPt(p) for p in rh.CurvePoints(r)[:-1]])
                    else:
                        return line.new_ref(ref, [fromPt(p) for p in rh.CurvePoints(r)])
                elif obj.IsCircle(Rhino.RhinoMath.ZeroTolerance):
                    return circle.new_ref(ref, fromPt(rh.CircleCenterPoint(r)), rh.CircleRadius(r))
                elif rh.IsCurveClosed(r):
                    return closed_spline.new_ref(ref, [fromPt(p) for p in rh.CurvePoints(r)])
                else:
                    return spline.new_ref(ref, [fromPt(p) for p in rh.CurvePoints(r)])
            elif rh.IsObject(r) and rh.IsObjectSolid(r):
                return solid(native_ref(r))
            elif rh.IsSurface(r) or rh.IsPolysurface(r):
                return surface(native_ref(r))
            else:
                raise RuntimeError("{0}: Unknown Rhino object {1}".format('shape_from_ref', r))
                """

#
rhino"public Point3d[] GetPosition(string prompt)"

select_position(prompt::String, b::RH) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = RHGetPosition(connection(b), prompt)
      length(ans) > 0 ? ans[1] : nothing
    end
  end

select_positions(prompt::String, b::RH) =
  let ps = Loc[]
      p = nothing
    @info "$(prompt) on the $(b) backend."
    while ((p = RHGetPosition(connection(b), prompt)) |> length) > 0
        push!(ps, p...)
    end
    ps
  end

rhino"public Guid[] GetPoint(string prompt)"
rhino"public Guid[] GetPoints(string prompt)"

# HACK: The next operations should receive a set of shapes to avoid re-creating already existing shapes

select_point(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, RHGetPoint)

select_points(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, RHGetPoints)

rhino"public Guid[] GetCurve(string prompt)"
rhino"public Guid[] GetCurves(string prompt)"

select_curve(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, RHGetCurve)

select_curves(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, RHGetCurves)

rhino"public Guid[] GetSurface(string prompt)"
rhino"public Guid[] GetSurfaces(string prompt)"

select_surface(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, RHGetSurface)

select_surfaces(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, RHGetSurfaces)

rhino"public Guid[] GetSolid(string prompt)"
rhino"public Guid[] GetSolids(string prompt)"

select_solid(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, RHGetSolid)

select_solids(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, RHGetSolids)

rhino"public Guid[] GetShape(string prompt)"
rhino"public Guid[] GetShapes(string prompt)"

select_shape(prompt::String, b::RH) =
  select_one_with_prompt(prompt, b, RHGetShape)

select_shapes(prompt::String, b::RH) =
  select_many_with_prompt(prompt, b, RHGetShapes)

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
    print("captured_shapes(autocad, [")
    for s in ss
      print(ref(s).value)
      print(", ")
    end
    println("])")
  end


rhino"public Guid[] GetAllShapes()"
rhino"public Guid[] GetAllShapesInLayer(String name)"

all_shapes(b::RH) =
  [shape_from_ref(r) for r in RHGetAllShapes(connection(b))]

all_shapes_in_layer(layer, b::RH) =
  [shape_from_ref(r) for r in RHGetAllShapesInLayer(connection(b), layer)]
