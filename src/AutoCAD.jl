export autocad

#=
We need to ensure the AutoCAD plugin is properly installed.
For AutoCAD, there are a few places where plugins can be installed:

A plug-in can be deployed by placing it in one of the ApplicationPlugins or ApplicationAddins folders on a local drive.

General Installation folder
%PROGRAMFILES%\Autodesk\ApplicationPlugins
All Users Profile folders
%ALLUSERSPROFILE%\Autodesk\ApplicationPlugins
User Profile folders
%APPDATA%\Autodesk\ApplicationPlugins

When the SECURELOAD system variable is set to 1 or 2,
the program is restricted to loading and executing files that contain code from trusted locations;
trusted locations are specified by the TRUSTEDPATHS system variable.

=#
dlls = ["KhepriBase.dll", "KhepriAutoCAD.dll"]
bundle_name = "Khepri.bundle"
bundle_dll_folder = joinpath(bundle_name, "Contents")
xml_name = "PackageContents.xml"
bundle_xml = joinpath(bundle_name, xml_name)
local_plugins = joinpath(dirname(dirname(abspath(@__FILE__))), "Plugins", "AutoCAD")
local_khepri_plugin = joinpath(local_plugins, bundle_name)
local_khepri_plugin_dll_folder = joinpath(local_plugins, bundle_dll_folder)

autocad_general_plugins = joinpath(dirname(joinpath(ENV["CommonProgramFiles"])), "Autodesk", "ApplicationPlugins")
autocad_allusers_plugins = joinpath(ENV["ALLUSERSPROFILE"], "Autodesk", "ApplicationPlugins")
autocad_user_plugins = joinpath(ENV["APPDATA"], "Autodesk", "ApplicationPlugins")
autocad_khepri_plugin = joinpath(autocad_user_plugins, bundle_name)
autocad_khepri_plugin_dll_folder = joinpath(autocad_user_plugins, bundle_dll_folder)

update_file_if_needed(path) =
    let local_path = joinpath(local_khepri_plugin, path)
        autocad_path = joinpath(autocad_khepri_plugin, path)
        if ! isfile(autocad_path) || mtime(autocad_path) < mtime(local_path)
            isfile(autocad_path) && rm(autocad_path)
            cp(local_path, autocad_path)
        end
    end

update_plugin() =
    begin
        # Do we have the bundle folder?
        isdir(autocad_khepri_plugin) || mkpath(autocad_khepri_plugin)
        update_file_if_needed(xml_name)
        isdir(autocad_khepri_plugin_dll_folder) || mkpath(autocad_khepri_plugin_dll_folder)
        map(dlls) do dll
            update_file_if_needed(joinpath("Contents", dll))
        end
    end

checked_plugin = false

check_plugin() =
    begin
        global checked_plugin
        if ! checked_plugin
            @info("Checking plugin...")
            try
                update_plugin()
                @info("done.")
                checked_plugin = true
            catch exc
                if isa(exc, Base.IOError)
                    @error("Please, close AutoCAD and retry.")
                else
                    throw(exc)
                end
            end
        end
    end


#app = AutoCAD()
#doc = app.ActiveDocument
#doc.SendCommand("(command "._NETLOAD" "{0}") ".format(join(dirname(dirname(abspath(__file__))),
#                                                      "Khepri", "KhepriAutoCAD", "KhepriAutoCAD", "bin", "x64", "Debug", "KhepriAutoCAD.dll")).replace("\\","/"))
#db = doc.ModelSpace
#util = doc.Utility

macro acad_str(str)
    rpc("ACAD", str)
end

# We need some additional Encoders
encode_Entity = encode_int
decode_Entity = decode_int_or_error
encode_ObjectId = encode_int
decode_ObjectId = decode_int_or_error
encode_ObjectId_array = encode_int_array
decode_ObjectId_array = decode_int_array

acad"public int DeleteAll()"
acad"public void SetLengthUnit(String unit)"
acad"public void SetView(Point3d position, Point3d target, double lens, bool perspective, string style)"
acad"public void View(Point3d position, Point3d target, double lens)"
acad"public void ViewTop()"
acad"public Point3d ViewCamera()"
acad"public Point3d ViewTarget()"
acad"public double ViewLens()"
acad"public byte Sync()"
acad"public byte Disconnect()"
acad"public void Delete(ObjectId id)"
acad"public void DeleteMany(ObjectId[] ids)"
acad"public ObjectId Copy(ObjectId id)"

acad"public Entity Point(Point3d p)"
acad"public Point3d PointPosition(Entity ent)"

acad"public Entity PolyLine(Point3d[] pts)"
acad"public Point3d[] LineVertices(ObjectId id)"

acad"public Entity Spline(Point3d[] pts)"
acad"public Entity InterpSpline(Point3d[] pts, Vector3d tan0, Vector3d tan1)"
acad"public Entity ClosedPolyLine(Point3d[] pts)"
acad"public Entity ClosedSpline(Point3d[] pts)"
acad"public Entity InterpClosedSpline(Point3d[] pts)"

acad"public Entity Circle(Point3d c, Vector3d n, double r)"
acad"public Point3d CircleCenter(Entity ent)"
acad"public Vector3d CircleNormal(Entity ent)"
acad"public double CircleRadius(Entity ent)"

acad"public Entity Ellipse(Point3d c, Vector3d n, Vector3d majorAxis, double radiusRatio)"

acad"public Entity Arc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)"
acad"public Point3d ArcCenter(Entity ent)"
acad"public Vector3d ArcNormal(Entity ent)"
acad"public double ArcRadius(Entity ent)"
acad"public double ArcStartAngle(Entity ent)"
acad"public double ArcEndAngle(Entity ent)"

acad"public ObjectId JoinCurves(ObjectId[] ids)"

acad"public Entity Text(string str, Point3d corner, Vector3d vx, Vector3d vy, double height)"
acad"public Entity SurfaceFromCurve(Entity curve)"
acad"public Entity SurfaceCircle(Point3d c, Vector3d n, double r)"
acad"public Entity SurfaceEllipse(Point3d c, Vector3d n, Vector3d majorAxis, double radiusRatio)"
acad"public Entity SurfaceArc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)"
acad"public Entity SurfaceClosedPolyLine(Point3d[] pts)"
acad"public ObjectId[] SurfaceFromCurves(ObjectId[] ids)"
acad"public ObjectId[] CurvesFromSurface(ObjectId id)"
acad"public Entity Sphere(Point3d c, double r)"
acad"public Entity Torus(Point3d c, Vector3d vz, double majorRadius, double minorRadius)"
acad"public Entity ConeFrustum(Point3d bottom, double base_radius, Point3d top, double top_radius)"
acad"public Entity Cylinder(Point3d bottom, double radius, Point3d top)"
acad"public Entity Cone(Point3d bottom, double radius, Point3d top)"
acad"public Entity Box(Frame3d frame, double dx, double dy, double dz)"
acad"public Entity CenteredBox(Frame3d frame, double dx, double dy, double dz)"
acad"public ObjectId IrregularPyramidMesh(Point3d[] pts, Point3d apex)"
acad"public ObjectId IrregularPyramid(Point3d[] pts, Point3d apex)"
acad"public ObjectId IrregularPyramidFrustum(Point3d[] bpts, Point3d[] tpts)"
acad"public Entity MeshFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN)"
acad"public Entity SurfaceFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level)"
acad"public Entity SolidFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level, double thickness)"
acad"public ObjectId Thicken(ObjectId obj, double thickness)"
acad"public ObjectId NurbSurfaceFrom(ObjectId id)"
acad"public ObjectId Extrude(ObjectId profileId, Vector3d dir)"
acad"public ObjectId Sweep(ObjectId pathId, ObjectId profileId, double rotation, double scale)"
acad"public ObjectId Loft(ObjectId[] profilesIds, ObjectId[] guidesIds, bool ruled, bool closed)"
acad"public ObjectId Unite(ObjectId objId0, ObjectId objId1)"
acad"public ObjectId Intersect(ObjectId objId0, ObjectId objId1)"
acad"public ObjectId Subtract(ObjectId objId0, ObjectId objId1)"
acad"public void Slice(ObjectId id, Point3d p, Vector3d n)"
acad"public ObjectId Revolve(ObjectId profileId, Point3d p, Vector3d n, double startAngle, double amplitude)"
acad"public void Transform(ObjectId id, Frame3d frame)"
acad"public void Move(ObjectId id, Vector3d v)"
acad"public void Scale(ObjectId id, Point3d p, double s)"
acad"public void Rotate(ObjectId id, Point3d p, Vector3d n, double a)"
acad"public ObjectId Mirror(ObjectId id, Point3d p, Vector3d n, bool copy)"
acad"public Point3d[] GetPoint(string prompt)"
acad"public ObjectId[] GetAllShapes()"
acad"public ObjectId[] GetAllShapesInLayer(ObjectId layerId)"
acad"public Point3d[] BoundingBox(ObjectId[] ids)"
acad"public void ZoomExtents()"
acad"public ObjectId CreateLayer(string name)"
acad"public void SetLayerColor(ObjectId id, byte r, byte g, byte b)"
acad"public void SetShapeColor(ObjectId id, byte r, byte g, byte b)"
acad"public ObjectId CurrentLayer()"
acad"public void SetCurrentLayer(ObjectId id)"
acad"public ObjectId ShapeLayer(ObjectId objId)"
acad"public void SetShapeLayer(ObjectId objId, ObjectId layerId)"
acad"public void SetSystemVariableInt(string name, int value)"
acad"public int Command(string cmd)"
acad"public void DisableUpdate()"
acad"public void EnableUpdate()"

acad"public bool IsPoint(Entity e)"
acad"public bool IsCircle(Entity e)"
acad"public bool IsPolyLine(Entity e)"
acad"public bool IsSpline(Entity e)"
acad"public bool IsInterpSpline(Entity e)"
acad"public bool IsClosedPolyLine(Entity e)"
acad"public bool IsClosedSpline(Entity e)"
acad"public bool IsInterpClosedSpline(Entity e)"
acad"public bool IsEllipse(Entity e)"
acad"public bool IsArc(Entity e)"
acad"public bool IsText(Entity e)"
acad"public byte ShapeCode(ObjectId id)"

acad"public BIMLevel FindOrCreateLevelAtElevation(double elevation)"
acad"public BIMLevel UpperLevel(BIMLevel currentLevel, double addedElevation)"
acad"public double GetLevelElevation(BIMLevel level)"

acad"public FloorFamily FloorFamilyInstance(double totalThickness, double coatingThickness)"
acad"public Entity LightweightPolyLine(Point2d[] pts, double[] angles, double elevation)"
acad"public Entity SurfaceLightweightPolyLine(Point2d[] pts, double[] angles, double elevation)"
acad"public ObjectId CreatePathFloor(Point2d[] pts, double[] angles, BIMLevel level, FloorFamily family)"

acad"public ObjectId CreateBlockFromShapes(String baseName, ObjectId[] ids)"
acad"public ObjectId CreateBlockInstance(ObjectId id, Frame3d frame)"
acad"public ObjectId CreateInstanceFromBlockNamed(String name, Frame3d frame)"
acad"public ObjectId CreateInstanceFromBlockNamedAtRotated(String name, Point3d c, double angle)"
acad"public ObjectId CreateRectangularTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)"
acad"public ObjectId Table(Point3d c, double angle, ObjectId family)"
acad"public ObjectId CreateChairFamily(double length, double width, double height, double seat_height, double thickness)"
acad"public ObjectId Chair(Point3d c, double angle, ObjectId family)"
acad"public ObjectId CreateRectangularTableAndChairsFamily(ObjectId tableFamily, ObjectId chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)"
acad"public ObjectId TableAndChairs(Point3d c, double angle, ObjectId family)"
acad"public ObjectId CreateAlignedDimension(Point3d p0, Point3d p1, Point3d p, double scale, String mark)"
acad"public String TextString(Entity ent)"
acad"public Point3d TextPosition(Entity ent)"
acad"public double TextHeight(Entity ent)"
acad"public String MTextString(Entity ent)"
acad"public Point3d MTextPosition(Entity ent)"
acad"public double MTextHeight(Entity ent)"

acad"public void SaveAs(String pathname, String format)"


abstract type ACADKey end
const ACADId = Int
const ACADIds = Vector{ACADId}
const ACADRef = GenericRef{ACADKey, ACADId}
const ACADRefs = Vector{ACADRef}
const ACADEmptyRef = EmptyRef{ACADKey, ACADId}
const ACADUniversalRef = UniversalRef{ACADKey, ACADId}
const ACADNativeRef = NativeRef{ACADKey, ACADId}
const ACADUnionRef = UnionRef{ACADKey, ACADId}
const ACADSubtractionRef = SubtractionRef{ACADKey, ACADId}
const ACAD = Socket_Backend{ACADKey, ACADId}

void_ref(b::ACAD) = ACADNativeRef(-1)

create_ACAD_connection() =
    begin
        check_plugin()
        create_backend_connection("AutoCAD", 11000)
    end

const autocad = ACAD(LazyParameter(TCPSocket, create_ACAD_connection))

#current_backend(autocad)


backend_stroke_color(b::ACAD, path::Path, color::RGB) =
    let r = backend_stroke(b, path)
        ACADSetShapeColor(connection(b), r, color.r, color.g, color.b)
        r
    end

backend_stroke(b::ACAD, path::CircularPath) =
    ACADCircle(connection(b), path.center, vz(1, path.center.cs), path.radius)
backend_stroke(b::ACAD, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        ACADClosedPolyLine(connection(b), [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)])
    end
backend_stroke(b::ACAD, path::ArcPath) =
    backend_stroke_arc(b, path.center, path.radius, path.start_angle, path.amplitude)

backend_stroke(b::ACAD, path::OpenPolygonalPath) =
  	ACADPolyLine(connection(b), path.vertices)
backend_stroke(b::ACAD, path::ClosedPolygonalPath) =
    ACADClosedPolyLine(connection(b), path.vertices)
backend_fill(b::ACAD, path::ClosedPolygonalPath) =
    ACADSurfaceClosedPolyLine(connection(b), path.vertices)
backend_fill(b::ACAD, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        ACADSurfaceClosedPolyLine(connection(b), [c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)])
    end
backend_fill_curves(b::ACAD, refs::ACADIds) = ACADSurfaceFromCurves(connection(b), refs)
backend_fill_curves(b::ACAD, ref::ACADId) = ACADSurfaceFromCurves(connection(b), [ref])

backend_stroke_line(b::ACAD, vs) = ACADPolyLine(connection(b), vs)

backend_stroke_arc(b::ACAD, center::Loc, radius::Real, start_angle::Real, amplitude::Real) =
    let end_angle = start_angle + amplitude
        if end_angle > start_angle
            ACADArc(connection(b), center, vz(1, center.cs), radius, start_angle, end_angle)
        else
            ACADArc(connection(b), center, vz(1, center.cs), radius, end_angle, start_angle)
        end
    end
backend_stroke_unite(b::ACAD, refs) = ACADJoinCurves(connection(b), refs)



realize(b::ACAD, s::EmptyShape) =
  ACADEmptyRef()
realize(b::ACAD, s::UniversalShape) =
  ACADUniversalRef()
realize(b::ACAD, s::Point) =
  ACADPoint(connection(b), s.position)
realize(b::ACAD, s::Line) =
  ACADPolyLine(connection(b), s.vertices)
realize(b::ACAD, s::Spline) =
  if (s.v0 == false) && (s.v1 == false)
    #ACADSpline(connection(b), s.points)
    ACADInterpSpline(connection(b),
                     s.points,
                     s.points[2]-s.points[1],
                     s.points[end]-s.points[end-1])
  elseif (s.v0 != false) && (s.v1 != false)
    ACADInterpSpline(connection(b), s.points, s.v0, s.v1)
  else
    ACADInterpSpline(connection(b),
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1)
  end
realize(b::ACAD, s::ClosedSpline) =
  ACADInterpClosedSpline(connection(b), s.points)
realize(b::ACAD, s::Circle) =
  ACADCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
realize(b::ACAD, s::Arc) =
  if s.radius == 0
    ACADPoint(connection(b), s.center)
  elseif s.amplitude == 0
    ACADPoint(connection(b), s.center + vpol(s.radius, s.start_angle, s.center.cs))
  elseif abs(s.amplitude) >= 2*pi
    ACADCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
  else
    end_angle = s.start_angle + s.amplitude
    if end_angle > s.start_angle
      ACADArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle)
    else
      ACADArc(connection(b), s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle)
    end
  end

realize(b::ACAD, s::Ellipse) =
  if s.radius_x > s.radius_y
    ACADEllipse(connection(b), s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x)
  else
    ACADEllipse(connection(b), s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y)
  end
realize(b::ACAD, s::EllipticArc) =
  error("Finish this")

realize(b::ACAD, s::Polygon) =
  ACADClosedPolyLine(connection(b), s.vertices)
realize(b::ACAD, s::RegularPolygon) =
  ACADClosedPolyLine(connection(b), regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
realize(b::ACAD, s::Rectangle) =
  ACADClosedPolyLine(connection(b), [s.c, add_x(s.c, s.dx), add_xy(s.c, s.dx, s.dy), add_y(s.c, s.dy)])
realize(b::ACAD, s::SurfaceCircle) =
  ACADSurfaceCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
realize(b::ACAD, s::SurfaceArc) =
    #ACADSurfaceArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, s.start_angle + s.amplitude)
    if s.radius == 0
        ACADPoint(connection(b), s.center)
    elseif s.amplitude == 0
        ACADPoint(connection(b), s.center + vpol(s.radius, s.start_angle, s.center.cs))
    elseif abs(s.amplitude) >= 2*pi
        ACADSurfaceCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
    else
        end_angle = s.start_angle + s.amplitude
        if end_angle > s.start_angle
            ACADSurfaceFromCurves(connection(b),
                [ACADArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle),
                 ACADPolyLine(connection(b), [add_pol(s.center, s.radius, end_angle),
                                              add_pol(s.center, s.radius, s.start_angle)])])
        else
            ACADSurfaceFromCurves(connection(b),
                [ACADArc(connection(b), s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle),
                 ACADPolyLine(connection(b), [add_pol(s.center, s.radius, s.start_angle),
                                              add_pol(s.center, s.radius, end_angle)])])
        end
    end

#realize(b::ACAD, s::SurfaceElliptic_Arc) = ACADCircle(connection(b),
#realize(b::ACAD, s::SurfaceEllipse) = ACADCircle(connection(b),
realize(b::ACAD, s::SurfacePolygon) =
  ACADSurfaceClosedPolyLine(connection(b), s.vertices)
realize(b::ACAD, s::SurfaceRegularPolygon) =
  ACADSurfaceClosedPolyLine(connection(b), regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
realize(b::ACAD, s::SurfaceRectangle) =
  ACADSurfaceClosedPolyLine(connection(b), [s.c, add_x(s.c, s.dx), add_xy(s.c, s.dx, s.dy), add_y(s.c, s.dy)])
realize(b::ACAD, s::Surface) =
  let #ids = map(r->ACADNurbSurfaceFrom(connection(b),r), ACADSurfaceFromCurves(connection(b), collect_ref(s.frontier)))
      ids = ACADSurfaceFromCurves(connection(b), collect_ref(s.frontier))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::ACAD, s::Shape2D) =
    map(shape_from_ref, ACADCurvesFromSurface(connection(b), ref(s).value))

# Iterating over curves and surfaces

acad"public double[] CurveDomain(Entity ent)"
acad"public double CurveLength(Entity ent)"
acad"public Frame3d CurveFrameAt(Entity ent, double t)"
acad"public Frame3d CurveFrameAtLength(Entity ent, double l)"

backend_map_division(b::ACAD, f::Function, s::Shape1D, n::Int) =
    let conn = connection(b)
        r = ref(s).value
        (t1, t2) = ACADCurveDomain(conn, r)
        map_division(t1, t2, n) do t
            f(ACADCurveFrameAt(conn, r, t))
        end
    end


acad"public Vector3d RegionNormal(Entity ent)"
acad"public Point3d RegionCentroid(Entity ent)"
acad"public double[] SurfaceDomain(Entity ent)"
acad"public Frame3d SurfaceFrameAt(Entity ent, double u, double v)"

backend_surface_domain(b::ACAD, s::Shape2D) =
    tuple(ACADSurfaceDomain(connection(b), ref(s).value)...)

backend_map_division(b::ACAD, f::Function, s::Shape2D, nu::Int, nv::Int) =
    let conn = connection(b)
        r = ref(s).value
        (u1, u2, v1, v2) = ACADSurfaceDomain(conn, r)
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(ACADSurfaceFrameAt(conn, r, u, v))
            end
        end
    end

# The previous method cannot be applied to meshes in AutoCAD, which are created by surface_grid

backend_map_division(b::ACAD, f::Function, s::SurfaceGrid, nu::Int, nv::Int) =
    let (u1, u2, v1, v2) = ACADSurfaceDomain(conn, r)
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(ACADSurfaceFrameAt(conn, r, u, v))
            end
        end
    end

realize(b::ACAD, s::Text) =
  ACADText(connection(b), s.str, s.c, vx(1, s.c.cs), vy(1, s.c.cs), s.h)

realize(b::ACAD, s::Sphere) =
  ACADSphere(connection(b), s.center, s.radius)
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
  ACADIrregularPyramid(connection(b), s.cbs, s.ct)
realize(b::ACAD, s::RegularPrism) =
  let cbs = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
    ACADIrregularPyramidFrustum(connection(b),
                                   cbs,
                                   map(p -> add_z(p, s.h), cbs))
  end
realize(b::ACAD, s::IrregularPyramidFustrum) =
    ACADIrregularPyramidFrustum(connection(b), s.cbs, s.cts)

realize(b::ACAD, s::IrregularPrism) =
  ACADIrregularPyramidFrustum(connection(b),
                              s.cbs,
                              map(p -> (p + s.v), s.cbs))
realize(b::ACAD, s::RightCuboid) =
  ACADCenteredBox(connection(b), s.cb, s.width, s.height, s.h)
realize(b::ACAD, s::Box) =
  ACADBox(connection(b), s.c, s.dx, s.dy, s.dz)
realize(b::ACAD, s::Cone) =
  ACADCone(connection(b), add_z(s.cb, s.h), s.r, s.cb)
realize(b::ACAD, s::ConeFrustum) =
  ACADConeFrustum(connection(b), s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt)
realize(b::ACAD, s::Cylinder) =
  ACADCylinder(connection(b), s.cb, s.r, s.cb + vz(s.h, s.cb.cs))
#realize(b::ACAD, s::Circle) = ACADCircle(connection(b),

realize(b::Backend, s::Extrusion) =
  backend_extrusion(b, s.profile, s.v)

backend_extrusion(b::Backend, p::Point, v::Vec) =
  realize_and_delete_shapes(line([p.position, p.position + v], backend=b), [p])

backend_extrusion(b::Backend, s::Shape, v::Vec) =
    and_mark_deleted(
        map_ref(s) do r
            ACADExtrude(connection(b), r, v)
        end,
        s)

realize(b::Backend, s::Sweep) =
  backend_sweep(b, s.path, s.profile, s.rotation, s.scale)

backend_sweep(b::ACAD, path::Shape, profile::Shape, rotation::Real, scale::Real) =
  map_ref(profile) do profile_r
    map_ref(path) do path_r
      ACADSweep(connection(b), path_r, profile_r, rotation, scale)
    end
  end

realize(b::ACAD, s::Revolve) =
  and_delete_shape(
    map_ref(s.profile) do r
      ACADRevolve(connection(b), r, s.p, s.n, s.start_angle, s.amplitude)
    end,
    s.profile)

backend_loft_points(b::Backend, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  let f = (ruled ? (closed ? polygon : line) : (closed ? closed_spline : spline))
    and_delete_shapes(ref(f(map(point_position, profiles), backend=b)),
                      vcat(profiles, rails))
  end

backend_loft_curves(b::ACAD, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  and_delete_shapes(ACADLoft(connection(b),
                             collect_ref(profiles),
                             collect_ref(rails),
                             ruled, closed),
                    vcat(profiles, rails))

backend_loft_surfaces(b::ACAD, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
    backend_loft_curves(b, profiles, rails, ruled, closed)

backend_loft_curve_point(b::ACAD, profile::Shape, point::Shape) =
    and_delete_shapes(ACADLoft(connection(b),
                               vcat(collect_ref(profile), collect_ref(point)),
                               [],
                               true, false),
                      [profile, point])

backend_loft_surface_point(b::ACAD, profile::Shape, point::Shape) =
    backend_loft_curve_point(b, profile, point)

unite_ref(b::ACAD, r0::ACADNativeRef, r1::ACADNativeRef) =
    ensure_ref(b, ACADUnite(connection(b), r0.value, r1.value))

intersect_ref(b::ACAD, r0::ACADNativeRef, r1::ACADNativeRef) =
    ensure_ref(b, ACADIntersect(connection(b), r0.value, r1.value))

subtract_ref(b::ACAD, r0::ACADNativeRef, r1::ACADNativeRef) =
    ensure_ref(b, ACADSubtract(connection(b), r0.value, r1.value))

slice_ref(b::ACAD, r::ACADNativeRef, p::Loc, v::Vec) =
    (ACADSlice(connection(b), r.value, p, v); r)

slice_ref(b::ACAD, r::ACADUnionRef, p::Loc, v::Vec) =
    map(r->slice_ref(b, r, p, v), r.values)

unite_refs(b::ACAD, refs::Vector{<:ACADRef}) =
    ACADUnionRef(tuple(refs...))


realize(b::ACAD, s::IntersectionShape) =
  foldl((r0,r1)->intersect_ref(b,r0,r1), ACADUniversalRef(), map(ref, s.shapes))

realize(b::ACAD, s::Slice) =
  slice_ref(b, ref(s.shape), s.p, s.n)

realize(b::ACAD, s::Move) =
  let r = map_ref(s.shape) do r
            ACADMove(connection(b), r, s.v)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Transform) =
  let r = map_ref(s.shape) do r
            ACADTransform(connection(b), r, s.xform)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Scale) =
  let r = map_ref(s.shape) do r
            ACADScale(connection(b), r, s.p, s.s)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Rotate) =
  let r = map_ref(s.shape) do r
            ACADRotate(connection(b), r, s.p, s.v, s.angle)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Mirror) =
  and_delete_shape(map_ref(s.shape) do r
                    ACADMirror(connection(b), r, s.p, s.n, false)
                   end,
                   s.shape)

realize(b::ACAD, s::UnionMirror) =
  let r0 = ref(s.shape),
      r1 = map_ref(s.shape) do r
            ACADMirror(connection(b), r, s.p, s.n, true)
          end
    UnionRef((r0,r1))
  end

realize(b::ACAD, s::SurfaceGrid) =
    ACADSurfaceFromGrid(
        connection(b),
        size(s.points,1),
        size(s.points,2),
        reshape(s.points,:),
        s.closed_u,
        s.closed_v,
        2)

realize(b::ACAD, s::Thicken) =
  and_delete_shape(
    map_ref(s.shape) do r
      ACADThicken(connection(b), r, s.thickness)
    end,
    s.shape)

# backend_frame_at
backend_frame_at(b::ACAD, s::Circle, t::Real) = add_pol(s.center, s.radius, t)

backend_frame_at(b::ACAD, c::Shape1D, t::Real) = ACADCurveFrameAt(connection(b), ref(c).value, t)

#backend_frame_at(b::ACAD, s::Surface, u::Real, v::Real) =
    #What should we do with v?
#    backend_frame_at(b, s.frontier[1], u)

#backend_frame_at(b::ACAD, s::SurfacePolygon, u::Real, v::Real) =

backend_frame_at(b::ACAD, s::Shape2D, u::Real, v::Real) = ACADSurfaceFrameAt(connection(b), ref(s).value, u, v)

# BIM
backend_get_family(b::ACAD, f::TableFamily) =
    ACADCreateRectangularTableFamily(connection(b), f.length, f.width, f.height, f.top_thickness, f.leg_thickness)
backend_get_family(b::ACAD, f::ChairFamily) =
    ACADCreateChairFamily(connection(b), f.length, f.width, f.height, f.seat_height, f.thickness)
backend_get_family(b::ACAD, f::TableChairFamily) =
    ACADCreateRectangularTableAndChairsFamily(connection(b),
        ref(f.table_family), ref(f.chair_family),
        f.table_family.length, f.table_family.width,
        f.chairs_top, f.chairs_bottom, f.chairs_right, f.chairs_left,
        f.spacing)

backend_rectangular_table(b::ACAD, c, angle, family) =
    ACADTable(connection(b), c, angle, ref(family))

backend_chair(b::ACAD, c, angle, family) =
    ACADChair(connection(b), c, angle, ref(family))

backend_rectangular_table_and_chairs(b::ACAD, c, angle, family) =
    ACADTableAndChairs(connection(b), c, angle, ref(family))

backend_slab(b::ACAD, profile, thickness) =
    map_ref(b,
            r->ACADExtrude(connection(b), r, vz(thickness)),
            ensure_ref(b, backend_fill(b, profile)))

#Beams are aligned along the top axis.
realize(b::ACAD, s::Beam) =
    let o = loc_from_o_phi(s.cb, s.angle)
        ACADCenteredBox(connection(b), add_y(o, -s.family.height/2), s.family.width, s.family.height, s.h)
    end
#    ACADCenteredBox(connection(b), s.cb, vx(1, s.cb.cs), vy(1, s.cb.cs), s.family.width, s.family.height, s.h)

#Columns are aligned along the center axis.
realize(b::ACAD, s::Column) =
    let o = loc_from_o_phi(s.cb, s.angle)
        ACADCenteredBox(connection(b), o, s.family.width, s.family.height, s.h)
    end

backend_wall(b::ACAD, path, height, thickness) =
    let conn = connection(b)
        ACADThicken(conn,
                    ACADExtrude(conn,
                                backend_stroke(b, path),
                                vz(height)),
                    thickness)
    end
#=
backend_wall(b::ACAD, path, height, thickness) =
    let conn = connection(b)
        ACADSweep(conn,
                  backend_stroke(b, path),
                  ACADPolyLine(conn, [xy(thickness/-2,0), xy(thickness/+2,0)]),
                  0,
                  1)
    end
=#
###
#=
sweep_fractions(conn, b, verts, thickness) =
    begin
        ACADSweep(conn,
                  ACADPolyLine(conn, verts[1:2]),
                  ACADPolyLine(conn, [xy(thickness/-2,0), xy(thickness/+2,0)]),
                  0,
                  1)
        if verts.length >= 2
            sweep_fractions(conn, b, verts[1:end], thickness)
        end
    end

#
sweep_fractions(b, verts, thickness) =
    let p = verts[1]
        q = verts[2]
        o = loc_from_o_phi(p, pol_phi(q-p))
        if length(verts) == 2
            stroke(rectangular_path(add_y(o, thickness/-2), distance(p, q), thickness), b)
        else
            sweep_fractions(b, verts[2:end], thickness)
        end
    end

backend_wall(b::ACAD, path, height, thickness) =
    backend_wall_path(b, path, height, thickness)

backend_wall_path(b::ACAD, path::RectangularPath, height, thickness) =
    stroke(path, b)
backend_wall_path(b::ACAD, path::OpenPolygonalPath, height, thickness) =
    sweep_fractions(b, path.vertices, thickness)
=#

############################################

backend_bounding_box(b::ACAD, shapes::Shapes) =
  ACADBoundingBox(connection(b), collect_ref(shapes))


backend_name(b::ACAD) = "AutoCAD"

Base.view(camera::Loc, target::Loc, lens::Real, b::ACAD) =
  ACADView(connection(b), camera, target, lens)

get_view(b::ACAD) =
  let c = connection(b)
    ACADViewCamera(c), ACADViewTarget(c), ACADViewLens(c)
  end

zoom_extents(b::ACAD) = ACADZoomExtents(connection(b))

view_top(b::ACAD) = ACADViewTop(connection(b))

backend_delete_shapes(b::ACAD, shapes::Shapes) =
  ACADDeleteMany(connection(b), collect_ref(shapes))

delete_all_shapes(b::ACAD) = ACADDeleteAll(connection(b))
set_length_unit(unit::String, b::ACAD) = ACADSetLengthUnit(connection(b), unit)

# Dimensions

const ACADDimensionStyles = Dict(:architectural => "_ARCHTICK", :mechanical => "")

dimension(p0::Loc, p1::Loc, p::Loc, scale::Real, style::Symbol, b::ACAD=current_backend()) =
    ACADCreateAlignedDimension(connection(b), p0, p1, p,
        scale,
        ACADDimensionStyles[style])

dimension(p0::Loc, p1::Loc, sep::Real, scale::Real, style::Symbol, b::ACAD=current_backend()) =
    let v = p1 - p0
        angle = pol_phi(v)
        dimension(p0, p1, add_pol(p0, sep, angle + pi/2), scale, style, b)
    end

# Layers
ACADLayer = Int

current_layer(b::ACAD=current_backend())::ACADLayer =
  ACADCurrentLayer(connection(b))

current_layer(layer::ACADLayer, b::ACAD=current_backend()) =
  ACADSetCurrentLayer(connection(b), layer)

create_layer(name::String, b::ACAD=current_backend()) =
  ACADCreateLayer(connection(b), name)

create_layer(name::String, color::RGB, b::ACAD=current_backend()) =
  let layer = ACADCreateLayer(connection(b), name)
    ACADSetLayerColor(connection(b), layer, color.r, color.g, color.b)
    layer
  end

# Blocks

realize(b::ACAD, s::Block) =
    ACADCreateBlockFromShapes(connection(b), s.name, collect_ref(s.shapes))

backend_create_block(name::String, shapes::Shapes, b::ACAD=current_backend()) =
    ACADCreateBlockFromShapes(connection(b), name, collect_ref(shapes))

realize(b::ACAD, s::BlockInstance) =
    ACADCreateBlockInstance(
        connection(b),
        collect_ref(s.block)[1],
        center_scaled_cs(s.loc, s.scale, s.scale, s.scale))

#=

# Manual process
@time for i in 1:1000 for r in 1:10 circle(x(i*10), r) end end

# Create block...
Khepri.create_block("Foo", [circle(radius=r) for r in 1:10])

# ...and instantiate it
@time for i in 1:1000 Khepri.instantiate_block("Foo", x(i*10), 0) end

=#

# Lights
acad"public Entity SpotLight(Point3d position, double hotspot, double falloff, Point3d target)"
acad"public Entity IESLight(String webFile, Point3d position, Point3d target, Vector3d rotation)"

backend_spotlight(b::ACAD, loc::Loc, dir::Vec, hotspot::Real, falloff::Real) =
    ACADSpotLight(connection(b), loc, hotspot, falloff, loc + dir)

backend_ieslight(b::ACAD, file::String, loc::Loc, dir::Vec, alpha::Real, beta::Real, gamma::Real) =
    ACADIESLight(connection(b), file, loc, loc + dir, vxyz(alpha, beta, gamma))



prompt_position(prompt::String, b::ACAD) =
  let ans = ACADGetPoint(connection(b), prompt)
    length(ans) > 0 && ans[1]
  end

shape_from_ref(r, b::ACAD=current_backend()) =
    let c = connection(b)
        code = ACADShapeCode(c, r)
        ref = LazyRef(b, ACADNativeRef(r))
        if code == 1 # Point
            point(ACADPointPosition(c, r),
                  backend=b, ref=ref)
        elseif code == 2
            circle(maybe_loc_from_o_vz(ACADCircleCenter(c, r), ACADCircleNormal(c, r)),
                   ACADCircleRadius(c, r),
                   backend=b, ref=ref)
        elseif 3 <= code <= 6
            line(ACADLineVertices(c, r),
                 backend=b, ref=ref)
        elseif code == 7
            spline([xy(0,0)], false, false, #HACK obtain interpolation points
                   backend=b, ref=ref)
        elseif code == 9
            let start_angle = mod(ACADArcStartAngle(c, r), 2pi)
                end_angle = mod(ACADArcEndAngle(c, r), 2pi)
                if end_angle > start_angle
                    arc(maybe_loc_from_o_vz(ACADArcCenter(c, r), ACADArcNormal(c, r)),
                        ACADArcRadius(c, r), start_angle, end_angle - start_angle,
                        backend=b, ref=ref)
                else
                    arc(maybe_loc_from_o_vz(ACADArcCenter(c, r), ACADArcNormal(c, r)),
                        ACADArcRadius(c, r), end_angle, start_angle - end_angle,
                        backend=b, ref=ref)
                end
            end
        elseif code == 10
            let str = ACADTextString(c, r)
                height = ACADTextHeight(c, r)
                loc = ACADTextPosition(c, r)
                text(str, loc, height, backend=b, ref=ref)
            end
        elseif code == 11
            let str = ACADMTextString(c, r)
                height = ACADMTextHeight(c, r)
                loc = ACADMTextPosition(c, r)
                text(str, loc, height, backend=b, ref=ref)
            end
        elseif code == 12 || code == 13
            surface(Shapes1D[], backend=b, ref=ref)
        elseif code == 50
            block_instance(block("To be finished!"))
        elseif code == 70
            block_instance(block("A viewport to be finished!"))
        elseif 103 <= code <= 106
            polygon(ACADLineVertices(c, r),
                    backend=b, ref=ref)
        else
            block_instance(block("Unknown shape. To be finished!"))
            #error("Unknown shape with code $(code)")
        end
    end

all_shapes(b::ACAD) =
  [shape_from_ref(r, b) for r in ACADGetAllShapes(connection(b))]

all_shapes_in_layer(layer, b::ACAD) =
  [shape_from_ref(r, b) for r in ACADGetAllShapesInLayer(connection(b), layer)]

disable_update(b::ACAD) =
    ACADDisableUpdate(connection(b))

enable_update(b::ACAD) =
    ACADEnableUpdate(connection(b))
# Render

acad"public void Render(int width, int height, string path, int levels, double exposure)"
#render exposure: [-3, +3] -> [-6, 21]
convert_render_exposure(b::ACAD, v::Real) = -4.05*v + 8.8
#render quality: [-1, +1] -> [+1, +50]
convert_render_quality(b::ACAD, v::Real) = round(Int, 25.5 + 24.5*v)

render_view(name::String, b::ACAD=current_backend()) =
    ACADRender(connection(b),
               render_width(), render_height(),
               prepare_for_saving_file(render_pathname(name)),
               convert_render_quality(b, render_quality()),
               convert_render_exposure(b, render_exposure()))


save_as(pathname::String, format::String, b::ACAD) =
    ACADSaveAs(connection(b), pathname, format)
