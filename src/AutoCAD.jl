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

env(name) = Sys.iswindows() ? ENV[name] : ""

autocad_general_plugins = joinpath(dirname(env("CommonProgramFiles")), "Autodesk", "ApplicationPlugins")
autocad_allusers_plugins = joinpath(env("ALLUSERSPROFILE"), "Autodesk", "ApplicationPlugins")
autocad_user_plugins = joinpath(env("APPDATA"), "Autodesk", "ApplicationPlugins")
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
      @info("Checking AutoCAD plugin...")
      for i in 1:10
        try
          update_plugin()
          @info("done.")
          checked_plugin = true
          return
        catch exc
          if isa(exc, Base.IOError) && i < 10
            @error("The AutoCAD plugin is outdated! Please, close AutoCAD.")
            sleep(5)
          else
            throw(exc)
          end
        end
      end
    end
  end

# We need some additional Encoders
encode_Entity = encode_int
decode_Entity = decode_int_or_error
encode_ObjectId = encode_int
decode_ObjectId = decode_int_or_error
encode_ObjectId_array = encode_int_array
decode_ObjectId_array = decode_int_array

acad_api = @remote_functions :CS """
public void SetLengthUnit(String unit)
public void SetView(Point3d position, Point3d target, double lens, bool perspective, string style)
public void View(Point3d position, Point3d target, double lens)
public void ViewTop()
public Point3d ViewCamera()
public Point3d ViewTarget()
public double ViewLens()
public void SetSkyFromDateLocation(int year, int month, int day, int hour, int minute, double latitude, double longitude, double meridian)
public byte Sync()
public byte Disconnect()
public ObjectId Copy(ObjectId id)
public Entity Point(Point3d p)
public Point3d PointPosition(Entity ent)
public Entity PolyLine(Point3d[] pts)
public Point3d[] LineVertices(ObjectId id)
public Entity Spline(Point3d[] pts)
public Entity InterpSpline(Point3d[] pts, Vector3d tan0, Vector3d tan1)
public Entity ClosedPolyLine(Point3d[] pts)
public Entity ClosedSpline(Point3d[] pts)
public Entity InterpClosedSpline(Point3d[] pts)
public Point3d[] SplineInterpPoints(Entity ent)
public Vector3d[] SplineTangents(Entity ent)
public Entity Circle(Point3d c, Vector3d n, double r)
public Point3d CircleCenter(Entity ent)
public Vector3d CircleNormal(Entity ent)
public double CircleRadius(Entity ent)
public Entity Ellipse(Point3d c, Vector3d n, Vector3d majorAxis, double radiusRatio)
public Entity Arc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)
public Point3d ArcCenter(Entity ent)
public Vector3d ArcNormal(Entity ent)
public double ArcRadius(Entity ent)
public double ArcStartAngle(Entity ent)
public double ArcEndAngle(Entity ent)
public ObjectId JoinCurves(ObjectId[] ids)
public Entity Text(string str, Point3d corner, Vector3d vx, Vector3d vy, double height)
public Entity SurfaceFromCurve(Entity curve)
public Entity SurfaceCircle(Point3d c, Vector3d n, double r)
public Entity SurfaceEllipse(Point3d c, Vector3d n, Vector3d majorAxis, double radiusRatio)
public Entity SurfaceArc(Point3d c, Vector3d n, double radius, double startAngle, double endAngle)
public Entity SurfaceClosedPolyLine(Point3d[] pts)
public ObjectId[] SurfaceFromCurves(ObjectId[] ids)
public ObjectId[] CurvesFromSurface(ObjectId id)
public Entity Sphere(Point3d c, double r)
public Entity Torus(Point3d c, Vector3d vz, double majorRadius, double minorRadius)
public Entity ConeFrustum(Point3d bottom, double base_radius, Point3d top, double top_radius)
public Entity Cylinder(Point3d bottom, double radius, Point3d top)
public Entity Cone(Point3d bottom, double radius, Point3d top)
public Entity Box(Frame3d frame, double dx, double dy, double dz)
public Entity CenteredBox(Frame3d frame, double dx, double dy, double dz)
public ObjectId IrregularPyramidMesh(Point3d[] pts, Point3d apex)
public ObjectId IrregularPyramid(Point3d[] pts, Point3d apex)
public ObjectId IrregularPyramidFrustum(Point3d[] bpts, Point3d[] tpts)
public ObjectId Thicken(ObjectId obj, double thickness)
public ObjectId NurbSurfaceFrom(ObjectId id)
public ObjectId Extrude(ObjectId profileId, Vector3d dir)
public ObjectId Sweep(ObjectId pathId, ObjectId profileId, double rotation, double scale)
public ObjectId Loft(ObjectId[] profilesIds, ObjectId[] guidesIds, bool ruled, bool closed)
public ObjectId Unite(ObjectId objId0, ObjectId objId1)
public ObjectId Intersect(ObjectId objId0, ObjectId objId1)
public ObjectId Subtract(ObjectId objId0, ObjectId objId1)
public void Slice(ObjectId id, Point3d p, Vector3d n)
public ObjectId Revolve(ObjectId profileId, Point3d p, Vector3d n, double startAngle, double amplitude)
public void Transform(ObjectId id, Frame3d frame)
public void Move(ObjectId id, Vector3d v)
public void Scale(ObjectId id, Point3d p, double s)
public void Rotate(ObjectId id, Point3d p, Vector3d n, double a)
public ObjectId Mirror(ObjectId id, Point3d p, Vector3d n, bool copy)
public Point3d[] BoundingBox(ObjectId[] ids)
public void ZoomExtents()
public ObjectId CreateLayer(string name, bool active, byte r, byte g, byte b)
public void SetLayerColor(ObjectId id, byte r, byte g, byte b)
public void SetShapeColor(ObjectId id, byte r, byte g, byte b)
public ObjectId CurrentLayer()
public void SetCurrentLayer(ObjectId id)
public ObjectId ShapeLayer(ObjectId objId)
public void SetShapeLayer(ObjectId objId, ObjectId layerId)
public void SetSystemVariableInt(string name, int value)
public int Command(string cmd)
public void DisableUpdate()
public void EnableUpdate()
public bool IsPoint(Entity e)
public bool IsCircle(Entity e)
public bool IsPolyLine(Entity e)
public bool IsSpline(Entity e)
public bool IsInterpSpline(Entity e)
public bool IsClosedPolyLine(Entity e)
public bool IsClosedSpline(Entity e)
public bool IsInterpClosedSpline(Entity e)
public bool IsEllipse(Entity e)
public bool IsArc(Entity e)
public bool IsText(Entity e)
public byte ShapeCode(ObjectId id)
public BIMLevel FindOrCreateLevelAtElevation(double elevation)
public BIMLevel UpperLevel(BIMLevel currentLevel, double addedElevation)
public double GetLevelElevation(BIMLevel level)
public FloorFamily FloorFamilyInstance(double totalThickness, double coatingThickness)
public Entity LightweightPolyLine(Point2d[] pts, double[] angles, double elevation)
public Entity SurfaceLightweightPolyLine(Point2d[] pts, double[] angles, double elevation)
public ObjectId CreatePathFloor(Point2d[] pts, double[] angles, BIMLevel level, FloorFamily family)
public ObjectId CreateBlockFromShapes(String baseName, ObjectId[] ids)
public ObjectId CreateBlockInstance(ObjectId id, Frame3d frame)
public ObjectId CreateInstanceFromBlockNamed(String name, Frame3d frame)
public ObjectId CreateInstanceFromBlockNamedAtRotated(String name, Point3d c, double angle)
public ObjectId CreateRectangularTableFamily(double length, double width, double height, double top_thickness, double leg_thickness)
public ObjectId Table(Point3d c, double angle, ObjectId family)
public ObjectId CreateChairFamily(double length, double width, double height, double seat_height, double thickness)
public ObjectId Chair(Point3d c, double angle, ObjectId family)
public ObjectId CreateRectangularTableAndChairsFamily(ObjectId tableFamily, ObjectId chairFamily, double tableLength, double tableWidth, int chairsOnTop, int chairsOnBottom, int chairsOnRight, int chairsOnLeft, double spacing)
public ObjectId TableAndChairs(Point3d c, double angle, ObjectId family)
public ObjectId CreateAlignedDimension(Point3d p0, Point3d p1, Point3d p, double scale, String mark)
public String TextString(Entity ent)
public Point3d TextPosition(Entity ent)
public double TextHeight(Entity ent)
public String MTextString(Entity ent)
public Point3d MTextPosition(Entity ent)
public double MTextHeight(Entity ent)
public void SaveAs(String pathname, String format)
public double[] CurveDomain(Entity ent)
public double CurveLength(Entity ent)
public Frame3d CurveFrameAt(Entity ent, double t)
public Frame3d CurveFrameAtLength(Entity ent, double l)
public Point3d[] CurvePointsAt(Entity ent, double[] ts)
public Vector3d[] CurveTangentsAt(Entity ent, double[] ts)
public Vector3d[] CurveNormalsAt(Entity ent, double[] ts)
public Vector3d RegionNormal(Entity ent)
public Point3d RegionCentroid(Entity ent)
public double[] SurfaceDomain(Entity ent)
public Frame3d SurfaceFrameAt(Entity ent, double u, double v)
public Entity MeshFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN)
public int[] PolygonMeshData(Entity e)
public Point3dCollection MeshVertices(ObjectId id)
public Entity SurfaceFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level)
public Entity SolidFromGrid(int m, int n, Point3d[] pts, bool closedM, bool closedN, int level, double thickness)
public void DeleteAll()
public void DeleteAllInLayer(ObjectId layerId)
public void Delete(ObjectId id)
public void DeleteMany(ObjectId[] ids)
public Entity SpotLight(Point3d position, double hotspot, double falloff, Point3d target)
public Entity IESLight(String webFile, Point3d position, Point3d target, Vector3d rotation)
public Point3d[] GetPosition(string prompt)
public ObjectId[] GetPoint(string prompt)
public ObjectId[] GetPoints(string prompt)
public ObjectId[] GetCurve(string prompt)
public ObjectId[] GetCurves(string prompt)
public ObjectId[] GetSurface(string prompt)
public ObjectId[] GetSurfaces(string prompt)
public ObjectId[] GetSolid(string prompt)
public ObjectId[] GetSolids(string prompt)
public ObjectId[] GetShape(string prompt)
public ObjectId[] GetShapes(string prompt)
public ObjectId[] GetPreSelectedShapes()
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
public void SelectShapes(ObjectId[] ids)
public void Render(int width, int height, string path, int levels, double exposure)
"""

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
const ACAD = SocketBackend{ACADKey, ACADId}

void_ref(b::ACAD) = ACADNativeRef(-1)

create_ACAD_connection() =
    begin
        check_plugin()
        create_backend_connection("AutoCAD", 11000)
    end

const autocad = ACAD(LazyParameter(TCPSocket, create_ACAD_connection), acad_api)

backend_name(b::ACAD) = "AutoCAD"

#=

Default families

=#


abstract type ACADFamily <: Family end

struct ACADLayerFamily <: ACADFamily
  name::String
  color::RGB
  ref::Parameter{Any}
end

acad_layer_family(name, color::RGB=rgb(1,1,1)) =
  ACADLayerFamily(name, color, Parameter{Any}(nothing))

backend_get_family_ref(b::ACAD, f::Family, af::ACADLayerFamily) =
  backend_create_layer(b, af.name, true, af.color)

set_backend_family(default_wall_family(), autocad, acad_layer_family("Wall"))
set_backend_family(default_slab_family(), autocad, acad_layer_family("Slab"))
set_backend_family(default_roof_family(), autocad, acad_layer_family("Roof"))
set_backend_family(default_beam_family(), autocad, acad_layer_family("Beam"))
set_backend_family(default_column_family(), autocad, acad_layer_family("Column"))
set_backend_family(default_door_family(), autocad, acad_layer_family("Door"))
set_backend_family(default_panel_family(), autocad, acad_layer_family("Panel"))

set_backend_family(default_table_family(), autocad, acad_layer_family("Table"))
set_backend_family(default_chair_family(), autocad, acad_layer_family("Chair"))
set_backend_family(default_table_chair_family(), autocad, acad_layer_family("TableChairs"))

set_backend_family(default_curtain_wall_family(), autocad, acad_layer_family("CurtainWall"))
set_backend_family(default_curtain_wall_family().panel, autocad, acad_layer_family("CurtainWall-Panel"))
set_backend_family(default_curtain_wall_family().boundary_frame, autocad, acad_layer_family("CurtainWall-Boundary"))
set_backend_family(default_curtain_wall_family().transom_frame, autocad, acad_layer_family("CurtainWall-Transom"))
set_backend_family(default_curtain_wall_family().mullion_frame, autocad, acad_layer_family("CurtainWall-Mullion"))
#current_backend(autocad)

set_backend_family(default_truss_node_family(), autocad, acad_layer_family("TrussNode"))
set_backend_family(default_truss_bar_family(), autocad, acad_layer_family("TrussBar"))

set_backend_family(fixed_truss_node_family, autocad, acad_layer_family("FixedTrussNode"))
set_backend_family(free_truss_node_family, autocad, acad_layer_family("FreeTrussNode"))

#use_family_in_layer(b::ACAD) = true

backend_stroke_color(b::ACAD, path::Path, color::RGB) =
    let r = backend_stroke(b, path)
        @remote(b, SetShapeColor(r, color.r, color.g, color.b))
        r
    end

backend_stroke(b::ACAD, path::CircularPath) =
    @remote(b, Circle(path.center, vz(1, path.center.cs), path.radius))
backend_stroke(b::ACAD, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        @remote(b, ClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)]))
    end
backend_stroke(b::ACAD, path::ArcPath) =
    backend_stroke_arc(b, path.center, path.radius, path.start_angle, path.amplitude)

backend_stroke(b::ACAD, path::OpenPolygonalPath) =
  	@remote(b, PolyLine(path.vertices))
backend_stroke(b::ACAD, path::ClosedPolygonalPath) =
    backend_polygon(b, path.vertices)
backend_polygon(b::ACAD, vs::Locs) =
  @remote(b, ClosedPolyLine(vs))
backend_fill(b::ACAD, path::ClosedPolygonalPath) =
    @remote(b, SurfaceClosedPolyLine(path.vertices))
backend_fill(b::ACAD, path::RectangularPath) =
    let c = path.corner,
        dx = path.dx,
        dy = path.dy
        @remote(b, SurfaceClosedPolyLine([c, add_x(c, dx), add_xy(c, dx, dy), add_y(c, dy)]))
    end
backend_stroke(b::ACAD, path::OpenSplinePath) =
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
backend_stroke(b::ACAD, path::ClosedSplinePath) =
    @remote(b, InterpClosedSpline(path.vertices))
backend_fill(b::ACAD, path::ClosedSplinePath) =
    backend_fill_curves(b, @remote(b, InterpClosedSpline(path.vertices)))

backend_fill_curves(b::ACAD, refs::ACADIds) = @remote(b, SurfaceFromCurves(refs))
backend_fill_curves(b::ACAD, ref::ACADId) = @remote(b, SurfaceFromCurves([ref]))

backend_stroke_arc(b::ACAD, center::Loc, radius::Real, start_angle::Real, amplitude::Real) =
  let end_angle = start_angle + amplitude
    if end_angle > start_angle
      @remote(b, Arc(center, vz(1, center.cs), radius, start_angle, end_angle))
    else
      @remote(b, Arc(center, vz(1, center.cs), radius, end_angle, start_angle))
    end

  end
backend_stroke_unite(b::ACAD, refs) = @remote(b, JoinCurves(refs))



realize(b::ACAD, s::EmptyShape) =
  ACADEmptyRef()
realize(b::ACAD, s::UniversalShape) =
  ACADUniversalRef()
realize(b::ACAD, s::Point) =
  @remote(b, Point(s.position))
realize(b::ACAD, s::Line) =
  @remote(b, PolyLine(s.vertices))
realize(b::ACAD, s::Spline) = # This should be merged with opensplinepath
  if (s.v0 == false) && (s.v1 == false)
    #@remote(b, Spline(s.points))
    @remote(b, InterpSpline(
                     s.points,
                     s.points[2]-s.points[1],
                     s.points[end]-s.points[end-1]))
  elseif (s.v0 != false) && (s.v1 != false)
    @remote(b, InterpSpline(s.points, s.v0, s.v1))
  else
    @remote(b, InterpSpline(
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1))
  end
realize(b::ACAD, s::ClosedSpline) =
  @remote(b, InterpClosedSpline(s.points))
realize(b::ACAD, s::Circle) =
  @remote(b, Circle(s.center, vz(1, s.center.cs), s.radius))
realize(b::ACAD, s::Arc) =
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

realize(b::ACAD, s::Ellipse) =
  if s.radius_x > s.radius_y
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x))
  else
    @remote(b, Ellipse(s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y))
  end
realize(b::ACAD, s::EllipticArc) =
  error("Finish this")

realize(b::ACAD, s::Polygon) =
  @remote(b, ClosedPolyLine(s.vertices))
realize(b::ACAD, s::RegularPolygon) =
  @remote(b, ClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))
realize(b::ACAD, s::Rectangle) =
  @remote(b, ClosedPolyLine(
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)]))
realize(b::ACAD, s::SurfaceCircle) =
  @remote(b, SurfaceCircle(s.center, vz(1, s.center.cs), s.radius))
realize(b::ACAD, s::SurfaceArc) =
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

realize(b::ACAD, s::SurfaceEllipse) =
  if s.radius_x > s.radius_y
    @remote(b, SurfaceEllipse(s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x))
  else
    @remote(b, SurfaceEllipse(s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y))
  end


backend_surface_polygon(b::ACAD, vs::Locs) =
  @remote(b, SurfaceClosedPolyLine(vs))
realize(b::ACAD, s::Surface) =
  let #ids = map(r->@remote(b, NurbSurfaceFrom(r)), @remote(b, SurfaceFromCurves(collect_ref(s.frontier))))
      ids = @remote(b, SurfaceFromCurves(collect_ref(s.frontier)))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::ACAD, s::Shape2D) =
    map(c -> shape_from_ref(c, b), @remote(b, CurvesFromSurface(ref(s).value)))

# Iterating over curves and surfaces


old_backend_map_division(b::ACAD, f::Function, s::Shape1D, n::Int) =
  let r = ref(s).value,
      (t1, t2) = @remote(b, CurveDomain(r))
    map_division(t1, t2, n) do t
      f(@remote(b, CurveFrameAt(r, t)))
    end
  end

# For low level access:

backend_map_division(b::ACAD, f::Function, s::Shape1D, n::Int) =
  let r = ref(s).value,
      (t1, t2) = @remote(b, CurveDomain(r)),
      ti = division(t1, t2, n),
      ps = @remote(b, CurvePointsAt(r, ti)),
      ts = @remote(b, CurveTangentsAt(r, ti)),
      #ns = @remote(b, CurveNormalsAt(r, ti)),
      frames = rotation_minimizing_frames(@remote(b, CurveFrameAt(r, t1)), ps, ts)
    map(f, frames)
  end

#=
rotation_minimizing_frames(u0, xs, ts) =
  let ri = in_world(vy(1, u0.cs)),
      new_frames = [loc_from_o_vx_vy(xs[1], ri, cross(ts[1], ri))]
    for i in 1:length(xs)-1
      let xi = xs[i],
          xii = xs[i+1],
          ti = ts[i],
          tii = ts[i+1],
          v1 = xii - xi,
          c1 = dot(v1, v1),
          ril = ri - v1*(2/c1*dot(v1,ri)),
          til = ti - v1*(2/c1*dot(v1,ti)),
          v2 = tii - til,
          c2 = dot(v2, v2),
          rii = ril - v2*(2/c2*dot(v2, ril)),
          sii = cross(tii, rii),
          uii = loc_from_o_vx_vy(xii, rii, sii)
        push!(new_frames, uii)
        ri = rii
      end
    end
    new_frames
  end
=#

#


backend_surface_domain(b::ACAD, s::Shape2D) =
    tuple(@remote(b, SurfaceDomain(ref(s).value))...)

backend_map_division(b::ACAD, f::Function, s::Shape2D, nu::Int, nv::Int) =
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


backend_map_division(b::ACAD, f::Function, s::SurfaceGrid, nu::Int, nv::Int) =
let conn = connection(b)
    r = ref(s).value
    (u1, u2, v1, v2) = @remote(b, SurfaceDomain(r))
    map_division(u1, u2, nu) do u
        map_division(v1, v2, nv) do v
            f(@remote(b, SurfaceFrameAt(r, u, v)))
        end
    end
end

realize(b::ACAD, s::Text) =
  @remote(b, Text(
    s.str, s.corner, vx(1, s.corner.cs), vy(1, s.corner.cs), s.height))

realize(b::ACAD, s::Sphere) =
  @remote(b, Sphere(s.center, s.radius))

realize(b::ACAD, s::Torus) =
  @remote(b, Torus(s.center, vz(1, s.center.cs), s.re, s.ri))

backend_pyramid(b::ACAD, bs::Locs, t::Loc) =
  @remote(b, IrregularPyramid(bs, t))
backend_pyramid_frustum(b::ACAD, bs::Locs, ts::Locs) =
  @remote(b, IrregularPyramidFrustum(bs, ts))

backend_right_cuboid(b::ACAD, cb, width, height, h, material) =
  @remote(b, CenteredBox(cb, width, height, h))
realize(b::ACAD, s::Box) =
  @remote(b, Box(s.c, s.dx, s.dy, s.dz))
realize(b::ACAD, s::Cone) =
  @remote(b, Cone(add_z(s.cb, s.h), s.r, s.cb))
realize(b::ACAD, s::ConeFrustum) =
  @remote(b, ConeFrustum(s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt))
backend_cylinder(b::ACAD, cb::Loc, r::Real, h::Real) =
  @remote(b, Cylinder(cb, r, add_z(cb, h)))

backend_extrusion(b::ACAD, s::Shape, v::Vec) =
    and_mark_deleted(
        map_ref(s) do r
            @remote(b, Extrude(r, v))
        end,
        s)

backend_sweep(b::ACAD, path::Shape, profile::Shape, rotation::Real, scale::Real) =
  and_mark_deleted(
    map_ref(profile) do profile_r
      map_ref(path) do path_r
        @remote(b, Sweep(path_r, profile_r, rotation, scale))
      end
  end, [profile, path])

backend_revolve_point(b::ACAD, profile::Shape, p::Loc, n::Vec, start_angle::Real, amplitude::Real) =
  realize(b, arc(loc_from_o_vz(p, n), distance(profile, p), start_angle, amplitude))
backend_revolve_curve(b::ACAD, profile::Shape, p::Loc, n::Vec, start_angle::Real, amplitude::Real) =
  acad_revolution(b, profile, p, n, start_angle, amplitude)
backend_revolve_surface(b::ACAD, profile::Shape, p::Loc, n::Vec, start_angle::Real, amplitude::Real) =
  acad_revolution(b, profile, p, n, start_angle, amplitude)

acad_revolution(b::ACAD, profile::Shape, p::Loc, n::Vec, start_angle::Real, amplitude::Real) =
  and_delete_shape(
    map_ref(profile) do r
      @remote(b, Revolve(r, p, n, start_angle, amplitude))
    end,
    profile)

backend_loft_curves(b::ACAD, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  and_delete_shapes(@remote(b, Loft(
                             collect_ref(profiles),
                             collect_ref(rails),
                             ruled, closed)),
                    vcat(profiles, rails))

backend_loft_surfaces(b::ACAD, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
    backend_loft_curves(b, profiles, rails, ruled, closed)

backend_loft_curve_point(b::ACAD, profile::Shape, point::Shape) =
    and_delete_shapes(@remote(b, Loft(
                               vcat(collect_ref(profile), collect_ref(point)),
                               [],
                               true, false)),
                      [profile, point])

backend_loft_surface_point(b::ACAD, profile::Shape, point::Shape) =
    backend_loft_curve_point(b, profile, point)

unite_ref(b::ACAD, r0::ACADNativeRef, r1::ACADNativeRef) =
    ensure_ref(b, @remote(b, Unite(r0.value, r1.value)))

intersect_ref(b::ACAD, r0::ACADNativeRef, r1::ACADNativeRef) =
    ensure_ref(b, @remote(b, Intersect(r0.value, r1.value)))

subtract_ref(b::ACAD, r0::ACADNativeRef, r1::ACADNativeRef) =
    ensure_ref(b, @remote(b, Subtract(r0.value, r1.value)))

slice_ref(b::ACAD, r::ACADNativeRef, p::Loc, v::Vec) =
    (@remote(b, Slice(r.value, p, v)); r)

slice_ref(b::ACAD, r::ACADUnionRef, p::Loc, v::Vec) =
    ACADUnionRef(map(r->slice_ref(b, r, p, v), r.values))

unite_refs(b::ACAD, refs::Vector{<:ACADRef}) =
    ACADUnionRef(tuple(refs...))

realize(b::ACAD, s::IntersectionShape) =
  let r = foldl(intersect_ref(b), map(ref, s.shapes),
                init=ACADUniversalRef())
    mark_deleted(s.shapes)
    r
  end

realize(b::ACAD, s::Slice) =
  slice_ref(b, ref(s.shape), s.p, s.n)

realize(b::ACAD, s::Move) =
  let r = map_ref(s.shape) do r
            @remote(b, Move(r, s.v))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Transform) =
  let r = map_ref(s.shape) do r
            @remote(b, Transform(r, s.xform))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Scale) =
  let r = map_ref(s.shape) do r
            @remote(b, Scale(r, s.p, s.s))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Rotate) =
  let r = map_ref(s.shape) do r
            @remote(b, Rotate(r, s.p, s.v, s.angle))
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::ACAD, s::Mirror) =
  and_mark_deleted(map_ref(s.shape) do r
                    @remote(b, Mirror(r, s.p, s.n, false))
                   end,
                   s.shape)

realize(b::ACAD, s::UnionMirror) =
  let r0 = ref(s.shape),
      r1 = map_ref(s.shape) do r
            @remote(b, Mirror(r, s.p, s.n, true))
          end
    UnionRef((r0,r1))
  end

backend_surface_grid(b::ACAD, points, closed_u, closed_v, smooth_u, smooth_v) =
    @remote(b, SurfaceFromGrid(
        size(points,2),
        size(points,1),
        reshape(points,:),
        closed_u,
        closed_v,
        # Autocad does not allow us to distinguish smoothness along different dimensions
        smooth_u && smooth_v ? 2 : 0))

realize(b::ACAD, s::Thicken) =
  and_mark_deleted(
    map_ref(s.shape) do r
      @remote(b, Thicken(r, s.thickness))
    end,
    s.shape)

# backend_frame_at
backend_frame_at(b::ACAD, s::Circle, t::Real) = add_pol(s.center, s.radius, t)

backend_frame_at(b::ACAD, c::Shape1D, t::Real) = @remote(b, CurveFrameAt(ref(c).value, t))
backend_curve_length(b::ACAD, c::Shape1D) = @remote(b, CurveLength(ref(c).value))
backend_frame_at_length(b::ACAD, c::Shape1D, t::Real) = @remote(b, CurveFrameAtLength(ref(c).value, t))

#backend_frame_at(b::ACAD, s::Surface, u::Real, v::Real) =
    #What should we do with v?
#    backend_frame_at(b, s.frontier[1], u)

#backend_frame_at(b::ACAD, s::SurfacePolygon, u::Real, v::Real) =

backend_frame_at(b::ACAD, s::Shape2D, u::Real, v::Real) = @remote(b, SurfaceFrameAt(ref(s).value, u, v))

# BIM
realize(b::ACAD, f::TableFamily) =
    @remote(b, CreateRectangularTableFamily(f.length, f.width, f.height, f.top_thickness, f.leg_thickness))
realize(b::ACAD, f::ChairFamily) =
    @remote(b, CreateChairFamily(f.length, f.width, f.height, f.seat_height, f.thickness))
realize(b::ACAD, f::TableChairFamily) =
    @remote(b, CreateRectangularTableAndChairsFamily(
        realize(b, f.table_family), realize(b, f.chair_family),
        f.table_family.length, f.table_family.width,
        f.chairs_top, f.chairs_bottom, f.chairs_right, f.chairs_left,
        f.spacing))

backend_rectangular_table(b::ACAD, c, angle, family) =
    @remote(b, Table(c, angle, realize(b, family)))

backend_chair(b::ACAD, c, angle, family) =
    @remote(b, Chair(c, angle, realize(b, family)))

backend_rectangular_table_and_chairs(b::ACAD, c, angle, family) =
    @remote(b, TableAndChairs(c, angle, realize(b, family)))

backend_slab(b::ACAD, profile, holes, thickness, family) =
  let slab(profile) = map_ref(b, r -> @remote(b, Extrude(r, vz(thickness))),
                              ensure_ref(b, backend_fill(b, profile))),
      main_body = slab(profile),
      holes_bodies = map(slab, holes)
    foldl((r0, r1)->subtract_ref(b, r0, r1), holes_bodies, init=main_body)
  end

#=
realize_beam_profile(b::ACAD, s::Union{Beam,FreeColumn,Column}, profile::CircularPath, cb::Loc, length::Real) =
  @remote(b, Cylinder(cb, profile.radius, add_z(cb, length)))

realize_beam_profile(b::ACAD, s::Union{Beam,Column}, profile::RectangularPath, cb::Loc, length::Real) =
  let profile_u0 = profile.corner,
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
      # need to test whether it is rotation on center or on axis
      o = loc_from_o_phi(c, s.angle)
    @remote(b, CenteredBox(add_y(o, -profile.dy/2), profile.dx, profile.dy, length))
  end

#Columns are aligned along the center axis.
realize_beam_profile(b::ACAD, s::FreeColumn, profile::RectangularPath, cb::Loc, length::Real) =
  let profile_u0 = profile.corner,
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
      # need to test whether it is rotation on center or on axis
      o = loc_from_o_phi(c, s.angle)
    @remote(b, CenteredBox(o, profile.dx, profile.dy, length))
  end
=#

backend_wall(b::ACAD, path, height, l_thickness, r_thickness, family) =
  @remote(b,
    Thicken(@remote(b,
      Extrude(backend_stroke(b, offset(path, (l_thickness - r_thickness)/2)),
              vz(height))),
      r_thickness + l_thickness))

backend_panel(b::ACAD, bot::Locs, top::Locs, family) =
  @remote(b, IrregularPyramidFrustum(bot, top))

############################################

backend_bounding_box(b::ACAD, shapes::Shapes) =
  @remote(b, BoundingBox(collect_ref(shapes)))

set_view(camera::Loc, target::Loc, lens::Real, aperture::Real, b::ACAD) =
  @remote(b, View(camera, target, lens))

get_view(b::ACAD) =
  @remote(b, ViewCamera()), @remote(b, ViewTarget()), @remote(b, ViewLens())

zoom_extents(b::ACAD) = @remote(b, ZoomExtents())

view_top(b::ACAD) = @remote(b, ViewTop())

backend_realistic_sky(b::ACAD, date, latitude, longitude, meridian, turbidity, withsun) =
  @remote(b, SetSkyFromDateLocation(year(date), month(date), day(date),
                                    hour(date), minute(date),
                                    latitude, longitude, meridian))


backend_delete_shapes(b::ACAD, shapes::Shapes) =
  @remote(b, DeleteMany(collect_ref(shapes)))

delete_all_shapes(b::ACAD) =
  @remote(b, DeleteAll())

set_length_unit(unit::String, b::ACAD) = @remote(b, SetLengthUnit(unit))

# Dimensions

const ACADDimensionStyles = Dict(:architectural => "_ARCHTICK", :mechanical => "")

dimension(p0::Loc, p1::Loc, p::Loc, scale::Real, style::Symbol, b::ACAD=current_backend()) =
    @remote(b, CreateAlignedDimension(p0, p1, p,
        scale,
        ACADDimensionStyles[style]))

dimension(p0::Loc, p1::Loc, sep::Real, scale::Real, style::Symbol, b::ACAD=current_backend()) =
    let v = p1 - p0
        angle = pol_phi(v)
        dimension(p0, p1, add_pol(p0, sep, angle + pi/2), scale, style, b)
    end

# Layers
ACADLayer = Int

current_layer(b::ACAD)::ACADLayer =
  @remote(b, CurrentLayer())

current_layer(layer::ACADLayer, b::ACAD) =
  @remote(b, SetCurrentLayer(layer))

backend_create_layer(b::ACAD, name::String, active::Bool, color::RGB) =
  let to255(x) = round(UInt8, x*255)
    @remote(b, CreateLayer(name, true, to255(red(color)), to255(green(color)), to255(blue(color))))
  end

delete_all_shapes_in_layer(layer::ACADLayer, b::ACAD) =
  @remote(b, DeleteAllInLayer(layer))

switch_to_layer(to, b::ACAD) =
    if to != from
      set_layer_active(to, true)
      set_layer_active(from, false)
      current_layer(to)
    end

# Materials
ACADMaterial = Int

current_material(b::ACAD)::ACADMaterial =
  -1 #@remote(b, CurrentMaterial())

current_material(material::ACADMaterial, b::ACAD) =
  -1 #@remote(b, SetCurrentMaterial(material))

get_material(name::String, b::ACAD) =
  -1 #@remote(b, CreateMaterial(name))

# Blocks

realize(b::ACAD, s::Block) =
    @remote(b, CreateBlockFromShapes(s.name, collect_ref(s.shapes)))

realize(b::ACAD, s::BlockInstance) =
    @remote(b, CreateBlockInstance(
        collect_ref(s.block)[1],
        center_scaled_cs(s.loc, s.scale, s.scale, s.scale)))

#=

# Manual process
@time for i in 1:1000 for r in 1:10 circle(x(i*10), r) end end

# Create block...
Khepri.create_block("Foo", [circle(radius=r) for r in 1:10])

# ...and instantiate it
@time for i in 1:1000 Khepri.instantiate_block("Foo", x(i*10), 0) end

=#

# Lights
backend_pointlight(b::ACAD, loc::Loc, color::RGB, range::Real, intensity::Real) =
  # HACK: Fix this
  @remote(b, SpotLight(loc, intensity, range, loc+vz(-1)))

backend_spotlight(b::ACAD, loc::Loc, dir::Vec, hotspot::Real, falloff::Real) =
    @remote(b, SpotLight(loc, hotspot, falloff, loc + dir))

backend_ieslight(b::ACAD, file::String, loc::Loc, dir::Vec, alpha::Real, beta::Real, gamma::Real) =
    @remote(b, IESLight(file, loc, loc + dir, vxyz(alpha, beta, gamma)))

# User Selection

shape_from_ref(r, b::ACAD) =
    let c = connection(b),
        code = @remote(b, ShapeCode(r)),
        ref = LazyRef(b, ACADNativeRef(r))
        if code == 1 # Point
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
            let tans = @remote(b, SplineTangents(r))
                if length(tans[1]) < 1e-20 && length(tans[2]) < 1e-20
                    closed_spline(@remote(b, SplineInterpPoints(r))[1:end-1],
                                  backend=b, ref=ref)
                else
                    spline(@remote(b, SplineInterpPoints(r)), tans[1], tans[2],
                           backend=b, ref=ref)
                end
            end
        elseif code == 9
            let start_angle = mod(@remote(b, ArcStartAngle(r)), 2pi),
                end_angle = mod(@remote(b, ArcEndAngle(r)), 2pi)
                arc(maybe_loc_from_o_vz(@remote(b, ArcCenter(r)), @remote(b, ArcNormal(r))),
                    @remote(b, ArcRadius(r)), start_angle, mod(end_angle - start_angle, 2pi),
                    backend=b, ref=ref)
            #=    if end_angle > start_angle
                    arc(maybe_loc_from_o_vz(@remote(b, ArcCenter(r)), @remote(b, ArcNormal(r))),
                        @remote(b, ArcRadius(r)), start_angle, end_angle - start_angle,
                        backend=b, ref=ref)
                else
                    arc(maybe_loc_from_o_vz(@remote(b, ArcCenter(r)), @remote(b, ArcNormal(r))),
                        @remote(b, ArcRadius(r)), end_angle, start_angle - end_angle,
                        backend=b, ref=ref)
                end=#
            end
        elseif code == 10
            let str = @remote(b, TextString(r)),
                height = @remote(b, TextHeight(r)),
                loc = @remote(b, TextPosition(r))
                text(str, loc, height, backend=b, ref=ref)
            end
        elseif code == 11
            let str = @remote(b, MTextString(r)),
                height = @remote(b, MTextHeight(r)),
                loc = @remote(b, MTextPosition(r))
                text(str, loc, height, backend=b, ref=ref)
            end
        elseif code == 16
            let pts = @remote(b, MeshVertices(r)),
                (type, n, m, n_closed, m_closed) = @remote(b, PolygonMeshData(r))
                surface_grid(reshape(pts, (n, m)), n_closed == 1, m_closed == 1, ref=ref)
            end
        elseif 12 <= code <= 14
            surface(Shapes1D[], backend=b, ref=ref)
        elseif 103 <= code <= 106
            polygon(@remote(b, LineVertices(r)),
                    backend=b, ref=ref)
        elseif code == 107
            closed_spline(@remote(b, SplineInterpPoints(r))[1:end-1],
                          backend=b, ref=ref)
        else
            #unknown(backend=b, ref=ref)
            unknown(r, backend=b, ref=LazyRef(b, ACADNativeRef(r), 0, 0)) # To force copy
            #error("Unknown shape with code $(code)")
        end
    end
#

#=
In case we need to realize an Unknown shape, we just copy it
=#

realize(b::ACAD, s::Unknown) =
    @remote(b, Copy(s.baseref))



select_position(prompt::String, b::ACAD) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = @remote(b, GetPosition(prompt))
      length(ans) > 0 ? ans[1] : nothing
    end
  end

select_positions(prompt::String, b::ACAD) =
  let sel() =
    let p = select_position(prompt, b)
      if p == nothing
        []
      else
        [p, sel()...]
      end
    end
    sel()
  end


# HACK: The next operations should receive a set of shapes to avoid re-creating already existing shapes

select_point(prompt::String, b::ACAD) =
  select_one_with_prompt(prompt, b, @get_remote b GetPoint)

select_points(prompt::String, b::ACAD) =
  select_many_with_prompt(prompt, b, @get_remote b GetPoints)

select_curve(prompt::String, b::ACAD) =
  select_one_with_prompt(prompt, b, @get_remote b GetCurve)

select_curves(prompt::String, b::ACAD) =
  select_many_with_prompt(prompt, b, @get_remote b GetCurves)


select_surface(prompt::String, b::ACAD) =
  select_one_with_prompt(prompt, b, @get_remote b GetSurface)

select_surfaces(prompt::String, b::ACAD) =
  select_many_with_prompt(prompt, b, @get_remote b GetSurfaces)


select_solid(prompt::String, b::ACAD) =
  select_one_with_prompt(prompt, b, @get_remote b GetSolid)

select_solids(prompt::String, b::ACAD) =
  select_many_with_prompt(prompt, b, @get_remote b GetSolids)


select_shape(prompt::String, b::ACAD) =
  select_one_with_prompt(prompt, b, @get_remote b GetShape)

select_shapes(prompt::String, b::ACAD) =
  select_many_with_prompt(prompt, b, @get_remote b GetShapes)


captured_shape(b::ACAD, handle) =
  shape_from_ref(@remote(b, GetShapeFromHandle(handle)), b)
#
captured_shapes(b::ACAD, handles) =
  map(handles) do handle
      shape_from_ref(@remote(b, GetShapeFromHandle(handle)), b)
  end

generate_captured_shape(s::Shape, b::ACAD) =
    println("captured_shape(autocad, $(@remote(b, GetHandleFromShape(ref(s).value))))")

generate_captured_shapes(ss::Shapes, b::ACAD) =
  begin
    print("captured_shapes(autocad, [")
    for s in ss
      print(@remote(b, GetHandleFromShape(ref(s).value)))
      print(", ")
    end
    println("])")
  end

# Register for notification


register_for_changes(s::Shape, b::ACAD) =
    let conn = connection(b)
        @remote(b, RegisterForChanges(ref(s).value))
        @remote(b, DetectCancel())
        s
    end

unregister_for_changes(s::Shape, b::ACAD) =
    let conn = connection(b)
        @remote(b, UnregisterForChanges(ref(s).value))
        @remote(b, UndetectCancel())
        s
    end

waiting_for_changes(s::Shape, b::ACAD) =
    ! @remote(b, WasCanceled())

changed_shape(ss::Shapes, b::ACAD) =
    let conn = connection(b)
        changed = []
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
all_shapes(b::ACAD) =
  Shape[shape_from_ref(r, b)
        for r in filter(r -> @remote(b, ShapeCode(r)) != 0, @remote(b, GetAllShapes()))]

all_shapes_in_layer(layer, b::ACAD) =
  Shape[shape_from_ref(r, b) for r in @remote(b, GetAllShapesInLayer(layer))]

highlight_shape(s::Shape, b::ACAD) =
  @remote(b, SelectShapes(collect_ref(s)))

highlight_shapes(ss::Shapes, b::ACAD) =
  @remote(b, SelectShapes(collect_ref(ss)))

pre_selected_shapes_from_set(ss::Shapes) =
  length(ss) == 0 ? [] : pre_selected_shapes_from_set(ss, backend(ss[1]))

# HACK: This must be implemented for all backends
pre_selected_shapes_from_set(ss::Shapes, b::Backend) = []

pre_selected_shapes_from_set(ss::Shapes, b::ACAD) =
  let refs = map(id -> @remote(b, GetHandleFromShape(id)), @remote(b, GetPreSelectedShapes()))
    filter(s -> @remote(b, GetHandleFromShape(ref(s).value)) in refs, ss)
  end
disable_update(b::ACAD) =
  @remote(b, DisableUpdate())

enable_update(b::ACAD) =
  @remote(b, EnableUpdate())

# Render

#render exposure: [-3, +3] -> [-6, 21]
convert_render_exposure(b::ACAD, v::Real) = -4.05*v + 8.8
#render quality: [-1, +1] -> [+1, +50]
convert_render_quality(b::ACAD, v::Real) = round(Int, 25.5 + 24.5*v)

render_view(path::String, b::ACAD) =
    @remote(b, Render(
               render_width(), render_height(),
               path,
               convert_render_quality(b, render_quality()),
               convert_render_exposure(b, render_exposure())))

export mentalray_render_view
mentalray_render_view(name::String) =
    let b = current_backend()
        @remote(b, SetSystemVariableInt("SKYSTATUS", 2)) # skystatus:background-and-illumination
        @remote(b, Command("._-render P _R $(render_width()) $(render_height()) _yes $(prepare_for_saving_file(render_pathname(name)))\n"))
    end

save_as(pathname::String, format::String, b::ACAD) =
    @remote(b, SaveAs(pathname, format))


export autocad_command
autocad_command(s::String) =
  @remote(autocad, Command("$(s)\n"))
