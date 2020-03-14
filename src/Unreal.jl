export unreal, fast_unreal,
       unreal_material_family

# We need some additional Encoders
encode_Actor = encode_int
decode_Actor = decode_int_or_error
encode_Actor_array = encode_int_array
decode_Actor_array = decode_int_array

abstract type UnrealKey end
const UnrealId = Int
const UnrealIds = Vector{UnrealId}
const UnrealRef = GenericRef{UnrealKey, UnrealId}
const UnrealRefs = Vector{UnrealRef}
const UnrealEmptyRef = EmptyRef{UnrealKey, UnrealId}
const UnrealUniversalRef = UniversalRef{UnrealKey, UnrealId}
const UnrealNativeRef = NativeRef{UnrealKey, UnrealId}
const UnrealUnionRef = UnionRef{UnrealKey, UnrealId}
const UnrealSubtractionRef = SubtractionRef{UnrealKey, UnrealId}
const Unreal = SocketBackend{UnrealKey, UnrealId}

void_ref(b::Unreal) = UnrealNativeRef(-1)

create_Unreal_connection() =
    begin
        #check_plugin()
        create_backend_connection("Unreal", 11002)
    end

unreal_functions = @remote_functions :CPP """
  public Actor Sphere(Vector3 center, float radius)
"""

const unreal = Unreal(LazyParameter(TCPSocket, create_Unreal_connection),
                      unreal_functions)

backend_name(b::Unreal) = "Unreal"

realize(b::Unreal, s::EmptyShape) =
  UnrealEmptyRef()
realize(b::Unreal, s::UniversalShape) =
  UnrealUniversalRef()

#=
unreal"public void SetApplyMaterials(bool apply)"
unreal"public void SetApplyColliders(bool apply)"

#=
realize(b::Unreal, s::Point) =
  UnrealPoint(connection(b), s.position)
realize(b::Unreal, s::Line) =
  UnrealPolyLine(connection(b), s.vertices)
realize(b::Unreal, s::Spline) =
  if (s.v0 == false) && (s.v1 == false)
    #UnrealSpline(connection(b), s.points)
    UnrealInterpSpline(connection(b),
                     s.points,
                     s.points[2]-s.points[1],
                     s.points[end]-s.points[end-1])
  elseif (s.v0 != false) && (s.v1 != false)
    UnrealInterpSpline(connection(b), s.points, s.v0, s.v1)
  else
    UnrealInterpSpline(connection(b),
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1)
  end
realize(b::Unreal, s::ClosedSpline) =
  UnrealInterpClosedSpline(connection(b), s.points)
realize(b::Unreal, s::Circle) =
  UnrealCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
realize(b::Unreal, s::Arc) =
  if s.radius == 0
    UnrealPoint(connection(b), s.center)
  elseif s.amplitude == 0
    UnrealPoint(connection(b), s.center + vpol(s.radius, s.start_angle, s.center.cs))
  elseif abs(s.amplitude) >= 2*pi
    UnrealCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
  else
    end_angle = s.start_angle + s.amplitude
    if end_angle > s.start_angle
      UnrealArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle)
    else
      UnrealArc(connection(b), s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle)
    end
  end

realize(b::Unreal, s::Ellipse) =
  if s.radius_x > s.radius_y
    UnrealEllipse(connection(b), s.center, vz(1, s.center.cs), vxyz(s.radius_x, 0, 0, s.center.cs), s.radius_y/s.radius_x)
  else
    UnrealEllipse(connection(b), s.center, vz(1, s.center.cs), vxyz(0, s.radius_y, 0, s.center.cs), s.radius_x/s.radius_y)
  end
realize(b::Unreal, s::EllipticArc) =
  error("Finish this")

realize(b::Unreal, s::Polygon) =
  UnrealClosedPolyLine(connection(b), s.vertices)
realize(b::Unreal, s::RegularPolygon) =
  UnrealClosedPolyLine(connection(b), regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
realize(b::Unreal, s::Rectangle) =
  UnrealClosedPolyLine(
    connection(b),
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)])
realize(b::Unreal, s::SurfaceCircle) =
  UnrealSurfaceCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
realize(b::Unreal, s::SurfaceArc) =
    #UnrealSurfaceArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, s.start_angle + s.amplitude)
    if s.radius == 0
        UnrealPoint(connection(b), s.center)
    elseif s.amplitude == 0
        UnrealPoint(connection(b), s.center + vpol(s.radius, s.start_angle, s.center.cs))
    elseif abs(s.amplitude) >= 2*pi
        UnrealSurfaceCircle(connection(b), s.center, vz(1, s.center.cs), s.radius)
    else
        end_angle = s.start_angle + s.amplitude
        if end_angle > s.start_angle
            UnrealSurfaceFromCurves(connection(b),
                [UnrealArc(connection(b), s.center, vz(1, s.center.cs), s.radius, s.start_angle, end_angle),
                 UnrealPolyLine(connection(b), [add_pol(s.center, s.radius, end_angle),
                                              add_pol(s.center, s.radius, s.start_angle)])])
        else
            UnrealSurfaceFromCurves(connection(b),
                [UnrealArc(connection(b), s.center, vz(1, s.center.cs), s.radius, end_angle, s.start_angle),
                 UnrealPolyLine(connection(b), [add_pol(s.center, s.radius, s.start_angle),
                                              add_pol(s.center, s.radius, end_angle)])])
        end
    end

#realize(b::Unreal, s::SurfaceElliptic_Arc) = UnrealCircle(connection(b),
#realize(b::Unreal, s::SurfaceEllipse) = UnrealCircle(connection(b),
=#

unreal"public Actor SurfacePolygon(Vector3[] ps)"

realize(b::Unreal, s::SurfacePolygon) =
  UnrealSurfacePolygon(connection(b), reverse(s.vertices))

backend_fill(b::Unreal, path::ClosedPolygonalPath) =
  UnrealSurfacePolygon(connection(b), path.vertices)

#=
realize(b::Unreal, s::SurfaceRegularPolygon) =
  UnrealSurfaceClosedPolyLine(connection(b), regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed))
realize(b::Unreal, s::SurfaceRectangle) =
  UnrealSurfaceClosedPolyLine(
    connection(b),
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)])
realize(b::Unreal, s::Surface) =
  let #ids = map(r->UnrealNurbSurfaceFrom(connection(b),r), UnrealSurfaceFromCurves(connection(b), collect_ref(s.frontier)))
      ids = UnrealSurfaceFromCurves(connection(b), collect_ref(s.frontier))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::Unreal, s::Shape2D) =
    map(shape_from_ref, UnrealCurvesFromSurface(connection(b), ref(s).value))
=#
backend_fill(b::Unreal, path::ClosedPathSequence) =
  backend_fill(b, convert(ClosedPolygonalPath, path))

#=
# Iterating over curves and surfaces

Unreal"public double[] CurveDomain(Entity ent)"
Unreal"public double CurveLength(Entity ent)"
Unreal"public Frame3d CurveFrameAt(Entity ent, double t)"
Unreal"public Frame3d CurveFrameAtLength(Entity ent, double l)"
=#

backend_map_division(b::Unreal, f::Function, s::Shape1D, n::Int) =
  let (t1, t2) = curve_domain(s)
    map_division(t1, t2, n) do t
      f(frame_at(s, t))
    end
  end
#=
Unreal"public Vector3d RegionNormal(Entity ent)"
Unreal"public Point3d RegionCentroid(Entity ent)"
Unreal"public double[] SurfaceDomain(Entity ent)"
Unreal"public Frame3d SurfaceFrameAt(Entity ent, double u, double v)"

backend_surface_domain(b::Unreal, s::Shape2D) =
    tuple(UnrealSurfaceDomain(connection(b), ref(s).value)...)

backend_map_division(b::Unreal, f::Function, s::Shape2D, nu::Int, nv::Int) =
    let conn = connection(b)
        r = ref(s).value
        (u1, u2, v1, v2) = UnrealSurfaceDomain(conn, r)
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(UnrealSurfaceFrameAt(conn, r, u, v))
            end
        end
    end

# The previous method cannot be applied to meshes in AutoCAD, which are created by surface_grid

backend_map_division(b::Unreal, f::Function, s::SurfaceGrid, nu::Int, nv::Int) =
    let (u1, u2, v1, v2) = UnrealSurfaceDomain(conn, r)
        map_division(u1, u2, nu) do u
            map_division(v1, v2, nv) do v
                f(UnrealSurfaceFrameAt(conn, r, u, v))
            end
        end
    end
=#

unreal"public Actor Text(string txt, Vector3 position, Vector3 vx, Vector3 vy, string fontName, int fontSize)"

realize(b::Unreal, s::Text) =
  UnrealText(
    connection(b),
    s.str, s.corner, vz(-1, s.corner.cs), vy(1, s.corner.cs), "Fonts/Inconsolata-Regular", s.height)

=#

realize(b::Unreal, s::Sphere) =
  @remote(b, Sphere(s.center, s.radius))

#=
#=
realize(b::Unreal, s::Torus) =
  UnrealTorus(connection(b), s.center, vz(1, s.center.cs), s.re, s.ri)
=#

unreal"public Actor Pyramid(Vector3[] ps, Vector3 q)"
unreal"public Actor PyramidFrustum(Vector3[] ps, Vector3[] qs)"

realize(b::Unreal, s::Cuboid) =
  UnrealPyramidFrustum(connection(b), [s.b0, s.b1, s.b2, s.b3], [s.t0, s.t1, s.t2, s.t3])

realize(b::Unreal, s::RegularPyramidFrustum) =
    UnrealPyramidFrustum(connection(b),
                        regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
                        regular_polygon_vertices(s.edges, add_z(s.cb, s.h), s.rt, s.angle, s.inscribed))

realize(b::Unreal, s::RegularPyramid) =
  UnrealPyramid(connection(b),
               regular_polygon_vertices(s.edges, s.cb, s.rb, s.angle, s.inscribed),
               add_z(s.cb, s.h))

realize(b::Unreal, s::IrregularPyramid) =
  UnrealPyramid(connection(b), s.bs, s.t)

realize(b::Unreal, s::RegularPrism) =
  let bs = regular_polygon_vertices(s.edges, s.cb, s.r, s.angle, s.inscribed)
    UnrealPyramidFrustum(connection(b),
                        bs,
                        map(p -> add_z(p, s.h), bs))
  end

realize(b::Unreal, s::IrregularPyramidFrustum) =
    UnrealPyramidFrustum(connection(b), s.bs, s.ts)

realize(b::Unreal, s::IrregularPrism) =
  UnrealPyramidFrustum(connection(b),
                      s.bs,
                      map(p -> (p + s.v), s.bs))

unreal"public Actor RightCuboid(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz, float angle)"

realize(b::Unreal, s::RightCuboid) =
  UnrealRightCuboid(connection(b), s.cb, vz(1, s.cb.cs), vx(1, s.cb.cs), s.height, s.width, s.h, s.angle)

unreal"public Actor Box(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz)"
#unreal"public Actor Box2(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz)"

realize(b::Unreal, s::Box) =
  UnrealBox(connection(b), s.c, vz(1, s.c.cs), vx(1, s.c.cs), s.dy, s.dx, s.dz)

realize(b::Unreal, s::Cone) =
  UnrealPyramid(connection(b), regular_polygon_vertices(64, s.cb, s.r), add_z(s.cb, s.h))

realize(b::Unreal, s::ConeFrustum) =
  UnrealPyramidFrustum(connection(b),
    regular_polygon_vertices(64, s.cb, s.rb),
    regular_polygon_vertices(64, s.cb + vz(s.h, s.cb.cs), s.rt))

unreal"public Actor Cylinder(Vector3 bottom, float radius, Vector3 top)"

realize(b::Unreal, s::Cylinder) =
  UnrealCylinder(connection(b), s.cb, s.r, s.cb + vz(s.h, s.cb.cs))

#=
backend_extrusion(b::Unreal, s::Shape, v::Vec) =
    and_mark_deleted(
        map_ref(s) do r
            UnrealExtrude(connection(b), r, v)
        end,
        s)

backend_sweep(b::Unreal, path::Shape, profile::Shape, rotation::Real, scale::Real) =
  map_ref(profile) do profile_r
    map_ref(path) do path_r
      UnrealSweep(connection(b), path_r, profile_r, rotation, scale)
    end
  end

realize(b::Unreal, s::Revolve) =
  and_delete_shape(
    map_ref(s.profile) do r
      UnrealRevolve(connection(b), r, s.p, s.n, s.start_angle, s.amplitude)
    end,
    s.profile)

backend_loft_curves(b::Unreal, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
  and_delete_shapes(UnrealLoft(connection(b),
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


backend_loft_surfaces(b::Unreal, profiles::Shapes, rails::Shapes, ruled::Bool, closed::Bool) =
    backend_loft_curves(b, profiles, rails, ruled, closed)

backend_loft_curve_point(b::Unreal, profile::Shape, point::Shape) =
    and_delete_shapes(UnrealLoft(connection(b),
                               vcat(collect_ref(profile), collect_ref(point)),
                               [],
                               true, false),
                      [profile, point])

backend_loft_surface_point(b::Unreal, profile::Shape, point::Shape) =
    backend_loft_curve_point(b, profile, point)

=#
unreal"public Actor Unite(Actor s0, Actor s1)"
unreal"public Actor Intersect(Actor s0, Actor s1)"
unreal"public Actor Subtract(Actor s0, Actor s1)"
unreal"public void SubtractFrom(Actor s0, Actor s1)"

unreal"public Actor Canonicalize(Actor s)"
###
unite_ref(b::Unreal, r0::UnrealNativeRef, r1::UnrealNativeRef) =
    ensure_ref(b, UnrealUnite(connection(b), r0.value, r1.value))

intersect_ref(b::Unreal, r0::UnrealNativeRef, r1::UnrealNativeRef) =
    ensure_ref(b, UnrealIntersect(connection(b), r0.value, r1.value))

subtract_ref(b::Unreal, r0::UnrealNativeRef, r1::UnrealNativeRef) =
    let r = UnrealSubtract(connection(b), r0.value, r1.value)
      UnrealDeleteMany(connection(b), [r0.value, r1.value])
      r
    end

#=
subtract_ref(b::Unreal, r0::UnrealNativeRef, r1::UnrealNativeRef) =
    begin
      UnrealSubtractFrom(connection(b), r0.value, r1.value)
      r0.value
    end
=#

#=
slice_ref(b::Unreal, r::UnrealNativeRef, p::Loc, v::Vec) =
    (UnrealSlice(connection(b), r.value, p, v); r)

slice_ref(b::Unreal, r::UnrealUnionRef, p::Loc, v::Vec) =
    map(r->slice_ref(b, r, p, v), r.values)

=#
unite_refs(b::Unreal, refs::Vector{<:UnrealRef}) =
    UnrealUnionRef(tuple(refs...))

#
realize(b::Unreal, s::UnionShape) =
  let r = foldl((r0,r1)->unite_ref(b,r0,r1), map(ref, s.shapes),
                init=UnrealEmptyRef())
    delete_shapes(s.shapes)
    #UnrealCanonicalize(connection(b), r.value)
    r
  end

realize(b::Unreal, s::IntersectionShape) =
  let r = foldl((r0,r1)->intersect_ref(b,r0,r1), map(ref, s.shapes),
                init=UnrealUniversalRef())
    delete_shapes(s.shapes)
    r
  end

realize(b::Unreal, s::Slice) =
  slice_ref(b, ref(s.shape), s.p, s.n)

unreal"public void Move(Actor s, Vector3 v)"
unreal"public void Scale(Actor s, Vector3 p, float scale)"
unreal"public void Rotate(Actor s, Vector3 p, Vector3 n, float a)"

realize(b::Unreal, s::Move) =
  let r = map_ref(s.shape) do r
            UnrealMove(connection(b), r, s.v)
            r
          end
    mark_deleted(s.shape)
    r
  end
#=
realize(b::Unreal, s::Transform) =
  let r = map_ref(s.shape) do r
            UnrealTransform(connection(b), r, s.xform)
            r
          end
    mark_deleted(s.shape)
    r
  end
=#
realize(b::Unreal, s::Scale) =
  let r = map_ref(s.shape) do r
            UnrealScale(connection(b), r, s.p, s.s)
            r
          end
    mark_deleted(s.shape)
    r
  end

realize(b::Unreal, s::Rotate) =
  let r = map_ref(s.shape) do r
            UnrealRotate(connection(b), r, s.p, s.v, s.angle)
            r
          end
    mark_deleted(s.shape)
    r
  end

#=
realize(b::Unreal, s::Mirror) =
  and_delete_shape(map_ref(s.shape) do r
                    UnrealMirror(connection(b), r, s.p, s.n, false)
                   end,
                   s.shape)

realize(b::Unreal, s::UnionMirror) =
  let r0 = ref(s.shape),
      r1 = map_ref(s.shape) do r
            UnrealMirror(connection(b), r, s.p, s.n, true)
          end
    UnionRef((r0,r1))
  end
=#
unreal"public Actor SurfaceFromGrid(int m, int n, Vector3[] pts, bool closedM, bool closedN, int level)"

realize(b::Unreal, s::SurfaceGrid) =
    UnrealSurfaceFromGrid(
        connection(b),
        size(s.points,1),
        size(s.points,2),
        reshape(s.points,:),
        s.closed_u,
        s.closed_v,
        2)
#=
realize(b::Unreal, s::Thicken) =
  and_delete_shape(
    map_ref(s.shape) do r
      UnrealThicken(connection(b), r, s.thickness)
    end,
    s.shape)

# backend_frame_at
backend_frame_at(b::Unreal, s::Circle, t::Real) = add_pol(s.center, s.radius, t)

backend_frame_at(b::Unreal, c::Shape1D, t::Real) = UnrealCurveFrameAt(connection(b), ref(c).value, t)

#backend_frame_at(b::Unreal, s::Surface, u::Real, v::Real) =
    #What should we do with v?
#    backend_frame_at(b, s.frontier[1], u)

#backend_frame_at(b::Unreal, s::SurfacePolygon, u::Real, v::Real) =

backend_frame_at(b::Unreal, s::Shape2D, u::Real, v::Real) = UnrealSurfaceFrameAt(connection(b), ref(s).value, u, v)

=#

# BIM

unreal"public Actor LoadResource(String name)"
unreal"public Material LoadMaterial(String name)"
unreal"public void SetCurrentMaterial(Material material)"

unreal"public Actor InstantiateResource(Actor family, Vector3 pos, Vector3 vx, Vector3 vy, float scale)"
unreal"public Actor InstantiateBIMElement(Actor family, Vector3 pos, float angle)"

# Families

abstract type UnrealFamily <: Family end

struct UnrealMaterialFamily <: UnrealFamily
  name::String
  parameter_map::Dict{Symbol,String}
  ref::Parameter{Any}
end

unreal_material_family(name, pairs...) = UnrealMaterialFamily(name, Dict(pairs...), Parameter{Any}(nothing))
backend_get_family_ref(b::Unreal, f::Family, uf::UnrealMaterialFamily) = UnrealLoadMaterial(connection(b), uf.name)

struct UnrealResourceFamily <: UnrealFamily
  name::String
  parameter_map::Dict{Symbol,String}
  ref::Parameter{Any}
end

unreal_resource_family(name, pairs...) = UnrealResourceFamily(name, Dict(pairs...), Parameter{Any}(nothing))
backend_get_family_ref(b::Unreal, f::Family, uf::UnrealResourceFamily) = UnrealLoadResource(connection(b), uf.name)


set_backend_family(default_wall_family(), unreal, unreal_material_family("Materials/Plaster/Plaster1"))
set_backend_family(default_slab_family(), unreal, unreal_material_family("Materials/Concrete/Concrete2"))
set_backend_family(default_roof_family(), unreal, unreal_material_family("Materials/Concrete/Concrete2"))
set_backend_family(default_beam_family(), unreal, unreal_material_family("Materials/Metal/Aluminum"))
set_backend_family(default_column_family(), unreal, unreal_material_family("Materials/Concrete/Concrete2"))
set_backend_family(default_door_family(), unreal, unreal_material_family("Materials/Wood/InteriorWood2"))
set_backend_family(default_panel_family(), unreal, unreal_material_family("Materials/Glass/Glass"))

unreal"public Actor Window(Vector3 position, Quaternion rotation, float dx, float dy, float dz)"
unreal"public Actor Shelf(Vector3 position, int rowLength, int lineLength, float cellWidth, float cellHeight, float cellDepth)"


set_backend_family(default_table_family(), unreal, unreal_resource_family("Prefabs/Tables/ModernTable/ModernTable"))
set_backend_family(default_chair_family(), unreal, unreal_resource_family("Prefabs/Chairs/ModernChair/ModernChair"))
set_backend_family(default_table_chair_family(), unreal, unreal_resource_family("Prefabs/TablesChairs/ModernTableChair/ModernTableChair"))

backend_rectangular_table(b::Unreal, c, angle, family) =
    UnrealInstantiateBIMElement(connection(b), realize(b, family), c, -angle)

backend_chair(b::Unreal, c, angle, family) =
    UnrealInstantiateBIMElement(connection(b), realize(b, family), c, -angle)

backend_rectangular_table_and_chairs(b::Khepri.Unreal, c, angle, family) =
    UnrealInstantiateBIMElement(connection(b), realize(b, family), c, -angle)

unreal"public Actor Slab(Vector3[] contour, Vector3[][] holes, float h, Material material)"

backend_slab(b::Unreal, profile, holes, thickness, family) =
  let bot_vs = path_vertices(profile)
      c = connection(b)
    UnrealSlab(c, bot_vs, map(path_vertices, holes), thickness, realize(b, family))
  end

unreal"public Actor BeamRectSection(Vector3 position, Vector3 vx, Vector3 vy, float dx, float dy, float dz, float angle, Material material)"
unreal"public Actor BeamCircSection(Vector3 bot, float radius, Vector3 top, Material material)"

realize(b::Unreal, s::Beam) =
  let profile = s.family.profile
      profile_u0 = profile.corner
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
    UnrealBeamRectSection(
      connection(b),
      c, vz(1, c.cs), vx(1, c.cs),
      profile.dy, profile.dx, s.h, -s.angle,
      realize(b, s.family))
  end

#Columns are aligned along the center axis.
realize(b::Unreal, s::FreeColumn) =
  let profile = s.family.profile
      profile_u0 = profile.corner
      c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2)
    UnrealBeamRectSection(
      connection(b),
      c, vz(1, c.cs), vx(1, c.cs),
      profile.dy, profile.dx, s.h, -s.angle,
      realize(b, s.family))
  end

realize(b::Unreal, s::Column) =
    let profile = s.family.profile,
        profile_u0 = profile.corner,
        c = add_xy(s.cb, profile_u0.x + profile.dx/2, profile_u0.y + profile.dy/2),
        base_height = s.bottom_level.height,
        height = s.top_level.height - base_height
        # need to test whether it is rotation on center or on axis
      UnrealBeamRectSection(connection(b),
        c, vz(1, c.cs), vx(1, c.cs), profile.dy, profile.dx, height, -s.angle,
        realize(b, s.family))
    end

unreal"public Actor Panel(Vector3[] pts, Vector3 n, Material material)"

realize(b::Unreal, s::Panel) =
  let #p1 = s.vertices[1],
      #p2 = s.vertices[2],
      #p3 = s.vertices[3],
      #n = vz(s.family.thickness, cs_from_o_vx_vy(p1, p2-p1, p3-p1))
      verts = in_world.(s.vertices)
      n = vertices_normal(verts)*(s.family.thickness/2)
    UnrealPanel(
      connection(b),
      map(p -> p - n, verts),
      n*2,
      realize(b, s.family))
  end

sweep_fractions(b, verts, height, thickness) =
  let p = add_z(verts[1], height/2),
      q = add_z(verts[2], height/2),
      (c, h) = position_and_height(p, q),
      s = UnrealNativeRef(UnrealRightCuboid(connection(b), c, vz(1, c.cs), vx(1, c.cs), height, thickness, h, 0))
    if length(verts) > 2
      (s, sweep_fractions(b, verts[2:end], height, thickness)...)
    else
      (s, )
    end
  end

backend_wall(b::Unreal, path, height, l_thickness, r_thickness, family) =
  path_length(path) < path_tolerance() ?
    UnrealEmptyRef() :
    let c = connection(b)
      UnrealSetCurrentMaterial(c, realize(b, family))
      backend_wall_path(
          b,
          path,
          height*0.999, #We reduce height just a bit to avoid Z-fighting
          l_thickness, r_thickness)
    end

backend_wall_path(b::Unreal, path::PolygonalPath, height, l_thickness, r_thickness) =
    UnrealUnionRef(sweep_fractions(b, path.vertices, height, l_thickness, r_thickness))

backend_wall_path(b::Unreal, path::RectangularPath, height, l_thickness, r_thickness) =
    backend_wall_path(b, convert(OpenPolygonalPath, path), height, l_thickness, r_thickness)


#=
realize(b::Radiance, w::Wall) =
  let w_base_height = w.bottom_level.height,
      w_height = w.top_level.height - w_base_height,
      r_thickness = r_thickness(w),
      l_thickness = l_thickness(w),
      w_path = translate(w.path, vz(w_base_height)),
      w_paths = subpaths(w_path),
      r_w_paths = subpaths(offset(w_path, r_thickness)),
      l_w_paths = subpaths(offset(w_path, l_thickness)),
      openings = [w.doors..., w.windows...],
      prevlength = 0
    for (w_seg_path, r_w_path, l_w_path) in zip(w_paths, r_w_paths, l_w_paths)
      let currlength = prevlength + path_length(w_seg_path),
          c_r_w_path = closed_path_for_height(r_w_path, w_height),
          c_l_w_path = closed_path_for_height(l_w_path, w_height)
        realize_pyramid_frustum(b, w, "wall", c_l_w_path, c_r_w_path, false)
        openings = filter(openings) do op
          if prevlength <= op.loc.x < currlength ||
             prevlength <= op.loc.x + op.family.width <= currlength # contained (at least, partially)
            let op_height = op.family.height,
                op_at_start = op.loc.x <= prevlength,
                op_at_end = op.loc.x + op.family.width >= currlength,
                op_path = subpath(w_path,
                                  max(prevlength, op.loc.x),
                                  min(currlength, op.loc.x + op.family.width)),
                r_op_path = offset(op_path, r_thickness),
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
              realize_pyramid_frustum(b, w, "wall", c_r_op_path, c_l_op_path, false)
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
        realize_polygon(b, w, "wall", c_l_w_path, false)
        realize_polygon(b, w, "wall", c_r_w_path, true)
      end
    end
    void_ref(b)
  end
=#

set_backend_family(default_curtain_wall_family().panel,
  unreal,
  unreal_material_family("Materials/Glass/GlassBlue"))
set_backend_family(default_curtain_wall_family().boundary_frame,
  unreal,
  unreal_material_family("Materials/Metal/Steel"))
set_backend_family(default_curtain_wall_family().transom_frame,
  unreal,
  unreal_material_family("Materials/Metal/Steel"))
set_backend_family(default_curtain_wall_family().mullion_frame,
  unreal,
  unreal_material_family("Materials/Metal/Steel"))

backend_curtain_wall(b::Unreal, s, path::Path, bottom::Real, height::Real, r_thickness::Real, l_thickness::Real, kind::Symbol) =
  backend_wall(b, translate(path, vz(bottom)), height, l_thickness, r_thickness, getproperty(s.family, kind))

############################################
#=
backend_bounding_box(b::Unreal, shapes::Shapes) =
  UnrealBoundingBox(connection(b), collect_ref(shapes))
=#

unreal"public void SetView(Vector3 position, Vector3 target, float lens)"
unreal"public Vector3 ViewCamera()"
unreal"public Vector3 ViewTarget()"
unreal"public float ViewLens()"

set_view(camera::Loc, target::Loc, lens::Real, b::Unreal) =
  let c = connection(b)
    UnrealSetView(c, camera, target, lens)
    interrupt_processing(c)
  end

get_view(b::Unreal) =
  let c = connection(b)
    UnrealViewCamera(c), UnrealViewTarget(c), UnrealViewLens(c)
  end

zoom_extents(b::Unreal) = UnrealZoomExtents(connection(b))

view_top(b::Unreal) = UnrealViewTop(connection(b))

unreal"public void DeleteAll()"
unreal"public void DeleteMany(Actor[] objs)"
delete_all_shapes(b::Unreal) = UnrealDeleteAll(connection(b))

backend_delete_shapes(b::Unreal, shapes::Shapes) =
  UnrealDeleteMany(connection(b), collect_ref(shapes))

set_length_unit(unit::String, b::Unreal) = nothing # Unused, for now

#=
# Dimensions

const UnrealDimensionStyles = Dict(:architectural => "_ARCHTICK", :mechanical => "")

dimension(p0::Loc, p1::Loc, p::Loc, scale::Real, style::Symbol, b::Unreal=current_backend()) =
    UnrealCreateAlignedDimension(connection(b), p0, p1, p,
        scale,
        UnrealDimensionStyles[style])

dimension(p0::Loc, p1::Loc, sep::Real, scale::Real, style::Symbol, b::Unreal=current_backend()) =
    let v = p1 - p0
        angle = pol_phi(v)
        dimension(p0, p1, add_pol(p0, sep, angle + pi/2), scale, style, b)
    end

=#

# Layers
# Experiment for multiple, simultaneous, alternative layers
# Layers
unreal"public Actor CreateParent(String name)"
unreal"public Actor CurrentParent()"
unreal"public Actor SetCurrentParent(Actor newParent)"
unreal"public void SetActive(Actor obj, bool state)"
unreal"public void DeleteAllInParent(Actor parent)"
unreal"public void SwitchToParent(Actor newParent)"

UnrealLayer = Int

current_layer(b::Unreal)::UnrealLayer =
  UnrealCurrentParent(connection(b))

current_layer(layer::UnrealLayer, b::Unreal) =
  UnrealSetCurrentParent(connection(b), layer)

create_layer(name::String, b::Unreal) =
  UnrealCreateParent(connection(b), name)

set_layer_active(layer::UnrealLayer, status, b::Unreal) =
  let c = connection(b)
    UnrealSetActive(c, layer, status)
    interrupt_processing(c)
  end

delete_all_shapes_in_layer(layer::UnrealLayer, b::Unreal) =
  UnrealDeleteAllInParent(connection(b), layer)

switch_to_layer(layer::UnrealLayer, b::Unreal) =
  UnrealSwitchToParent(connection(b), layer)

# Experiment to speed up things

canonicalize_layer(layer::UnrealLayer, b::Unreal) =
  UnrealCanonicalize(connection(b), layer)

# Materials

UnrealMaterial = Int

unreal"public Material LoadMaterial(String name)"
unreal"public void SetCurrentMaterial(Material material)"
unreal"public Material CurrentMaterial()"

current_material(b::Unreal)::UnrealMaterial =
  UnrealCurrentMaterial(connection(b))

current_material(material::UnrealMaterial, b::Unreal) =
  UnrealSetCurrentMaterial(connection(b), material)

get_material(name::String, b::Unreal) =
  UnrealLoadMaterial(connection(b), name)


# Blocks

unreal"public Actor CreateBlockInstance(Actor block, Vector3 position, Vector3 vx, Vector3 vy, float scale)"
unreal"public Actor CreateBlockFromShapes(String name, Actor[] objs)"

realize(b::Unreal, s::Block) =
  s.shapes == [] ?
    UnrealLoadResource(connection(b), s.name) :
    UnrealCanonicalize(connection(b), UnrealCreateBlockFromShapes(connection(b), s.name, collect_ref(s.shapes)))

realize(b::Unreal, s::BlockInstance) =
    UnrealCreateBlockInstance(
        connection(b),
        ref(s.block).value,
        s.loc, vy(1, s.loc.cs), vz(1, s.loc.cs), s.scale)
#=

# Manual process
@time for i in 1:1000 for r in 1:10 circle(x(i*10), r) end end

# Create block...
Khepri.create_block("Foo", [circle(radius=r) for r in 1:10])

# ...and instantiate it
@time for i in 1:1000 Khepri.instantiate_block("Foo", x(i*10), 0) end

=#

# Lights
unreal"public Actor PointLight(Vector3 position, Color color, float range, float intensity)"

backend_pointlight(b::Unreal, loc::Loc, color::RGB, range::Real, intensity::Real) =
    UnrealPointLight(connection(b), loc, color, range, intensity)
#=
backend_spotlight(b::Unreal, loc::Loc, dir::Vec, hotspot::Real, falloff::Real) =
    UnrealSpotLight(connection(b), loc, hotspot, falloff, loc + dir)

backend_ieslight(b::Unreal, file::String, loc::Loc, dir::Vec, alpha::Real, beta::Real, gamma::Real) =
    UnrealIESLight(connection(b), file, loc, loc + dir, vxyz(alpha, beta, gamma))

# User Selection
=#

shape_from_ref(r, b::Unreal) =
  let idx = findfirst(s -> r in collect_ref(s), collected_shapes())
    if isnothing(idx)
      let c = connection(b)
          unknown(r, backend=b, ref=LazyRef(b, UnrealNativeRef(r), 0, 0))
          #code = UnrealShapeCode(c, r),
          #ref = LazyRef(b, UnrealNativeRef(r))
          #error("Unknown shape with code $(code)")
      end
    else
      collected_shapes()[idx]
    end
  end
#
#=
Unreal"public Point3d[] GetPosition(string prompt)"

select_position(prompt::String, b::Unreal) =
  begin
    @info "$(prompt) on the $(b) backend."
    let ans = UnrealGetPosition(connection(b), prompt)
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

Unreal"public ObjectId[] GetPoint(string prompt)"

# HACK: The next operations should receive a set of shapes to avoid re-creating already existing shapes

select_point(prompt::String, b::Unreal) =
  select_with_prompt(prompt, b, UnrealGetPoint)

Unreal"public ObjectId[] GetCurve(string prompt)"

select_curve(prompt::String, b::Unreal) =
  select_with_prompt(prompt, b, UnrealGetCurve)

Unreal"public ObjectId[] GetSurface(string prompt)"

select_surface(prompt::String, b::Unreal) =
  select_with_prompt(prompt, b, UnrealGetSurface)

Unreal"public ObjectId[] GetSolid(string prompt)"

select_solid(prompt::String, b::Unreal) =
  select_with_prompt(prompt, b, UnrealGetSolid)

Unreal"public ObjectId[] GetShape(string prompt)"

select_shape(prompt::String, b::Unreal) =
  select_with_prompt(prompt, b, UnrealGetShape)

Unreal"public long GetHandleFromShape(Entity e)"
Unreal"public ObjectId GetShapeFromHandle(long h)"

captured_shape(b::Unreal, handle) =
  shape_from_ref(UnrealGetShapeFromHandle(connection(b), handle),
                 b)

generate_captured_shape(s::Shape, b::Unreal) =
    println("captured_shape(autocad, $(UnrealGetHandleFromShape(connection(b), ref(s).value)))")

# Register for notification

Unreal"public void RegisterForChanges(ObjectId id)"
Unreal"public void UnregisterForChanges(ObjectId id)"
Unreal"public ObjectId[] ChangedShape()"
Unreal"public void DetectCancel()"
Unreal"public void UndetectCancel()"
Unreal"public bool WasCanceled()"

register_for_changes(s::Shape, b::Unreal) =
    let conn = connection(b)
        UnrealRegisterForChanges(conn, ref(s).value)
        UnrealDetectCancel(conn)
        s
    end

unregister_for_changes(s::Shape, b::Unreal) =
    let conn = connection(b)
        UnrealUnregisterForChanges(conn, ref(s).value)
        UnrealUndetectCancel(conn)
        s
    end

waiting_for_changes(s::Shape, b::Unreal) =
    ! UnrealWasCanceled(connection(b))

changed_shape(ss::Shapes, b::Unreal) =
    let conn = connection(b)
        changed = []
        while length(changed) == 0 && ! UnrealWasCanceled(conn)
            changed =  UnrealChangedShape(conn)
            sleep(0.1)
        end
        if length(changed) > 0
            shape_from_ref(changed[1], b)
        else
            nothing
        end
    end

Unreal"public ObjectId[] GetAllShapes()"
Unreal"public ObjectId[] GetAllShapesInLayer(ObjectId layerId)"

# HACK: This should be filtered on the plugin, not here.
all_shapes(b::Unreal) =
    let c = connection(b)
        Shape[shape_from_ref(r, b)
              for r in filter(r -> UnrealShapeCode(c, r) != 0, UnrealGetAllShapes(c))]
    end

all_shapes_in_layer(layer, b::Unreal) =
    let c = connection(b)
        Shape[shape_from_ref(r, b) for r in UnrealGetAllShapesInLayer(c, layer)]
    end

disable_update(b::Unreal) =
    UnrealDisableUpdate(connection(b))

enable_update(b::Unreal) =
    UnrealEnableUpdate(connection(b))
# Render

=#
unreal"public void SetResolution(int width, int height)"
unreal"public void ScreenShot(String path)"

#render exposure: [-3, +3] -> [-6, 21]
convert_render_exposure(b::Unreal, v::Real) = -4.05*v + 8.8
#render quality: [-1, +1] -> [+1, +50]
convert_render_quality(b::Unreal, v::Real) = round(Int, 25.5 + 24.5*v)

render_view(path::String, b::Unreal) =
    let c = connection(b)
      UnrealSetResolution(c, render_width(), render_height())
      interrupt_processing(c)
      UnrealScreenShot(c, path)
    end

unreal"public void SelectActors(Actor[] objs)"

highlight_shape(s::Shape, b::Unreal) =
    UnrealSelectActors(connection(b), collect_ref(s))

highlight_shapes(ss::Shapes, b::Unreal) =
    UnrealSelectActors(connection(b), collect_ref(ss))


unreal"public void StartSelectingActor()"
unreal"public int SelectedActorId(bool existing)"

select_shape(prompt::String, b::Unreal) =
  select_one_with_prompt(prompt, b, (c, prompt) ->
    let s = -2 # Means not found
      UnrealStartSelectingActor(c)
      while s == -2
        sleep(0.1)
        s = UnrealSelectedActorId(c, true)
      end
      [s]
    end)
=#
