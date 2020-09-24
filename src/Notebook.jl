export notebook

####################################################
abstract type PlotKey end
const PlotId = Any #PlotlyJS.GenericTrace
const PlotRef = GenericRef{PlotKey, PlotId}
const PlotNativeRef = NativeRef{PlotKey, PlotId}
const PlotUnionRef = UnionRef{PlotKey, PlotId}
const PlotSubtractionRef = SubtractionRef{PlotKey, PlotId}

const PlotMaterial = String

struct PlotBackend{K,T} <: Backend{K,T}
  connection::LazyParameter{PlotlyJS.SyncPlot}
end

connection(b::PlotBackend{K,T}) where {K,T} = b.connection()

const PLOT = PlotBackend{PlotKey, PlotId}

has_boolean_ops(::Type{PLOT}) = HasBooleanOps{false}()
void_ref(b::PLOT) = PlotNativeRef(-1)

create_plot_connection() =
  let (width, height) = render_size(),
      layout =
        Layout(autosize=true, width=width, height=height,
               margin=attr(l=0, r=0, b=0, t=0),
               showlegend=false,
               scene=attr(
                 #camera=attr(
                #   center=attr(x=target.x,
                #               y=target.y,
                #               z=target.z),
                #   eye=attr(x=camera.x-target.x,
                #            y=camera.y-target.y,
                #            z=camera.z-target.z)),
                 aspectmode="data",
                 xaxis=attr(showticklabels=false),
                 yaxis=attr(showticklabels=false),
                 zaxis=attr(showticklabels=false))),
      s = PlotlyJS.SyncPlot(Plot(AbstractTrace[], layout))
    display(s)
    s
  end

notebook = PLOT(LazyParameter(PlotlyJS.SyncPlot, create_plot_connection))

reset_backend(b::PLOT) =
  let param = b.connection
    if ! isnothing(param.value)
      param.value = (copy(param.value))
      display(param.value)
    end
  end

new_backend(b::PLOT) =
  reset(b.connection)

#

set_view(camera::Loc, target::Loc, lens::Real, aperture::Real, b::Plot) =
  begin
    b.camera = camera
    b.target = target
    b.lens = lens
  end

get_view(b::Plot) =
  b.camera, b.target, b.lens

###################################
#
delete_all_shapes(b::PLOT) =
  let c = connection(b),
      idxs = collect(1:length(c.plot.data))
    if length(idxs) > 0
      deletetraces!(c, idxs...)
    end
  end

backend_delete_shapes(b::PLOT, shapes::Shapes) =
  let c = connection(b),
      idxs = [findfirst(==(ref(s).value), c.plot.data) for s in shapes]
    if length(idxs) > 0
      deletetraces!(c, idxs...)
    end
  end

#=
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
    @remote(b, ClosedPolyLine(path.vertices))
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

backend_stroke_line(b::ACAD, vs) = @remote(b, PolyLine(vs))

backend_stroke_arc(b::ACAD, center::Loc, radius::Real, start_angle::Real, amplitude::Real) =
  let end_angle = start_angle + amplitude
    @remote(b, Arc(center, vz(1, center.cs), radius, start_angle, end_angle))
  end
backend_stroke_unite(b::ACAD, refs) = @remote(b, JoinCurves(refs))



realize(b::ACAD, s::EmptyShape) =
  ACADEmptyRef()
realize(b::ACAD, s::UniversalShape) =
  ACADUniversalRef()
realize(b::ACAD, s::Point) =
  @remote(b, Point(s.position))

=#

backend_line(b::PLOT, vs::Locs) =
  let pts = map(in_world, vs),
      r = PlotlyJS.scatter3d(
         x=map(cx, pts),
         y=map(cy, pts),
         z=map(cz, pts),
         line_shape="linear",
         marker_size=2,
         autocolorscale=false,
         showscale=false,
         hoverinfo="skip")
    PlotlyJS.addtraces!(connection(b), r)
    PlotNativeRef(r)
  end

realize(b::PLOT, s::Line) =
  backend_line(p, s.vertices)

realize(b::PLOT, s::Spline) = # This should be merged with opensplinepath
  if (s.v0 == false) && (s.v1 == false)
    let pts = map_division(in_world, spline_path(s.points), length(s.points)*4),
        r = PlotlyJS.scatter3d(
           x=map(cx, pts),
           y=map(cy, pts),
           z=map(cz, pts),
           line_shape="spline",
           marker_size=2,
           autocolorscale=false,
           showscale=false,
           hoverinfo="skip")
      PlotlyJS.addtraces!(connection(b), r)
      PlotNativeRef(r)
    end
  elseif (s.v0 != false) && (s.v1 != false)
    @remote(b, InterpSpline(s.points, s.v0, s.v1))
  else
    @remote(b, InterpSpline(
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1))
  end

#=
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
=#

realize(b::PLOT, s::Polygon) =
  backend_line(p, [s.vertices..., s.vertices[1]])

#=
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


realize(b::ACAD, s::SurfacePolygon) =
  @remote(b, SurfaceClosedPolyLine(s.vertices))
realize(b::ACAD, s::SurfaceRegularPolygon) =
  @remote(b, SurfaceClosedPolyLine(regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))
realize(b::ACAD, s::SurfaceRectangle) =
  @remote(b, SurfaceClosedPolyLine(
    [s.corner,
     add_x(s.corner, s.dx),
     add_xy(s.corner, s.dx, s.dy),
     add_y(s.corner, s.dy)]))
realize(b::ACAD, s::Surface) =
  let #ids = map(r->@remote(b, NurbSurfaceFrom(r)), @remote(b, SurfaceFromCurves(collect_ref(s.frontier))))
      ids = @remote(b, SurfaceFromCurves(collect_ref(s.frontier)))
    foreach(mark_deleted, s.frontier)
    ids
  end
backend_surface_boundary(b::ACAD, s::Shape2D) =
    map(c -> shape_from_ref(c, b), @remote(b, CurvesFromSurface(ref(s).value)))


realize(b::Plot, s::Sphere) =
  let mat = get_material(b, s)
    write_Plot_object(buffer(b), "sphere", mat, in_world(s.center), s.radius)
    void_ref(b)
  end

realize(b::Plot, s::Torus) =
  let buf = buffer(b)
    write_Plot_object(buf, "torus", get_material(b, s), s.re, s.ri) do
      let p = in_world(s.center),
          t = s.center.cs.transform
        write_Plot_object(buf, "matrix", nothing,
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
  write_Plot_object(buffer(b), "prism", get_material(b, s), 0, norm(s.v)) do
    s.bs,
                              map(p -> (p + s.v), s.bs))
## FIXME: deal with the rotation angle
realize(b::ACAD, s::RightCuboid) =
  ACADCenteredBox(connection(b), s.cb, s.width, s.height, s.h)
=#

realize(b::Plot, s::Box) =
  let buf = buffer(b),
      bot = in_world(s.c),
      top = in_world(s.c + vxyz(s.dx, s.dy, s.dz, s.c.cs)),
      mat = get_material(b, s)
    write_Plot_object(buf, "box", mat, bot, top)
    void_ref(b)
  end

realize(b::Plot, s::Cone) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_Plot_object(buf, "cone", mat, bot, s.r, top, 0)
    void_ref(b)
  end

realize(b::Plot, s::ConeFrustum) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_Plot_object(buf, "cone", mat, bot, s.rb, top, s.rt)
    void_ref(b)
  end

realize(b::Plot, s::Cylinder) =
  let buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      mat = get_material(b, s)
    write_Plot_object(buf, "cylinder", mat, bot, top, s.r)
    void_ref(b)
  end

write_Plot_mesh(buf::IO, mat, points, closed_u, closed_v, smooth_u, smooth_v) =
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
    write_Plot_object(buf, "mesh2", mat) do
      write_Plot_object(buf, "vertex_vectors", nothing, si*sj, reshape(permutedims(pts), :)...)
      # Must understand how to handle smoothness along one direction
      if smooth_u && smooth_v
        write_Plot_object(buf, "normal_vectors", nothing, si*sj, reshape(permutedims(vcs), :)...)
      end
      write_Plot_object(buf, "face_indices", nothing, length(idxs), idxs...)
      # Must understand how to handle smoothness along one direction
      #write_Plot_object(buf, "normal_indices", nothing, length(idxs), idxs...)
    end
  end
=#

realize(b::PLOT, s::SurfaceGrid) =
  let mat = 1, #get_material(b, s)
      pts = map(in_world, s.points),
      pts = [pts[i,:] for i in 1:size(pts,1)],
      r = PlotlyJS.surface(
         x=map(r->map(cx, r), pts),
         y=map(r->map(cy, r), pts),
         z=map(r->map(cz, r), pts),
         autocolorscale=false,
         showscale=false,
         hoverinfo="skip",
         colorscale=[[0, "rgb(40,40,40)"], [1, "rgb(210,210,210)"]])
    PlotlyJS.addtraces!(connection(b), r)
    PlotNativeRef(r)
  end

#=
realize(b::Plot, s::SweepPath) =
  let vertices = in_world.(path_vertices(s.profile)),
      frames = map_division(identity, s.path, 20),
      buf = buffer(b),
      mat = get_material(b, s)
    write_Plot_mesh(
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
realize(b::PLOT, s::Thicken) =
  realize(b, s.shape)

realize(b::PLOT, s::EmptyShape) = void_ref(b)

realize(b::PLOT, s::UniversalShape) = void_ref(b)

#=
#=
realize(b::Plot, s::Move) =
    let r = map_ref(s.shape) do r
                PlotMove(connection(b), r, s.v)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::Plot, s::Scale) =
    let r = map_ref(s.shape) do r
                PlotScale(connection(b), r, s.p, s.s)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::Plot, s::Rotate) =
    let r = map_ref(s.shape) do r
                PlotRotate(connection(b), r, s.p, s.v, s.angle)
                r
            end
        mark_deleted(s.shape)
        r
    end

=#

realize(b::Plot, s::UnionShape) =
  let shapes = filter(! is_empty_shape, s.shapes)
    length(shapes) == 1 ?
      (ref(shapes[1]); delete_shape(shapes[1])) :
      write_Plot_object(buffer(b), "union", get_material(b, s)) do
        for ss in shapes
          ref(ss)
          delete_shape(ss)
        end
      end
    void_ref(b)
  end


realize(b::Plot, s::IntersectionShape) =
  write_Plot_object(buffer(b), "intersection", get_material(b, s)) do
    for ss in s.shapes
      ref(ss)
      delete_shape(ss)
    end
    void_ref(b)
  end

realize(b::Plot, s::SubtractionShape3D) =
  write_Plot_object(buffer(b), "difference", get_material(b, s)) do
    ref(s.shape)
    delete_shape(s.shape)
    for ss in s.shapes
      ref(ss)
      delete_shape(ss)
    end
    void_ref(b)
  end

# BIM

realize_box(b::Plot, mat, p, dx, dy, dz) =
  let buf = buffer(b),
      bot = in_world(p),
      top = in_world(add_xyz(p, dx, dy, dz))
    write_Plot_object(buf, "box", mat, bot, top)
    void_ref(b)
  end

realize_prism(b::Plot, top, bot, side, path::PathSet, h::Real) =
  # PathSets require a different approach
  let buf = buffer(b),
      bot_vss = map(path_vertices, path.paths),
      top_vss = map(path_vertices, translate(path, vz(h)).paths)
    write_Plot_polygons(buf, bot, map(reverse, bot_vss))
    write_Plot_polygons(buf, top, top_vss)
    for (bot_vs, top_vs) in zip(bot_vss, top_vss)
      for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
        write_Plot_polygon(buf, side, vs)
      end
    end
    void_ref(b)
  end

realize_pyramid_frustum(b::Plot, top, bot, side, bot_vs::Locs, top_vs::Locs, closed=true) =
  let buf = buffer(b)
    if closed
      write_Plot_polygon(buf, bot, reverse(bot_vs))
      write_Plot_polygon(buf, top, top_vs)
    end
    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
      write_Plot_polygon(buf, side, vs)
    end
    void_ref(b)
  end

backend_surface_polygon(b::Plot, mat, path::PathSet, acw=true) =
  acw ?
    write_Plot_polygons(buffer(b), mat, map(path_vertices, path.paths)) :
    write_Plot_polygons(buffer(b), mat, map(reverse ∘ path_vertices, path.paths))

backend_surface_polygon(b::Plot, mat, vs::Locs, acw=true) =
  let buf = buffer(b)
    write_Plot_polygon(buf, mat, acw ? vs : reverse(vs))
  end

write_Plot_polygon(io::IO, mat, vs) =
  write_Plot_object(io, "polygon", mat, length(vs)+1, vs..., vs[1])

write_Plot_polygons(io::IO, mat, vss) =
  write_Plot_object(io, "polygon", mat,
    mapreduce(length, +, vss) + length(vss),
    mapreduce(vs->[vs..., vs[1]], vcat, vss)...)

# Polygons with holes need a PathSets in Plot

subtract_paths(b::Plot, c_r_w_path::PathSet, c_l_w_path::PathSet, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path.paths..., c_r_op_path),
  path_set(c_l_w_path.paths..., c_l_op_path)

subtract_paths(b::Plot, c_r_w_path::Path, c_l_w_path::Path, c_r_op_path, c_l_op_path) =
  path_set(c_r_w_path, c_r_op_path),
  path_set(c_l_w_path, c_l_op_path)


#=
Plot families need to know the different kinds of materials
that go on each surface.
In some cases it might be the same material, but in others, such
as slabs, outside walls, etc, we will have different materials.
=#

const PlotMaterialFamily = BackendMaterialFamily{PlotMaterial}
Plot_material_family(mat::PlotMaterial) =
  PlotMaterialFamily(mat)

const PlotSlabFamily = BackendSlabFamily{PlotMaterial}
Plot_slab_family(top::PlotMaterial, bot::PlotMaterial=top, side::PlotMaterial=bot) =
  PlotSlabFamily(top, bot, side)

const PlotRoofFamily = BackendRoofFamily{PlotMaterial}
Plot_roof_family(top::PlotMaterial, bot::PlotMaterial=top, side::PlotMaterial=bot) =
  PlotRoofFamily(top, bot, side)

const PlotWallFamily = BackendWallFamily{PlotMaterial}
Plot_wall_family(right::PlotMaterial, left::PlotMaterial=right) =
  PlotWallFamily(right, left)

export Plot_material_family,
       Plot_slab_family,
       Plot_roof_family,
       Plot_wall_family,
       default_Plot_material

# Layers
current_layer(b::Plot) =
  default_Plot_material()

current_layer(layer, b::Plot) =
  default_Plot_material(layer)

backend_create_layer(b::Plot, name::String, active::Bool, color::RGB) =
  begin
    @assert active
    Plot_material(name, red=red(color), green=green(color), blue=blue(color))
  end

#=
create_ground_plane(shapes, material=default_Plot_ground_material()) =
  if shapes == []
    error("No shapes selected for analysis. Use add-Plot-shape!.")
  else
    let (p0, p1) = bounding_box(union(shapes)),
        (center, ratio) = (quad_center(p0, p1, p2, p3),
                  distance(p0, p4)/distance(p0, p2));
     ratio == 0 ?
      error("Couldn"t compute height. Use add-Plot-shape!.") :
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
realize(b::Plot, s::Beam) =
    ref(right_cuboid(s.p0, 0.2, 0.2, s.p1, 0))
=#
realize(b::Plot, s::Union{Door, Window}) =
  void_ref(b)

used_materials(b::Plot) =
  unique(map(f -> realize(s.family, b), b.shapes))

####################################################


export export_to_jupyter
export_to_jupyter(b::Jupyter=current_backend()) =
  let camera = in_world(b.camera),
      target = in_world(b.target),
      (width, height) = render_size(),
      layout =
    PlotlyJS.Layout(autosize=true, width=width, height=height,
                    margin=attr(l=0, r=0, b=0, t=0),
                    camera_center_x=target.x,
                    camera_center_y=target.y,
                    camera_center_z=target.z,
                    camera_eye_x=camera.x,
                    camera_eye_y=camera.y,
                    camera_eye_z=camera.z,
                    scene_aspectmode="data",
                    scene_xaxis_showticklabels=false,
                    scene_yaxis_showticklabels=false,
                    scene_zaxis_showticklabels=false,
                    )
    PlotlyJS.plot([Khepri.ref(s).value for s in b.shapes if !is_empty_shape(s)], layout)
  end


=#
