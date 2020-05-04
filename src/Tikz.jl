export tikz,
       tikz_output,
       display_tikz_output

tikz_e(out::IO, arg) =
  begin
    print(out, arg)
    println(out, ";")
  end

tikz_draw(out::IO, filled=false) = print(out, filled ? "\\fill " : "\\draw ")

tikz_number(out::IO, x::Real) =
  isinteger(x) ? print(out, x) : (abs(x) < 0.0001 ? print(out, 0) : print(out, round(x*10000.0)/10000.0))

tikz_cm(out::IO, x::Real) = begin
  tikz_number(out, x)
  print(out, "cm")
end

tikz_coord(out::IO, c::Loc) =
  begin
    print(out, "(")
    tikz_number(out, c.x)
    print(out, ",")
    tikz_number(out, c.y)
    if !iszero(c.z)
      print(out, ",")
      tikz_number(out, c.z)
    end
    print(out, ")")
  end

tikz_pgfpoint(out::IO, c::Loc) =
  begin
    if ! iszero(c.z)
      error("Can't handle 3D coords")
    end
    print(out, "\\pgfpoint{")
    tikz_cm(out, c.x)
    print(out, "}{")
    tikz_cm(out, c.y)
    print(out, "}")
  end

tikz_circle(out::IO, c::Loc, r::Real, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    tikz_coord(out, c)
    print(out, "circle(")
    tikz_cm(out, r)
    tikz_e(out, ")")
  end


tikz_point(out::IO, c::Loc) = tizk_circle(out, c, 0.01, true)

tikz_ellipse(out::IO, c::Loc, r0::Real, r1::Real, fi::Real, filled=false) =
  begin
    tikz_draw(out, filled)
    print(out, "[shift={")
    tikz_coord(out, c)
    print(out, "}]")
    print(out, "[rotate=")
    tikz_number(out, rad2deg(fi))
    print(out, "]")
    print(out, "(0,0)")
    print(out, "ellipse(")
    tikz_cm(out, r0)
    print(out, " and ")
    tikz_cm(out, r1)
    tikz_e(out, ")")
  end

tikz_arc(out::IO, c::Loc, r::Real, ai::Real, af::Real, filled=false) =
  begin
    tikz_draw(out, filled)
    if filled
      tikz_coord(out, c)
      print(out, "--")
    end
    tikz_coord(out, c+vpol(r, ai))
    print(out, "arc(")
    tikz_number(out, rad2deg(ai))
    print(out, ":")
    tikz_number(out, rad2deg(ai > af ? af+2*pi : af))
    print(out, ":")
    tikz_cm(out, r)
    print(out, ")")
    if filled
      tikz_e(out, "--cycle")    end
    println(out, ";")
  end

tikz_maybe_arc(out::IO, c::Loc, r::Real, ai::Real, da::Real, filled=false) =
  if iszero(r)
    tikz_point(out, c)
  elseif iszero(da)
    tikz_point(out, c + vpol(r, ai))
  elseif abs(da) >= 2*pi
    tikz_circle(out, c, r, filled)
  else
    let af = ai + da
      if af > ai
        tikz_arc(out, c, r, ai, af, filled)
      else
        tikz_arc(out, c, r, af, ai, filled)
      end
    end
  end

tikz_line(out::IO, pts::Locs) =
  begin
    tikz_draw(out, false)
    tikz_coord(out, first(pts))
    for pt in Iterators.drop(pts, 1)
      print(out, "--")
      tikz_coord(out, pt)
    end
    println(out, ";")
  end

tikz_dimension(out::IO, p::Loc, q::Loc, text::AbstractString) =
  begin
    print(out, "\\dimline{")
    tikz_coord(out, p)
    print(out, "}{")
    tikz_coord(out, q)
    print(out, "}{;")
    print(out, text)
    println(out, "};")
  end

tikz_closed_line(out::IO, pts::Locs, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    for pt in pts
      tikz_coord(out, pt)
      print(out, "--")
    end
    tikz_e(out, "cycle")
  end

tikz_spline(out::IO, pts::Locs, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    print(out, "plot [smooth,tension=1] coordinates {")
    for pt in pts
      tikz_coord(out, pt)
    end
    tikz_e(out, "}")
  end

tikz_closed_spline(out::IO, pts::Locs, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    print(out, "plot [smooth cycle,tension=1] coordinates {")
    for pt in pts
      tikz_coord(out, pt)
    end
    tikz_e(out, "}")
  end

# HACK we need to handle the starting and ending vectors
tikz_hobby_spline(out::IO, pts::Locs, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    print(out, "[hobby]")
    print(out, "plot coordinates {")
    for pt in pts
      tikz_coord(out, pt)
    end
    tikz_e(out, "}")
  end

# HACK we need to handle the starting and ending vectors
tikz_hobby_closed_spline(out::IO, pts::Locs, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    print(out, "[closed hobby]")
    print(out, "plot coordinates {")
    for pt in pts
      tikz_coord(out, pt)
    end
    tikz_e(out, "}")
  end

tikz_rectangle(out::IO, p::Loc, w::Real, h::Real, filled::Bool=false) =
  begin
    tikz_draw(out, filled)
    tikz_coord(out, p)
    print(out, "rectangle")
    tikz_coord(out, p+vxy(w, h))
    println(out, ";")
  end

# Assuming default Arial font for AutoCAD
tikz_text(out::IO, txt, p::Loc, h::Real) =
  let (scale_x, scale_y) = (3.7*h, 3.7*h)
    tikz_draw(out)
    print(out, "[anchor=base west]")
    tikz_coord(out, p)
    print(out, "node[font=\\fontfamily{phv}\\selectfont,outer sep=0pt,inner sep=0pt")
    print(out, ",xscale=")
    tikz_number(out, scale_x)
    print(out, ",yscale=")
    tikz_number(out, scale_y)
    print(out, "]{")
    print(out, txt)
    tikz_e(out, "}")
  end

# FINISH THIS
tikz_transform(out::IO, f::Function, c::Loc) =
  begin
    print(out, "\\begin{scope}")
    print(out, "[shift={")
    tikz_coord(out, c)
    print(out, "}]")
    print(out, "]")
    f(out)
    tikz_e(out, "\\end{scope}")
  end

tikz_set_view(out::IO, camera::Loc, target::Loc, lens::Real) =
  let (v, contents) = (camera-target, get_accumulated_tikz(out))
    print("\\tdplotsetmaincoords{")
    tikz_number(radians_>degrees(sph_psi(v)))
    print("}{")
    tikz_number(radians_>degrees(sph_phi(v))+90)
    println("}")
    println("\\begin{tikzpicture}[tdplot_main_coords]")
    write_bytes(contents, tikz_port)
    println("\\end{tikzpicture}")
  end


#

abstract type TikZKey end
const TikZId = Nothing
const TikZIds = Vector{TikZId}
const TikZRef = GenericRef{TikZKey, TikZId}
const TikZRefs = Vector{TikZRef}
const TikZEmptyRef = EmptyRef{TikZKey, TikZId}
const TikZUniversalRef = UniversalRef{TikZKey, TikZId}
const TikZNativeRef = NativeRef{TikZKey, TikZId}
const TikZUnionRef = UnionRef{TikZKey, TikZId}
const TikZSubtractionRef = SubtractionRef{TikZKey, TikZId}
const TikZ = IOBufferBackend{TikZKey, TikZId}

void_ref(b::TikZ) = TikZNativeRef(nothing)

create_TikZ_connection() = IOBuffer()

const tikz = TikZ(create_TikZ_connection())

get_accumulated_tikz(out::IO) = String(take!(out))

tikz_output(b::TikZ=current_backend()) = get_accumulated_tikz(connection(b))

display_tikz_output(b::TikZ=current_backend()) = print(tikz_output(b))

withTikZXForm(f, out, c) =
  if is_world_cs(c.cs)
    f(out, c)
  else
    tikz_transform(out,
      out -> f(out, u0(world_cs)),
      c)
   end

delete_all_shapes(b::TikZ) =
  truncate(connection(b), 0)

realize(b::TikZ, s::Point) =
  tikz_pgfpoint(connection(b), in_world(s.position))

realize(b::TikZ, s::Line) =
  tikz_line(connection(b), map(in_world, s.vertices))

realize(b::TikZ, s::Spline) =
  if (s.v0 == false) && (s.v1 == false)
    #TikZSpline(connection(b), s.points)
    tikz_hobby_spline(connection(b), map(in_world, s.points), false)
  elseif (s.v0 != false) && (s.v1 != false)
    TikZInterpSpline(connection(b), s.points, s.v0, s.v1)
  else
    TikZInterpSpline(connection(b),
                     s.points,
                     s.v0 == false ? s.points[2]-s.points[1] : s.v0,
                     s.v1 == false ? s.points[end-1]-s.points[end] : s.v1)
  end

realize(b::TikZ, s::ClosedSpline) =
  tikz_hobby_closed_spline(connection(b), map(in_world, s.points))

realize(b::TikZ, s::Circle) =
  withTikZXForm(connection(b), s.center) do out, c
    tikz_circle(out, c, s.radius)
  end
realize(b::TikZ, s::SurfaceCircle) =
  withTikZXForm(connection(b), s.center) do out, c
    tikz_circle(out, c, s.radius, true)
  end

realize(b::TikZ, s::Arc) =
  withTikZXForm(connection(b), s.center) do out, c
    tikz_maybe_arc(out, c, s.radius, s.start_angle, s.amplitude, false)
  end

realize(b::TikZ, s::SurfaceArc) =
  withTikZXForm(connection(b), s.center) do out, c
    tikz_maybe_arc(out, c, s.radius, s.start_angle, s.amplitude, true)
  end

realize(b::TikZ, s::Ellipse) =
  withTikZXForm(connection(b), s.center) do out, c
    tikz_ellipse(out, c, s.radius_x, s.radius_y, 0, false)
  end

realize(b::TikZ, s::SurfaceEllipse) =
  withTikZXForm(connection(b), s.center) do out, c
    tikz_ellipse(out, c, s.radius_x, s.radius_y, 0, true)
  end

realize(b::TikZ, s::EllipticArc) =
  error("Finish this")

#realize(b::TikZ, s::SurfaceElliptic_Arc) = TikZCircle(connection(b),

realize(b::TikZ, s::Polygon) =
  tikz_closed_line(connection(b), map(in_world, s.vertices))

realize(b::TikZ, s::SurfacePolygon) =
  tikz_closed_line(connection(b), map(in_world, s.vertices), true)

realize(b::TikZ, s::RegularPolygon) =
  tikz_closed_line(
    connection(b),
    map(in_world,
        regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)))

realize(b::TikZ, s::SurfaceRegularPolygon) =
  tikz_closed_line(
    connection(b),
    map(in_world,
        regular_polygon_vertices(s.edges, s.center, s.radius, s.angle, s.inscribed)),
    false)

realize(b::TikZ, s::Rectangle) =
  withTikZXForm(connection(b), s.corner) do out, c
    tikz_rectangle(out, c, s.dx, s.dy)
  end

realize(b::TikZ, s::SurfaceRectangle) =
  withTikZXForm(connection(b), s.corner) do out, c
    tikz_rectangle(out, c, s.dx, s.dy, true)
  end

realize(b::TikZ, s::Text) =
  withTikZXForm(connection(b), s.corner) do out, c
    tikz_text(out, s.str, c, s.height)
  end

#
realize(b::TikZ, s::SurfaceGrid) =
  let n = size(s.points,1),
      m = size(s.points,2)
    for i in 1:n
      tikz_hobby_spline(connection(b), map(in_world, s.points[i,:]), false)
    end
    for j in 1:m
      tikz_hobby_spline(connection(b), map(in_world, s.points[:,j]), false)
    end
  end

###
# Dimensioning

#backend_dimension(b::TikZ, pa, pb, sep, scale, style) =
#  tikz_dimension(connection(b), pa, pb, )
