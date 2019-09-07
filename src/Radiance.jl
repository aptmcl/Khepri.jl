export radiance,
       save_rad

write_rad_primitive(io::IO, modifier, typ, identifier, strings, ints, reals) =
  begin
    print_elems(elems) =
      begin
        print(io, length(elems))
        for e in elems print(io, " ", e) end
        println(io)
        # Some viewers are (incorrectly) restrictive regarding the RAD format
        #for i in 1:length(elems)
        #  print(io, " ", elems[i])
        #  if i%3 == 0
        #    println(io)
        #  end
        #end
      end
    println(io, modifier, " ", typ, " ", identifier)
    print_elems(strings)
    print_elems(ints)
    print_elems(reals)
    println(io)
  end

write_rad_polygon(io::IO, modifier, id, vertices) =
  begin
    with(current_backend, autocad, immediate_mode, true) do
      polygon(vertices)
    end
    println(io, modifier, " ", "polygon", " ", id)
    println(io, 0) #0 strings
    println(io, 0) #0 ints
    println(io, 3*length(vertices))
    for v in vertices println(io, " ", v.x, " ", v.y, " ", v.z) end
  end

write_rad_cone(io::IO, modifier, id, bot, bot_radius, top, top_radius) =
  write_rad_primitive(io, modifier, "cone", id, [], [],
    [bot.x, bot.y, bot.z,
     top.x, top.y, top.z,
     bot_radius, top_radius])

write_rad_cup(io::IO, modifier, id, bot, bot_radius, top, top_radius) =
  write_rad_primitive(io, modifier, "cup", id, [], [],
    [bot.x, bot.y, bot.z,
     top.x, top.y, top.z,
     bot_radius, top_radius])

write_rad_cylinder(io::IO, modifier, id, bot, radius, top) =
  write_rad_primitive(io, modifier, "cylinder", id, [], [],
    [bot.x, bot.y, bot.z,
     top.x, top.y, top.z,
     radius])

write_rad_tube(io::IO, modifier, id, bot, radius, top) =
  write_rad_primitive(io, modifier, "tube", id, [], [],
    [bot.x, bot.y, bot.z,
     top.x, top.y, top.z,
     radius])

write_rad_ring(io::IO, modifier, id, center, in_radius, out_radius, normal) =
   write_rad_primitive(io, modifier, "ring", id, [], [],
     [center.x, center.y, center.z,
      normal.x, normal.y, normal.z,
      in_radius, out_radius])

write_rad_quad(io::IO, modifier, id, sub_id, v0, v1, v2, v3) =
  begin
    println(io, modifier, " ", "polygon", " ", id, sub_id)
    println(io, 0) #0 strings
    println(io, 0) #0 ints
    println(io, 12)
    println(io, " ", v0.x, " ", v0.y, " ", v0.z)
    println(io, " ", v1.x, " ", v1.y, " ", v1.z)
    println(io, " ", v2.x, " ", v2.y, " ", v2.z)
    println(io, " ", v3.x, " ", v3.y, " ", v3.z)
  end

#=
Higher level primitives
=#

write_rad_box(io::IO, modifier, id, p0, l, w, h) =
  let p1 = p0 + vx(l),
      p2 = p0 + vxy(l, w),
      p3 = p0 + vy(w),
      p4 = p0 + vz(h),
      p5 = p4 + vx(l),
      p6 = p4 + vxy(l, w),
      p7 = p4 + vy(w)
    write_rad_quad(io, modifier, id, "face0", p0, p1, p5, p4)
    write_rad_quad(io, modifier, id, "face1", p1, p2, p6, p5)
    write_rad_quad(io, modifier, id, "face2", p2, p3, p7, p6)
    write_rad_quad(io, modifier, id, "face3", p3, p0, p4, p7)
    write_rad_quad(io, modifier, id, "face4", p3, p2, p1, p0)
    write_rad_quad(io, modifier, id, "face5", p2, p5, p6, p7)
  end

#

#=

We need to discretize paths so that we can extract the vertices
We will use some sort of tolerance to deal with curved paths

=#

abstract type RadianceKey end
const RadianceId = Int
const RadianceRef = GenericRef{RadianceKey, RadianceId}
const RadianceNativeRef = NativeRef{RadianceKey, RadianceId}
const RadianceUnionRef = UnionRef{RadianceKey, RadianceId}
const RadianceSubtractionRef = SubtractionRef{RadianceKey, RadianceId}

mutable struct RadianceBackend{K,T} <: Backend{K,T}
  buffer::LazyParameter{IOBuffer}
  count::Integer
  materials::Dict
end

const Radiance = RadianceBackend{RadianceKey, RadianceId}

void_ref(b::Radiance) = RadianceNativeRef(-1)

create_radiance_buffer() = IOBuffer()

const radiance = Radiance(LazyParameter(IOBuffer, create_radiance_buffer), 0, Dict())

buffer(b::Radiance) = b.buffer()
next_id(b::Radiance, s::Shape) =
    begin
        b.count += 1
        b.count -1
    end
next_modifier(b::Radiance, s::Shape) =
    get!(b.materials, isdefined(s, :family) ? s.family : missing, length(b.materials))

save_rad(path::String) =
    open(path, "w") do out
        write(out, String(take!(radiance.buffer())))
    end

#

current_backend(radiance)

delete_all_shapes(b::Radiance) =
  take!(radiance.buffer())




#=
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
  ACADIrregularPyramidFrustum(connection(b),
                              s.bs,
                              map(p -> (p + s.v), s.bs))
## FIXME: deal with the rotation angle
realize(b::ACAD, s::RightCuboid) =
  ACADCenteredBox(connection(b), s.cb, s.width, s.height, s.h)
realize(b::ACAD, s::Box) =
  ACADBox(connection(b), s.c, s.dx, s.dy, s.dz)
realize(b::ACAD, s::Cone) =
  ACADCone(connection(b), add_z(s.cb, s.h), s.r, s.cb)
realize(b::ACAD, s::ConeFrustum) =
  ACADConeFrustum(connection(b), s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt)
=#
realize(b::Radiance, s::Cylinder) =
  let bot_id = next_id(b, s),
      top_id = next_id(b, s),
      side_id = next_id(b, s),
      mod = next_modifier(b, s),
      kind = "cylinder",
      buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      normal = unitized(top-bot)
    write_rad_ring(buf, "mat$(kind)bot$(mod)", bot_id, bot, 0, s.r, -normal)
    write_rad_cylinder(buf, "mat$(kind)side$(mod)", side_id, bot, s.r, top)
    write_rad_ring(buf, "mat$(kind)top$(mod)", top_id, top, 0, s.r, normal)
    bot_id
  end

#=
realize(b::Radiance, s::EmptyShape) =
    EmptyRef{RadianceId}()
realize(b::Radiance, s::UniversalShape) =
    UniversalRef{RadianceId}()

realize(b::Radiance, s::Move) =
    let r = map_ref(s.shape) do r
                RadianceMove(connection(b), r, s.v)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::Radiance, s::Scale) =
    let r = map_ref(s.shape) do r
                RadianceScale(connection(b), r, s.p, s.s)
                r
            end
        mark_deleted(s.shape)
        r
    end

realize(b::Radiance, s::Rotate) =
    let r = map_ref(s.shape) do r
                RadianceRotate(connection(b), r, s.p, s.v, s.angle)
                r
            end
        mark_deleted(s.shape)
        r
    end

=#

# BIM

realize_pyramid_fustrum(b::Radiance, s::Shape, kind::String, bot_vs, top_vs) =
  let bot_id = next_id(b, s),
      top_id = next_id(b, s),
      mod = next_modifier(b, s),
      buf = buffer(b)
    write_rad_polygon(buf, "mat$(kind)top$(mod)", bot_id, reverse(bot_vs))
    write_rad_polygon(buf, "mat$(kind)bot$(mod)", top_id, top_vs)
    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
      write_rad_polygon(buffer(b), "mat$(kind)side$(mod)", next_id(b, s), vs)
    end
    bot_id
  end
end

realize(b::Radiance, s::Slab) =
  let base = vz(s.level.height + s.family.coating_thickness - s.family.thickness),
      thickness = vz(s.family.coating_thickness + s.family.thickness),
      path = s.contour
    for op in s.openings
      path = subtract_paths(path, op)
    end
    realize_pyramid_fustrum(
      b, s, "slab",
      path_vertices(translate(path, base)),
      path_vertices(translate(path, base + thickness)))
  end

#=
#FIXME define the family parameters for beams
realize(b::Radiance, s::Beam) =
    ref(right_cuboid(s.p0, 0.2, 0.2, s.p1, 0))

=#
realize(b::Radiance, s::Panel) =
  let p1 = s.vertices[1],
      p2 = s.vertices[2],
      p3 = s.vertices[3],
      n = vz(s.family.thickness/2, cs_from_o_vx_vy(p1, p2-p1, p3-p1))
    realize_pyramid_fustrum(
        b, s, "panel",
        map(p -> in_world(p - n), s.vertices),
        map(p -> in_world(p + n), s.vertices))
  end

closed_path_for_height(path, h) =
  let ps = path_vertices(path)
    closed_polygonal_path([ps..., reverse(map(p -> p+vz(h), ps))...])
  end

#=
One important restriction is that Radiance only supports _planar_ polygons.
This forces us to create multiple subpaths.
=#
realize(b::Radiance, w::Wall) =
  let w_base_height = w.bottom_level.height,
      w_height = w.top_level.height - w_base_height,
      w_path = translate(w.path, vz(w_base_height)),
      w_thickness = w.family.thickness,
      r_thickness = (1+w.offset)/2*w_thickness,
      l_thickness = (1-w.offset)/2*w_thickness,
      r_w_path = closed_path_for_height(offset(w_path, r_thickness), w_height),
      l_w_path = closed_path_for_height(offset(w_path, l_thickness), w_height),
      openings = [w.doors..., w.windows...]
    @assert length(path_vertices(w_path)) == 2 # fix this for bigger numbers
    for op in openings
      let op_base_height = op.loc.y,
          op_height = op.family.height,
          op_thickness = op.family.thickness,
          op_path = translate(subpath(w_path, op.loc.x, op.loc.x + op.family.width), vz(op_base_height)),
          r_op_path = closed_path_for_height(offset(op_path, +half_thickness), op_height),
          l_op_path = closed_path_for_height(offset(op_path, -half_thickness), op_height)
        r_w_path = subtract_paths(r_w_path, r_op_path)
        l_w_path = subtract_paths(l_w_path, l_op_path)
      end
    end
    realize_pyramid_fustrum(b, w, "wall", path_vertices(r_w_path), path_vertices(l_w_path))
end

#=

Upon daysim processing, we need to compute the useful daylight illumination, which is the fraction of time that a sensor is within a given range of illumination

To do that, we start by reading the .ill file generated

ill = CSV.read(path, delim=' ', datarow=1)

The file looks like this:

8760×136 DataFrame. Omitted printing of 109 columns
│ Row  │ Column1 │ Column2 │ Column3  │ Column4 │ Column5 │ Column6 │ Column7 │ Column8 │ Column9 │ Column10 │ Column11 │ Column12 │ Column13 │ Column14 │ Column15 │ Column16 │ Column17 │ Column18 │ Column19 │ Column20 │ Column21 │ Column22 │ Column23 │ Column24 │ Column25 │ Column26 │ Column27 │
│      │ Int64⍰  │ Int64⍰  │ Float64⍰ │ Missing │ Int64⍰  │ Int64⍰  │ Int64⍰  │ Int64⍰  │ Int64⍰  │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │ Int64⍰   │
├──────┼─────────┼─────────┼──────────┼─────────┼─────────┼─────────┼─────────┼─────────┼─────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤
│ 1    │ 1       │ 1       │ 0.5      │ missing │ 0       │ 0       │ 0       │ 0       │ 0       │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │
│ 2    │ 1       │ 1       │ 1.5      │ missing │ 0       │ 0       │ 0       │ 0       │ 0       │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │
│ 3    │ 1       │ 1       │ 2.5      │ missing │ 0       │ 0       │ 0       │ 0       │ 0       │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │
│ 4    │ 1       │ 1       │ 3.5      │ missing │ 0       │ 0       │ 0       │ 0       │ 0       │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │
│ 5    │ 1       │ 1       │ 4.5      │ missing │ 0       │ 0       │ 0       │ 0       │ 0       │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │
│ 6    │ 1       │ 1       │ 5.5      │ missing │ 0       │ 0       │ 0       │ 0       │ 0       │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │ 0        │
⋮

Then, for each sensor (i.e., starting at column 5) we compute the number of times that its illuminance was within predefined limits and we normalize for all possible times

=#

#=
#path = "C:\\Users\\aml\\Dropbox\\AML\\Projects\\rosetta\\tests\\xyz.ill"

path = "C:\\Users\\aml\\Downloads\\Geometry_0_0_0629.ill"

ill = CSV.read(path, delim=' ', datarow=1)

udi_in_range(df, min, max) =
  let occupied_hours_per_year = nrow(df)
    colwise(col->
      round(Int, count(i -> min < i < max, col)/occupied_hours_per_year*100),
      df[2:end,5:end])
  end

res = udi_in_range(ill, 299.9, 3000.1)
maximum(res)
minimum(res)

test = DataFrame([
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 0 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
0 0 0 0 0 20 10;
])

udi_in_range(test, -1, 5)
=#


#
const CIE_Overcast_Sky_rad_contents = """
!gensky 12 21 12.00 -c -a 42.300 -o 71.100 -m 75.000 -B 100

skyfunc glow sky_mat
0
0
4 1 1 1 0

sky_mat source sky
0
0
4 0 0 1 180

skyfunc glow ground_glow
0
0
4 1 1 1  0

ground_glow source ground
0
0
4 0 0 -1 180
"""

path_replace_suffix(path::String, suffix::String) =
  let (base, old_suffix) = splitext(path)
    base * suffix
  end

export_CIE_Overcast_Sky(path::String) =
  let skypath = path_replace_suffix(path, "_sky.rad")
    open(skypath, "w") do out
      write(out, CIE_Overcast_Sky_rad_contents)
    end
    skypath
  end
