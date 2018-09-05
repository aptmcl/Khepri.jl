export radiance,
       save_rad


write_rad_primitive(io::IO, modifier, typ, identifier, strings, ints, reals) =
    begin
        print_elems(elems) =
            begin
                print(io, length(elems))
                for e in elems print(io, " ", e) end
                println(io)
            end
        println(io, modifier, " ", typ, " ", identifier)
        print_elems(strings)
        print_elems(ints)
        print_elems(reals)
    end

write_rad_polygon(io::IO, modifier, id, vertices) =
    begin
        println(io, modifier, " ", "polygon", " ", id)
        println(io, 0) #0 strings
        println(io, 0) #0 ints
        println(io, 3*length(vertices))
        for v in vertices println(io, " ", v.x, " ", v.y, " ", v.z) end
    end

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

mutable struct IOBuffer_Backend{K,T} <: Backend{K,T}
  buffer::LazyParameter{IOBuffer}
  count::Integer
  materials::Dict
end

const Radiance = IOBuffer_Backend{RadianceKey, RadianceId}

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
    get!(b.materials, s.family, length(b.materials))

save_rad(path::String) =
    open(path, "w") do out
        write(out, String(take!(radiance.buffer())))
    end

#

current_backend(radiance)

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
    with(current_backend, autocad) do
        println(bot_vs)
        println(top_vs)
        irregular_pyramid_fustrum(bot_vs, top_vs)
    let bot_id = next_id(b, s)
        top_id = next_id(b, s)
        mod = next_modifier(b, s)
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
    let base = vz(s.level.height + s.family.coating_thickness - s.family.thickness)
        thickness = vz(s.family.coating_thickness + s.family.thickness)
        bot_vs = path_vertices(translate(s.contour, base))
        top_vs = path_vertices(translate(s.contour, base + thickness))
        realize_pyramid_fustrum(b, s, "slab", bot_vs, top_vs)
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

realize(b::Radiance, s::Wall) =
    let base_height = s.bottom_level.height,
        height = s.top_level.height - base_height
        realize_wall_path(b, s, s.path, base_height, height, s.family.thickness)
    end

realize_wall_path(b::Radiance, s::Wall, path::OpenPolygonalPath, base_height::Real, height::Real, thickness::Real) =
    let vs = path_vertices(path)
        half_thickness = thickness/2
        bot_vs = [offset(vs, half_thickness	)...,
                  reverse!(offset(vs, -half_thickness))...]
        top_vs = map(p -> add_z(p, height), bot_vs)
        realize_pyramid_fustrum(b, s, "wall", bot_vs, top_vs)
    end

realize_wall_path(b::Radiance, s::Wall, path::Path, base_height::Real, height::Real, thickness::Real) =
    realize_wall_path(b, s, convert(OpenPolygonalPath, path), base_height, height, thickness)
