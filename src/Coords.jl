import LinearAlgebra.cross, LinearAlgebra.dot, LinearAlgebra.norm

export Loc, Locs, LocOrZ,
       Vec, Vecs, VecOrZ,
       XYZ, xyz, cyl, sph,
       VXYZ, vxyz, vcyl, vsph,
       world_cs,
       current_cs,
       distance,
       u0,ux,uy,uz,uxy,uyz,uxz,uxyz,
       x,y,z,
       xy,xz,yz,pol,cyl,sph,
       cx,cy,cz,
       pol_rho, pol_phi,
       cyl_rho,cyl_phi,cyl_z,
       sph_rho,sph_phi,sph_psi,
       uvx,uvy,uvz,uvxy,uvyz,uvxz,uvxyz,
       vx,vy,vz,
       vxy,vxz,vyz,vpol,vcyl,vsph,
       add_x,add_y,add_z,add_xy,add_xz,add_yz,add_xyz,
       add_pol,add_cyl,add_sph,
       unitized, dot, cross,
       cs_from_o_vx_vy_vz,
       cs_from_o_vx_vy,
       cs_from_o_vz,
       cs_from_o_phi,
       loc_from_o_vx_vy,
       loc_from_o_vz,
       loc_from_o_phi,
       min_loc, max_loc,
       is_world_cs,
       in_cs, in_world,
       intermediate_loc,
       meta_program,
       translated_cs,
       scaled_cs,
       center_scaled_cs,
       translating_current_cs,
       regular_polygon_vertices,
       norm,
       angle_between,
       rotate_vector

#=
Some useful terminology:

Coordinate Space (also known as Reference Frame): an origin (position) and
three axis (directions) that allow the precise identification of locations and
translations (which includes directions and magnitude). There can be many reference
frames. One is considered the world reference frame. It is possible to interpret
the same location or the same translation regarding different reference frames.

Coordinate System: a way to assign meaning to numbers that represent locations
or translations relative to a reference frame. Different coordinate systems (e.g.,
rectangular, polar, cylindrical, and spherical) assign different meanings to the
three numbers. It is possible to convert between different coordinate systems.

Location: represented by a triple of numbers using a Coordinate System and a
Reference Frame.

Translation: represented by a triple of numbers using a Coordinate System and a
Reference Frame.

We need to frequently add translations to locations but it might happen that these
translations and locations have different reference frames, which cause surprising
results. To avoid this problem, we introduce the concept of void reference frame,
a reference frame which does not

=#

Vec4f = SVector{4,Float64}
Mat4f = SMatrix{4,4,Float64}

struct CS
  transform::Mat4f
end

translated_cs(cs::CS, x::Real, y::Real, z::Real) =
    CS(cs.transform * SMatrix{4,4,Float64}([
        1 0 0 x;
        0 1 0 y;
        0 0 1 z;
        0 0 0 1]))

scaled_cs(cs::CS, x::Real, y::Real, z::Real) =
    CS(cs.transform * SMatrix{4,4,Float64}([
        x 0 0 0;
        0 y 0 0;
        0 0 z 0;
        0 0 0 1]))

center_scaled_cs(cs::CS, x::Real, y::Real, z::Real) =
    let xt = cs.transform[4,1]
        yt = cs.transform[4,2]
        zt = cs.transform[4,3]
        translated_cs(
            scaled_cs(
                translated_cs(cs, -xt, -yt, -zt),
                x, y, z),
            xt, yt, zt)
    end

export rotated_around_p_v_cs
rotated_around_p_v_cs(cs::CS, a::Real, b::Real, c::Real, u::Real, v::Real, w::Real, phi::Real) =
  let u2 = u*u,
      v2 = v*v,
      w2 = w*w,
      cosT = cos(phi),
      oneMinusCosT = 1-cosT,
      sinT = sin(phi),
      m11 = u2 + (v2 + w2) * cosT,
      m12 = u*v * oneMinusCosT - w*sinT,
      m13 = u*w * oneMinusCosT + v*sinT,
      m14 = (a*(v2 + w2) - u*(b*v + c*w))*oneMinusCosT + (b*w - c*v)*sinT,
      m21 = u*v * oneMinusCosT + w*sinT,
      m22 = v2 + (u2 + w2) * cosT,
      m23 = v*w * oneMinusCosT - u*sinT,
      m24 = (b*(u2 + w2) - v*(a*u + c*w))*oneMinusCosT + (c*u - a*w)*sinT,
      m31 = u*w * oneMinusCosT - v*sinT,
      m32 = v*w * oneMinusCosT + u*sinT,
      m33 = w2 + (u2 + v2) * cosT,
      m34 = (c*(u2 + v2) - w*(a*u + b*v))*oneMinusCosT + (a*v - b*u)*sinT
    CS(cs.transform * SMatrix{4,4,Float64}([
          m11 m12 m13 m14;
          m21 m22 m23 m24;
          m31 m32 m33 m32;
          0 0 0 1]))
  end

rotate_vector(vector, axis, angle) =
  let s = sin(angle),
      c = cos(angle),
      x = (c + axis.x^2*(1 - c))*vector.x +
          (axis.x*axis.y*(1 - c) - axis.z*s)*vector.y +
          (axis.x*axis.z*(1 - c) + axis.y*s)*vector.z,
      y = (axis.y*axis.x*(1 - c) + axis.z*s)*vector.x +
          (c + axis.y^2*(1 - c))*vector.y +
          (axis.y*axis.z*(1 - c) - axis.x*s)*vector.z,
      z = (axis.z*axis.x*(1 - c) - axis.y*s)*vector.x +
          (axis.z*axis.y*(1 - c) + axis.x*s)*vector.y +
          (c + axis.z^2*(1 - c))*vector.z
    vxyz(x,y,z)
  end


global const world_cs = CS(Mat4f(I))
global const current_cs = Parameter(world_cs)
# Special cs for "transparent" vectors
global const null_cs = CS(Mat4f(I))

is_world_cs(cs::CS) = cs ===  world_cs

translating_current_cs(f, _dx::Real=0, _dy::Real=0, _dz::Real=0; dx::Real=_dx, dy::Real=_dy, dz::Real=_dz) =
    with(current_cs, translated_cs(current_cs(), dx, dy, dz)) do
        f()
    end

abstract type Loc end
abstract type Vec end

zero(::Type{<:Loc}) = u0()

#ideally, this should be Vector{Loc} but empty vectors of Loc are
#actually of type Vector{Any}
const Locs = Vector{<:Loc}
const Vecs = Vector{<:Vec}

#Base.==(cs0::CS, cs1::CS) = (cs0 === cs1) || (cs0.transform == cs1.transform)

#translation_matrix(x::Real, y::Real, z::Real) = CS(SMatrix{4,4,Float64}())

#Location

struct XYZ <: Loc
  x::Real
  y::Real
  z::Real
  cs::CS
  raw::Vec4f
end

# TODO add other fields, e.g.,
# Base.getproperty(p::XYZ, f::Symbol) = f === :rho ? pol_rho(p) : ...

# Basic conversions
# From tuples of Loc
convert(::Type{Locs}, ps::NTuple{N,Loc}) where {N} = collect(XYZ, ps)

# From arrays of Any. This looks like a failure in Julia type inference, particularly when
# an empty array is involved, e.g., line(vcat([xy(10,20), xy(30,40)], []))
convert(::Type{Locs}, ps::Vector{<:Any}) = collect(XYZ, ps)



show(io::IO, loc::XYZ) =
  print(io, "xyz($(loc.x),$(loc.y),$(loc.z)$(loc.cs == world_cs ? "" : ", ..."))")

#import Base.getfield, Base.Field
#getfield(p::XYZ, ::Field{:cyl_rho}) = hypot(p.x, p.y)

xyz(x::Real, y::Real, z::Real,cs::CS=current_cs()) =
  XYZ(x,y,z,cs,Vec4f(convert(Float64,x),convert(Float64,y),convert(Float64,z), 1.0))
xyz(s::Vec4f,cs::CS) =
  XYZ(s[1], s[2], s[3], cs, s)

scaled_cs(p::XYZ, x::Real, y::Real, z::Real) = xyz(p.x, p.y, p.z, scaled_cs(p.cs, x, y, z))
center_scaled_cs(p::XYZ, x::Real, y::Real, z::Real) = xyz(p.x/x, p.y/y, p.z/z, center_scaled_cs(p.cs, x, y, z))

cyl(rho::Real, phi::Real, z::Real, cs::CS=current_cs()) =
  xyz(rho*cos(phi), rho*sin(phi), z, cs)
add_cyl(p::Loc, rho::Real, phi::Real, z::Real) =
  p + vcyl(rho, phi, z, p.cs)

pol(rho::Real, phi::Real, cs::CS=current_cs()) =
  cyl(rho, phi, 0, cs)
add_pol(p::Loc, rho::Real, phi::Real) =
  p + vcyl(rho, phi, 0, p.cs)

sph(rho::Real, phi::Real, psi::Real, cs::CS=current_cs()) =
  let sin_psi = sin(psi)
    xyz(rho*cos(phi)*sin_psi, rho*sin(phi)*sin_psi, rho*cos(psi), cs)
  end
add_sph(p::Loc, rho::Real, phi::Real, psi::Real) =
  p + vsph(rho, phi, psi, p.cs)

# Vector
struct VXYZ <: Vec
    x::Real
    y::Real
    z::Real
    cs::CS
    raw::SVector{4,Float64}
end

show(io::IO, vec::VXYZ) =
  print(io, "vxyz($(vec.x),$(vec.y),$(vec.z)$(vec.cs == world_cs ? "" : ", ..."))")


vxyz(x::Real, y::Real, z::Real, cs::CS=current_cs()) =
  VXYZ(x,y,z,cs,Vec4f(convert(Float64,x),convert(Float64,y),convert(Float64,z), 0.0))
vxyz(s::Vec4f,cs::CS) = VXYZ(s[1], s[2], s[3], cs, s)

vcyl(rho::Real, phi::Real, z::Real, cs::CS=current_cs()) =
  vxyz(rho*cos(phi), rho*sin(phi), z, cs)
add_vcyl(v::Vec, rho::Real, phi::Real, z::Real) =
  v + vcyl(rho, phi, z, v.cs)

vpol(rho::Real, phi::Real, cs::CS=current_cs()) =
  vcyl(rho, phi, 0, cs)
add_vpol(v::Vec, rho::Real, phi::Real) =
  add_vcyl(v, rho, phi, 0)

vsph(rho::Real, phi::Real, psi::Real, cs::CS=current_cs()) =
  let sin_psi = sin(psi)
    vxyz(rho*cos(phi)*sin_psi, rho*sin(phi)*sin_psi, rho*cos(psi), cs)
  end
add_vsph(v::Vec, rho::Real, phi::Real, psi::Real) =
  v + vsph(rho, phi, psi, v.cs)

# Selectors

cx(p::Union{Loc,Vec}) = p.x
cy(p::Union{Loc,Vec}) = p.y
cz(p::Union{Loc,Vec}) = p.z

sph_rho(p::Union{Loc,Vec}) =
  let (x, y, z) = (p.x, p.y, p.z)
    sqrt(x*x + y*y + z*z)
  end
sph_phi(p::Union{Loc,Vec}) =
  let (x, y) = (p.x, p.y)
    0 == x == y ? 0 : mod(atan(y, x),2pi)
  end
sph_psi(p::Union{Loc,Vec}) =
  let (x, y, z) = (p.x, p.y, p.z)
    0 == x == y == z ? 0 : mod(atan(sqrt(x*x + y*y), z),2pi)
  end

cyl_rho(p::Union{Loc,Vec}) =
  let (x, y) = (p.x, p.y)
    sqrt(x*x + y*y)
  end
cyl_phi(p::Union{Loc,Vec}) = sph_phi(p)
cyl_z(p::Union{Loc,Vec}) = p.z

pol_rho = cyl_rho
pol_phi = cyl_phi

const min_norm = 1e-15

unitized(v::Vec) =
  let r = sqrt(sum(abs2, v.raw))
    @assert r > min_norm
    vxyz(v.raw./r, v.cs)
  end

in_cs(from_cs::CS, to_cs::CS) =
    to_cs === world_cs ?
        from_cs.transform :
        inv(to_cs.transform) * from_cs.transform

in_cs(p::Loc, cs::CS) =
  p.cs === cs ?
    p :
    cs === world_cs ?
      xyz(p.cs.transform * p.raw, world_cs) :
      xyz(inv(cs.transform) * p.cs.transform * p.raw, cs)

in_cs(p::Vec, cs::CS) =
  p.cs === cs ?
    p :
    cs === world_cs ?
      vxyz(p.cs.transform * p.raw, world_cs) :
      vxyz(inv(cs.transform) * p.cs.transform * p.raw, cs)

in_cs(p, q) = in_cs(p, q.cs)

in_world(p) = in_cs(p, world_cs)

export inverse_transformation
inverse_transformation(p::Loc) = xyz(0,0,0, CS(inv(translated_cs(p.cs, p.x, p.y, p.z).transform)))


cs_from_o_vx_vy_vz(o::Loc, ux::Vec, uy::Vec, uz::Vec) =
  CS(SMatrix{4,4,Float64}(ux.x, ux.y, ux.z, 0, uy.x, uy.y, uy.z, 0, uz.x, uz.y, uz.z, 0, o.x, o.y, o.z, 1))

LinearAlgebra.cross(v::Vec, w::Vec) = _cross(v.raw, in_cs(w, v.cs).raw, v.cs)
_cross(a::Vec4f, b::Vec4f, cs::CS) =
  vxyz(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1], cs)

LinearAlgebra.dot(v::Vec, w::Vec) = _dot(v.raw, in_cs(w, v.cs).raw)
_dot(a::Vec4f, b::Vec4f) =
  a[1]*b[1] + a[2]*b[2] + a[3]*b[3]

cs_from_o_vx_vy(o::Loc, vx::Vec, vy::Vec) =
  let o = in_world(o),
    vx = unitized(in_world(vx)),
    vz = unitized(cross(vx, in_world(vy)))
    cs_from_o_vx_vy_vz(o, vx, cross(vz,vx), vz)
  end

cs_from_o_vz(o::Loc, n::Vec) =
  let o = in_world(o),
      n = in_world(n),
      vx = vpol(1, sph_phi(n) + pi/2, o.cs),
      vy = unitized(cross(n, vx)),
      vz = unitized(n)
    cs_from_o_vx_vy_vz(o, vx, vy, vz)
  end

cs_from_o_phi(o::Loc, phi::Real) =
  let vx = in_world(vcyl(1, phi, 0, o.cs))
      vy = in_world(vcyl(1, phi + pi/2, 0, o.cs))
      vz = cross(vx, vy)
      o = in_world(o)
      cs_from_o_vx_vy_vz(o, vx, vy, vz)
  end

loc_from_o_vx_vy(o::Loc, vx::Vec, vy::Vec) = u0(cs_from_o_vx_vy(o, vx, vy))
loc_from_o_vz(o::Loc, vz::Vec) = u0(cs_from_o_vz(o, vz))
loc_from_o_phi(o::Loc, phi::Real) = u0(cs_from_o_phi(o, phi))

#To handle the common case
maybe_loc_from_o_vz(o::Loc, n::Vec) =
  let n = in_world(n)
    if n.x == 0 && n.y == 0
      o
    else
      loc_from_o_vz(o, n)
    end
  end

import Base.+, Base.-, Base.*, Base./, Base.length
#This is not needed!
#(+){T1,T2,T3,T4,T5,T6}(p::XYZ{T1,T2,T3},v::VXYZ{T4,T5,T6}) = xyz(p.x+v.x, p.y+v.y, p.z+v.z, p.raw+v.raw)


add_x(p::Loc, x::Real) = xyz(p.x+x, p.y, p.z, p.cs)
add_y(p::Loc, y::Real) = xyz(p.x, p.y+y, p.z, p.cs)
add_z(p::Loc, z::Real) = xyz(p.x, p.y, p.z+z, p.cs)
add_xy(p::Loc, x::Real, y::Real) = xyz(p.x+x, p.y+y, p.z, p.cs)
add_xz(p::Loc, x::Real, z::Real) = xyz(p.x+x, p.y, p.z+z, p.cs)
add_yz(p::Loc, y::Real, z::Real) = xyz(p.x, p.y+y, p.z+z, p.cs)
add_xyz(p::Loc, x::Real, y::Real, z::Real) = xyz(p.x+x, p.y+y, p.z+z, p.cs)

(+)(a::XYZ, b::VXYZ) = xyz(a.raw + in_cs(b, a.cs).raw, a.cs)
(+)(a::VXYZ, b::XYZ) = xyz(a.raw + in_cs(b, a.cs).raw, a.cs)
(+)(a::VXYZ, b::VXYZ) = vxyz(a.raw + in_cs(b, a.cs).raw, a.cs)
(-)(a::XYZ, b::VXYZ) = xyz(a.raw - in_cs(b, a.cs).raw, a.cs)
(-)(a::VXYZ, b::VXYZ) = vxyz(a.raw - in_cs(b, a.cs).raw, a.cs)
(-)(a::XYZ, b::XYZ) = vxyz(a.raw - in_cs(b, a.cs).raw, a.cs)
(-)(a::VXYZ) = vxyz(-a.raw, a.cs)
(*)(a::VXYZ, b::Real) = vxyz(a.raw * b, a.cs)
(*)(a::Real, b::VXYZ) = vxyz(a * b.raw, b.cs)
(/)(a::VXYZ, b::Real) = vxyz(a.raw / b, a.cs)

norm(v::Vec) = norm(v.raw)
length(v::Vec) = norm(v.raw)

min_loc(p::Loc, q::Loc) =
    xyz(min.(p.raw, in_cs(q, p.cs).raw), p.cs)
max_loc(p::Loc, q::Loc) =
    xyz(max.(p.raw, in_cs(q, p.cs).raw), p.cs)

distance(p::XYZ, q::XYZ) = norm((in_world(q)-in_world(p)).raw)

u0(cs=current_cs())   = xyz(0,0,0,cs)
ux(cs=current_cs())   = xyz(1,0,0,cs)
uy(cs=current_cs())   = xyz(0,1,0,cs)
uz(cs=current_cs())   = xyz(0,0,1,cs)
uxy(cs=current_cs())  = xyz(1,1,0,cs)
uyz(cs=current_cs())  = xyz(0,1,1,cs)
uxz(cs=current_cs())  = xyz(1,0,1,cs)
uxyz(cs=current_cs()) = xyz(1,1,1,cs)

x(x::Real=1,cs=current_cs()) = xyz(x,0,0,cs)
y(y::Real=1,cs=current_cs()) = xyz(0,y,0,cs)
z(z::Real=1,cs=current_cs()) = xyz(0,0,z,cs)
xy(x::Real=1,y::Real=1,cs=current_cs()) = xyz(x,y,0,cs)
yz(y::Real=1,z::Real=1,cs=current_cs()) = xyz(0,y,z,cs)
xz(x::Real=1,z::Real=1,cs=current_cs()) = xyz(x,0,z,cs)

uvx(cs=current_cs())   = vxyz(1,0,0,cs)
uvy(cs=current_cs())   = vxyz(0,1,0,cs)
uvz(cs=current_cs())   = vxyz(0,0,1,cs)
uvxy(cs=current_cs())  = vxyz(1,1,0,cs)
uvyz(cs=current_cs())  = vxyz(0,1,1,cs)
uvxz(cs=current_cs())  = vxyz(1,0,1,cs)
uvxyz(cs=current_cs()) = vxyz(1,1,1,cs)

vx(x::Real=1,cs=current_cs()) = vxyz(x,0,0,cs)
vy(y::Real=1,cs=current_cs()) = vxyz(0,y,0,cs)
vz(z::Real=1,cs=current_cs()) = vxyz(0,0,z,cs)
vxy(x::Real=1,y::Real=1,cs=current_cs()) = vxyz(x,y,0,cs)
vyz(y::Real=1,z::Real=1,cs=current_cs()) = vxyz(0,y,z,cs)
vxz(x::Real=1,z::Real=1,cs=current_cs()) = vxyz(x,0,z,cs)

position_and_height(p, q) = loc_from_o_vz(p, q - p), distance(p, q)

regular_polygon_vertices(edges::Integer=3, center::Loc=u0(), radius::Real=1, angle::Real=0, is_inscribed::Bool=true) = begin
  r = is_inscribed ? radius : radius/cos(pi/edges)
  [center + vpol(r, a, center.cs) for a in division(angle, angle + 2*pi, edges, false)]
end

intermediate_loc(p::Loc, q::Loc, f::Real=0.5) =
  if p.cs == q.cs
    p+(q-p)*f
  else
    o = intermediate_loc(in_world(p), in_world(q), f)
    v_x = in_world(vx(1, p.cs))*(1-f) + in_world(vx(1, q.cs))*f
    v_y = in_world(vy(1, p.cs))*(1-f) + in_world(vy(1, q.cs))*f
    loc_from_o_vx_vy(o, v_x, v_y)
  end

# Metaprogramming

meta_program(x::Any) = x # literals might be self evaluating
meta_program(x::Real) = round(x,sigdigits=8)
meta_program(x::Bool) = x
meta_program(x::Vector) = Expr(:vect, map(meta_program, x)...)
meta_program(p::Loc) =
    if cz(p) == 0
        Expr(:call, :xy, meta_program(cx(p)), meta_program(cy(p)))
    else
        Expr(:call, :xyz, meta_program(cx(p)), meta_program(cy(p)), meta_program(cz(p)))
    end
meta_program(v::Vec) =
    if cz(v) == 0
        Expr(:call, :vxy, meta_program(cx(v)), meta_program(cy(v)))
    else
        Expr(:call, :vxyz, meta_program(cx(v)), meta_program(cy(v)), meta_program(cz(v)))
    end

# Conversions
# We could accept some nice conversions
# convert(::Type{Loc}, t::Tuple{Real,Real,Real}) = xyz(t[1], t[2], t[3])


# Integration in standard protocols

# iteration for destructuring into components
iterate(v::Vec) = iterate(v.raw)
iterate(v::Vec, state) = iterate(v.raw, state)

iterate(v::Loc) = iterate(v.raw)
iterate(v::Loc, state) = iterate(v.raw, state)

# Utilities
export trig_center, trig_normal, quad_center, quad_normal, vertices_center, vertices_normal, iterate_quads

trig_center(p0, p1, p2) =
  xyz((p0.x+p1.x+p2.x)/3, (p0.y+p1.y+p2.y)/3, (p0.z+p1.z+p2.z)/3, p0.cs)

trig_normal(p0, p1, p2) =
  polygon_normal([p1 - p0, p2 - p1, p0 - p2])

quad_center(p0, p1, p2, p3) =
  intermediate_loc(intermediate_loc(p0, p2), intermediate_loc(p1, p3))

quad_normal(p0, p1, p2, p3) =
  let p0 = in_world(p0),
      p1 = in_world(p1),
      p2 = in_world(p2),
      p3 = in_world(p3)
    polygon_normal([p1 - p0, p2 - p1, p3 - p2, p0 - p3])
  end

vertices_center(pts) =
  let pts = map(in_world, pts),
      n=length(pts),
      xs=[cx(p) for p in pts],
      ys=[cy(p) for p in pts],
      zs=[cz(p) for p in pts]
    xyz(sum(xs)/n, sum(ys)/n, sum(zs)/n, world_cs)
  end

vertices_normal(ps) =
  let ps = map(in_world, ps)
    polygon_normal(p-q for (p,q) in zip(ps, drop(cycle(ps), 1)))
  end

polygon_normal(vs) =
  unitized(
    sum(
      cross(v0,v1)
      for (v0,v1) in zip(vs, drop(cycle(vs), 1))))

iterate_quads(f, ptss) =
  [[f(p0, p1, p2, p3)
    for (p0, p1, p2, p3)
    in zip(pts0[1:end-1], pts1[1:end-1], pts1[2:end], pts0[2:end])]
    for (pts0, pts1)
    in zip(ptss[1:end-1], ptss[2:end])]


# We need to implement smooth walks along curves
# Using just the Frenet frame is not adequate as
# changes in curvature cause it to suddenly change
# direction.

# The technique we use to solve this is based on
# rotation minimizing frames (RMF) presented in
# "Computation of Rotation Minimizing Frames"
# (Wenping Wang, Bert Jüttler, Dayue Zheng, and Yang Liu, 2008)
# Regarding the paper notation, t = vz, r = -vx, s = vy
export rotation_minimizing_frames

rotation_minimizing_frames(frames) =
  rotation_minimizing_frames(
    frames[1],
    [in_world(frame) for frame in frames],
    [in_world(vz(1, frame.cs)) for frame in frames])

#=
  let new_frames = [frames[1]]
    for x1 in drop(frames, 1)
      # Reflection of x0 tangent and axis onto x1
      # using reflection plane located between x0 and x1
      let x0 = new_frames[end],
          v1 = in_world(x1) - in_world(x0),
          c1 = dot(v1, v1),
          r0 = in_world(vx(1, x0.cs)),
          t0 = in_world(vz(1, x0.cs)),
          ril = r0 - v1*(2/c1*dot(v1,r0)),
          til = t0 - v1*(2/c1*dot(v1,t0)),
          # Reflection on a plane at x1, aligning the frame
          # tangent with the curve tangent
          t1 = in_world(vz(1, x1.cs)),
          v2 = t1 - til,
          c2 = dot(v2, v2),
          r1 = ril - v2*(2/c2*dot(v2, ril)),
          s1 = cross(t1, r1)
        push!(new_frames, loc_from_o_vx_vy(x1, r1, s1))
      end
    end
    new_frames
  end
=#

rotation_minimizing_frames(u0, xs, ts) =
  let ri = in_world(vy(1, u0.cs)),
      new_frames = [loc_from_o_vx_vy(xs[1], cross(ri, ts[1]), ri)]
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
          sii = cross(rii, tii),
          uii = loc_from_o_vx_vy(xii, sii, rii)
        push!(new_frames, uii)
        ri = rii
      end
    end
    new_frames
  end

#=
TODO: Consider adding vectors without frame of reference:

Δx(x) = vx(x)

p + Δx(5)
p + Δxy(1,2)
Addition would respect the frame of reference of p

Another hypotesis is to write

p + _x(5)
p + _xy(1,2)

or even

p + dx(5)
p + dxy(1,2)

For the moment, I'll choose this one
=#


# Broadcasting

broadcastable(p::Loc) = Ref(p)
broadcastable(v::Vec) = Ref(v)

# equality

Base.isequal(p::Loc, q::Loc) =
  let wp = in_world(p),
      wq = in_world(q)
    isequal(wp.x, wq.x) && isequal(wp.y, wq.y) && isequal(wp.z, wq.z)
  end

# angle between

angle_between(v1, v2) =
  let v1 = in_world(v1),
      v2 = in_world(v2)
    acos(dot(v1, v2)/(norm(v1)*norm(v2)))
  end

################################################################################
# To embed Khepri, it becomes useful to convert entities into 'raw' data

raw_point(v::Union{XYZ, VXYZ}) =
  let o = in_world(v)
    (float(o.x), float(o.y), float(o.z))
  end

raw_plane(v::XYZ) =
  let o = in_world(v),
      vx = in_world(vx(1, c.cs))
      vy = in_world(vy(1, c.cs))
    (float(o.x), float(o.y), float(o.z),
     float(vx.x), float(vx.y), float(vx.z),
     float(vy.x), float(vy.y), float(vy.z))
  end

################################################################################
export acw_vertices
acw_vertices(vs) =
  angle_between(vertices_normal(vs), vz()) < pi/4 ?
    vs : reverse(vs)
