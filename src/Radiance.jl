#=

A backend for Radiance

To visualize RAD files, check
https://www.ladybug.tools/spider-rad-viewer/rad-viewer/r7/rad-viewer.html

=#
export Radiance,
       radiance,
       save_rad,
       @save_rad,
       RadianceMaterial,
       radiance_material,
       radiance_plastic_material,
       radiance_metal_material,
       radiance_glass_material,
       radiance_light_material,
       radiance_dielectric_material

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
#    with(current_backend, autocad) do
#      polygon(vertices)
#      #surface_polygon(vertices)
#    end
    println(io, modifier, " ", "polygon", " ", id)
    println(io, 0) #0 strings
    println(io, 0) #0 ints
    println(io, 3*length(vertices))
    for v in vertices println(io, " ", v.x, " ", v.y, " ", v.z) end
  end

write_rad_sphere(io::IO, modifier, id, center, radius) =
  write_rad_primitive(io, modifier, "sphere", id, [], [],
    [center.x, center.y, center.z, radius])

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

write_rad_divided_quad(io::IO, modifier, id, sub_id, v0, v1, v2, v3) =
  begin
    println(io, modifier, " ", "polygon", " ", id, sub_id, 0)
    println(io, 0) #0 strings
    println(io, 0) #0 ints
    println(io, 9)
    println(io, " ", v0.x, " ", v0.y, " ", v0.z)
    println(io, " ", v1.x, " ", v1.y, " ", v1.z)
    println(io, " ", v2.x, " ", v2.y, " ", v2.z)
    println(io, modifier, " ", "polygon", " ", id, sub_id, 1)
    println(io, 0) #0 strings
    println(io, 0) #0 ints
    println(io, 9)
    println(io, " ", v0.x, " ", v0.y, " ", v0.z)
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
    write_rad_quad(io, modifier, id, "face5", p4, p5, p6, p7)
  end

#=
Materials
=#

struct RadianceMaterial
  name::String
  type::String
  red::Real
  green::Real
  blue::Real
  specularity::Union{Real, Nothing}
  roughness::Union{Real, Nothing}
  transmissivity::Union{Real, Nothing}
  transmitted_specular::Union{Real, Nothing}
end

radiance_string(m::RadianceMaterial) =
  isnothing(m.transmissivity) && isnothing(m.transmitted_specular) ?
    isnothing(m.specularity) && isnothing(m.roughness) ?
      "void $(m.type) $(m.name)\n0\n0\n3 $(m.red) $(m.green) $(m.blue)\n" :
      "void $(m.type) $(m.name)\n0\n0\n5 $(m.red) $(m.green) $(m.blue) $(m.specularity) $(m.roughness)\n" :
    "void $(m.type) $(m.name)\n0\n0\n7 $(m.red) $(m.green) $(m.blue) $(m.specularity) $(m.roughness) $(m.transmissivity) $(m.transmitted_specular)\n"

write_rad_material(io::IO, material::RadianceMaterial) =
  write(io, radiance_string(material))

#Base.show(io::IO, mat::RadianceMaterial) =
#  write_rad_material(io, mat)

radiance_material(name::String, type::String;
                  gray::Real=0.3,
                  red::Real=gray, green::Real=gray, blue::Real=gray,
                  specularity=nothing, roughness=nothing,
                  transmissivity=nothing, transmitted_specular=nothing) =
  RadianceMaterial(name, type,
                   red, green, blue,
                   specularity, roughness,
                   transmissivity, transmitted_specular)

radiance_light_material(name::String; args...) =
  radiance_material(name, "light"; args...)

radiance_plastic_material(name::String; args...) =
  radiance_material(name, "plastic"; specularity=0, roughness=0, args...)

radiance_metal_material(name::String; args...) =
  radiance_material(name, "metal"; specularity=0, roughness=0, args...)

radiance_glass_material(name::String; args...) =
  radiance_material(name, "glass"; args...)

radiance_dielectric_material(name::String; args...) =
  radiance_material(name, "dielectric"; args...)

#Some pre-defined materials
radiance_material_white = radiance_plastic_material("white", gray=1.0) #
radiance_generic_ceiling_70 = radiance_plastic_material("GenericCeiling_70", gray=0.7)
radiance_generic_ceiling_80 = radiance_plastic_material("GenericCeiling_80", gray=0.8)
radiance_generic_ceiling_90 = radiance_plastic_material("HighReflectanceCeiling_90", gray=0.9)
radiance_generic_floor_20 = radiance_plastic_material("GenericFloor_20", gray=0.2)
radiance_generic_interior_wall_50 = radiance_plastic_material("GenericInteriorWall_50", gray=0.5)
radiance_generic_interior_wall_70 = radiance_plastic_material("GenericInteriorWall_70", gray=0.7)
radiance_generic_furniture_50 = radiance_plastic_material("GenericFurniture_50", gray=0.5)
radiance_outside_facade_30 = radiance_plastic_material("OutsideFacade_30", gray=0.3)
radiance_outside_facade_35 = radiance_plastic_material("OutsideFacade_35", gray=0.35)
radiance_generic_glass_80 = radiance_glass_material("Glass_80", gray=0.8)
radiance_generic_metal = radiance_metal_material("SheetMetal_80", gray=0.8)

default_radiance_material = Parameter(radiance_generic_metal)

export radiance_material_white,
       radiance_generic_ceiling_70,
       radiance_generic_ceiling_80,
       radiance_generic_ceiling_90,
       radiance_generic_floor_20,
       radiance_generic_interior_wall_50,
       radiance_generic_interior_wall_70,
       radiance_generic_furniture_50,
       radiance_outside_facade_30,
       radiance_outside_facade_35,
       radiance_generic_glass_80,
       radiance_generic_metal,
       default_radiance_material

####################################################
# Sky models

const radiance_extra_sky_rad_contents = """
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
4 1 .8 .5 0
ground_glow source ground
0
0
4 0 0 -1 180
"""

radiance_utah_sky_string(date, latitude, longitude, meridian, turbidity, withsun) =
  let _2d(n) = lpad(n, 2, '0')
    # Positive longitude is east, negative is west but genutahsky considers the opposite
    "!genutahsky $(_2d(month(date))) $(_2d(day(date))) $(_2d(hour(date)+minute(date)/60)) -y $(year(date)) -t $(turbidity) -a $(latitude) -o $(360+longitude) -m $(360+meridian)" *
    "\n" *
    radiance_extra_sky_rad_contents
  end

radiance_cie_sky_string(date, latitude, longitude, meridian, turbidity, withsun) =
  let _2d(n) = lpad(n, 2, '0')
    # Positive longitude is east, negative is west but gensky considers the opposite
    "!gensky $(_2d(month(date))) $(_2d(day(date))) $(_2d(hour(date))):$(_2d(minute(date))) " *
    "$(withsun ? "+s" : "-s") -a $(latitude) -o $(-longitude) -m $(-meridian) -t $(turbidity)\n" *
    radiance_extra_sky_rad_contents
  end

radiance_cie_sky_string(altitude, azimuth, turbidity, withsun) =
  let _2d(n) = lpad(n, 2, '0')
    "!gensky -ang $(altitude) $(azimuth) " *
    "$(withsun ? "+s" : "-s") -t $(turbidity)\n" *
    radiance_extra_sky_rad_contents
  end

#=
Simulations need to be done on a temporary folder, so that we can have multiple
simulations running at the same time.
=#

radiance_simulation_path() =
 let (path, io) = mktemp(mktempdir(tempdir(), prefix="Radiance_"))
   close(io)
   path
 end

export radiance_folder
const radiance_folder = Parameter("C:/Program Files/Radiance/bin/")

radiance_folder("C:/DIVA/Radiance/bin_64/")

radiance_cmd(cmd::AbstractString) = radiance_folder() * cmd

##########################################
# oconv

radiance_oconv(radpath, matpath=missing, skypath=missing) =
  let octpath = path_replace_suffix(radpath, ".oct"),
      pipe = ismissing(matpath) ?
             `$(radiance_cmd("oconv")) $radpath` :
             ismissing(skypath) ?
            `$(radiance_cmd("oconv")) $matpath $radpath` :
            `$(radiance_cmd("oconv")) $matpath $skypath $radpath`
    run(pipeline(pipe, stdout=octpath))
    octpath
  end

##########################################
# rview and rpict share a lot of parameters that define the rendering quality
#=
rpict_abbrev_parameter = Dict(
  "ambient_bounces"=>"ab",
  "ambient_divisions"=>"ad",
  "ambient_supersamples"=>"as",
  "ambient_resolution"=>"ar",
  "ambient_accuracy"=>"aa",
  "pixel_sampling"=>"ps",
  "pixel_tolerance"=>"pt",
  "pixel_jitter"=>"pj",
  "direct_jitter"=>"dj",
  "direct_sampling"=>"ds",
  "direct_threshold"=>"dt",
  "direct_certainty"=>"dc",
  "direct_sec_relays"=>"dr",
  "direct_presamp_density"=>"dp",
  "specular_threshold"=>"st",
  "limit_reflections"=>"lr",
  "limit_weight"=>"lw",
  "specular_sampling"=>"ss"
)


rpict_number_parameters = {
    ("ambient_bounces", "ab", int, "default_value": 2, "values": [2, 3, 6]},
    ("ambient_divisions", "ad", int, "default_value": 512, "values": [512, 2048, 4096]},
    ("ambient_supersamples", "as", int, "default_value": 128, "values": [128, 2048, 4096]},
    ("ambient_resolution", "ar", int, "default_value": 16, "values": [16, 64, 128]},
    ("ambient_accuracy", "aa", float, "default_value": .25, "values": [.25, .2, .1]},
    ("pixel_sampling", "ps", int, "default_value": 8, "values": [8, 4, 2]},
    ("pixel_tolerance", "pt", float, "default_value": .15, "values": [.15, .10, .05]},
    ("pixel_jitter", "pj", float, "default_value": .6, "values": [.6, .9, .9]},
    ("direct_jitter", "dj", float, "default_value": 0, "values": [0, .5, 1]},
    ("direct_sampling", "ds", float, "default_value": .5, "values": [.5, .25, .05]},
    ("direct_threshold", "dt", float, "default_value": .5, "values": [.5, .25, .15]},
    ("direct_certainty", "dc", float, "default_value": .25, "values": [.25, .5, .75]},
    ("direct_sec_relays", "dr", int, "default_value": 0, "values": [0, 1, 3]},
    ("direct_presamp_density", "dp", int, "default_value": 64, "values": [64, 256, 512]},
    ("specular_threshold", "st", float, "default_value": .85, "values": [.85, .5, .15]},
    ("limit_reflections", "lr", int, "default_value": 4, "values": [4, 6, 8]},
    ("limit_weight", "lw", float, "default_value": .05, "values": [.05, .01, .005]},
    ("specular_sampling", "ss", "dscrip": "specular sampling", "type": float, "default_value": 0, "values": [0, .7, 1]}}
=#

##########################################
# rview

radiance_rview(octpath, camera, target, lens, light=(1,1,1)) =
  let p = camera,
      v = target-camera,
      (h_angle, v_angle) = view_angles(lens)
    run(`$(radiance_cmd("rvu"))
        -ab 2
        -vp $(p.x) $(p.y) $(p.z)
        -vd $(v.x) $(v.y) $(v.z)
        -vh $(h_angle)
        -vv $(v_angle)
        -av $(light[1]) $(light[2]) $(light[3])
        $octpath`,
        wait=false)
  end

##########################################
# rview

radiance_rpict(octpath, camera, target, lens) =
  let picpath = path_replace_suffix(octpath, ".pic"),
      p = camera,
      v = target-camera,
      (h_angle, v_angle) = view_angles(lens)
    run(pipeline(`$(radiance_cmd("rpict"))
        -vp $(p.x) $(p.y) $(p.z)
        -vd $(v.x) $(v.y) $(v.z)
        -vh $(h_angle)
        -vv $(v_angle)
        -x $(render_width())
        -y $(render_height())
        -aa 0.1 -ab 6 -ad 4096
        -dc 0.75 -st 0.15 -lw 0.005
        -as 4096 -ar 128 -lr 8
        -dt 0.15 -dr 3 -ds 0.05 -dp 512
        $octpath`,
        #`$(radiance_cmd("pfilt"))`,
        stdout=picpath))
    #run(`perl $(radiance_cmd("falsecolor.pl")) $picpath`, wait=false)
    run(`$(radiance_cmd("wxFalseColor")) $picpath`, wait=false)
  end

diva_render(octpath, camera, target, lens) =
  let ovepath = path_replace_suffix(octpath, ".overture"),
      ambpath = path_replace_suffix(octpath, ".amb"),
      unfpath = path_replace_suffix(octpath, ".unf"),
      picpath = path_replace_suffix(octpath, ".pic"),
      p = camera,
      v = target-camera,
      (h_angle, v_angle) = view_angles(lens),
      viewstr = "-vp $(p.x) $(p.y) $(p.z) -vd $(v.x) $(v.y) $(v.z) -vu 0 0 1",
      sizestr = "-vh $(h_angle) -vv $(v_angle)",
      resolutionstr = "-x $(2*render_width()) -y $(2*render_height())",
      basestr = "-vtv $(viewstr) $(sizestr) -vs 0 -vl 0 $(resolutionstr)",
      arg1str = "-ps 4 -pt .10 -pj .9 -dj .5 -ds .25 -dt .25 -dc .50 -dr 1 -dp 256 -st .50 -ab 3 -aa .2 -ar 256 -ad 2048 -as 1024 -lr 6 -lw .010",
      arg2str = "-ps 2 -pt .05 -pj .9 -dj .7 -ds .15 -dt .05 -dc .75 -dr 3 -dp 512 -st .15 -ab 4 -aa .1 -ar 512 -ad 2048 -as 1024 -lr 8 -lw .005",
      run1str = split("$(basestr) $(arg1str)", ' '),
      run2str = split("$(basestr) $(arg2str)", ' ')
    println(`$(radiance_cmd("rpict")) $(run1str) -af $(ambpath) $(octpath)`)
    @time run(pipeline(`$(radiance_cmd("rpict")) $(run1str) -af $(ambpath) $(octpath)`, stdout=ovepath))
    #rm(ovepath)
    println(`$(radiance_cmd("rpict")) $(run2str) -af $(ambpath) $(octpath)`)
    @time run(pipeline(`$(radiance_cmd("rpict")) $(run2str) -af $(ambpath) $(octpath)`, stdout=unfpath))
    println(`$(radiance_cmd("pfilt")) -r .6 -x /2 -y /2 $(unfpath)`)
    @time run(pipeline(`$(radiance_cmd("pfilt")) -r .6 -x /2 -y /2 $(unfpath)`, stdout=picpath))
    run(`$(radiance_cmd("wxFalseColor")) $(picpath)`, wait=false)
  end

#= From DIVA

rpict -t 15 -vtv -vp 9 19.3 4.4 -vd 0 -4.3 0.3 -vu 0 0 1 -vh 80.4462463224827 -vv 57.6922631046869 -vs 0 -vl 0 -af test01.amb -x 1600 -y 1200 -ps 4 -pt .10 -pj .9 -dj .5 -ds .25 -dt .25 -dc .5 -dr 1 -dp 256 -st .5 -ab 3 -aa .2 -ar 256 -ad 2048 -as 1024 -lr 6 -lw .01 test01.oct  1>test01_Perspective41.overture
del test01_Perspective41.overture
rpict -t 15 -vtv -vp 9 19.3 4.4 -vd 0 -4.3 0.3 -vu 0 0 1 -vh 80.4462463224827 -vv 57.6922631046869 -vs 0 -vl 0 -af test01.amb -x 1600 -y 1200 -ps 4 -pt .10 -pj .9 -dj .5 -ds .25 -dt .25 -dc .5 -dr 1 -dp 256 -st .5 -ab 3 -aa .2 -ar 256 -ad 2048 -as 1024 -lr 6 -lw .01 test01.oct  1>test01_Perspective41.unf
rpict -t 15 -vtv -vp 9 19.3 4.4 -vd 0 -4.3 0.3 -vu 0 0 1 -vh 80.4462463224827 -vv 57.6922631046869 -vs 0 -vl 0 -af test01.amb -x 1600 -y 1200 -ps 2 -pt .05 -pj .9 -dj .7 -ds .15 -dt .05 -dc .75 -dr 3 -dp 512 -st .15 -ab 4 -aa .1 -ar 512 -ad 2048 -as 1024 -lr 8 -lw .005 test01.oct  1>test01_Perspective41.unf



pfilt -r .6 -x /2 -y /2 test01_Perspective41.unf 1>test01_Perspective41.pic
del test01_Perspective41.unf
del test01.amb



rpict
-t 15  # report progress
-vtv #perpective view
-vp 9 19.3 4.4 # camera
-vd 0 -4.3 0.3 # target
-vu 0 0 1 # up
-vh 80.4462463224827
-vv 57.6922631046869
-vs 0
-vl 0
-af test01.amb
-x 1600 -y 1200
-ps 2
-pt .05
-pj .9
-dj .7
-ds .15
-dt .05
-dc .75
-dr 3
-dp 512
-st .15
-ab 4
-aa .1
-ar 512
-ad 2048
-as 1024
-lr 8
-lw .005
test01.oct  1>test01_Perspective41.overture

=#

# Test1
#=
room_rad = radiance_simulation_path()

open(room_rad, "w") do out
  write_rad_material(out, radiance_light_material("bright", gray=100))
  write_rad_material(out, radiance_plastic_material("red_plastic", red=.7, green=.05, blue=.05, specularity=.5, roughness=.5))
  write_rad_sphere(out, "bright", "fixture", xyz(2, 1, 1.5), .125)
  write_rad_sphere(out, "red_plastic", "ball", xyz(.7, 1.125, .625), .125)
end
radiance_rview(radiance_oconv(room_rad), xyz(2.25, .375, 1), xyz(2.25, .375, 1) + vxyz(-.25, .125, -.125), 50)
=#

#=
# Test2

open(room_rad, "w") do out
  write_rad_material(out, radiance_light_material("bright", gray=100))
  write_rad_material(out, radiance_plastic_material("red_plastic", red=.7, green=.05, blue=.05, specularity=.5, roughness=.5))
  write_rad_material(out, radiance_plastic_material("gray_paint", gray=.5))
  write_rad_sphere(out, "bright", "fixture", xyz(2, 1, 1.5), .125)
  write_rad_sphere(out, "red_plastic", "ball", xyz(.7, 1.125, .625), .125)
  write_rad_quad(out, "gray_paint", "room", "1", xyz(0,0,0), xyz(0,0,1.75), xyz(3,0,1.75), xyz(3,0,0))
  write_rad_quad(out, "gray_paint", "room", "2", xyz(0,0,0), xyz(0,2,0), xyz(0,2,1.75), xyz(0,0,1.75))
  write_rad_quad(out, "gray_paint", "room", "3", xyz(0,0,0), xyz(3,0,0), xyz(3,2,0), xyz(0,2,0))
  write_rad_quad(out, "gray_paint", "room", "4", xyz(3,2,1.75), xyz(0,2,1.75), xyz(0,2,0), xyz(3,2,0))
  write_rad_quad(out, "gray_paint", "room", "5", xyz(3,2,1.75), xyz(3,2,0), xyz(3,0,0), xyz(3,0,1.75))
  write_rad_quad(out, "gray_paint", "room", "6", xyz(3,2,1.75), xyz(3,0,1.75), xyz(0,0,1.75), xyz(0,2,1.75))
end

radiance_rview(radiance_oconv(room_rad), xyz(2.25, .375, 1), xyz(2.25, .375, 1) + vxyz(-.25, .125, -.125), 50, (.5, .5, .5))

# Test2

open(room_rad, "w") do out
  write_rad_material(out, radiance_light_material("bright", gray=100))
  write_rad_material(out, radiance_plastic_material("red_plastic", red=.7, green=.05, blue=.05, specularity=.5, roughness=.5))
  write_rad_material(out, radiance_plastic_material("gray_paint", gray=.5))
  write_rad_sphere(out, "bright", "fixture", xyz(2, 1, 1.5), .125)
  write_rad_sphere(out, "red_plastic", "ball", xyz(.7, 1.125, .625), .125)
  write_rad_quad(out, "gray_paint", "room", "1", xyz(0,0,0), xyz(0,0,1.75), xyz(3,0,1.75), xyz(3,0,0))
  write_rad_quad(out, "gray_paint", "room", "2", xyz(0,0,0), xyz(0,2,0), xyz(0,2,1.75), xyz(0,0,1.75))
  write_rad_quad(out, "gray_paint", "room", "3", xyz(0,0,0), xyz(3,0,0), xyz(3,2,0), xyz(0,2,0))
  write_rad_quad(out, "gray_paint", "room", "4", xyz(3,2,1.75), xyz(0,2,1.75), xyz(0,2,0), xyz(3,2,0))
  write_rad_quad(out, "gray_paint", "room", "5", xyz(3,2,1.75), xyz(3,2,0), xyz(3,0,0), xyz(3,0,1.75))
  write_rad_quad(out, "gray_paint", "room", "6", xyz(3,2,1.75), xyz(3,0,1.75), xyz(0,0,1.75), xyz(0,2,1.75))
  write_rad_material(out, radiance_plastic_material("blue_plastic", red=.1, green=.1, blue=.6, specularity=.05, roughness=.1))
  write_rad_quad(out, "blue_plastic", "box", "1540", xyz(0.98,0.88,0), xyz(0.98,0.88,0.5), xyz(0.5,0.75,0.5), xyz(0.5,0.75,0))
  write_rad_quad(out, "blue_plastic", "box", "4620", xyz(0.5,0.75,0.5), xyz(0.37,1.23,0.5), xyz(0.37,1.23,0), xyz(0.5,0.75,0))
  write_rad_quad(out, "blue_plastic", "box", "2310", xyz(0.37,1.23,0), xyz(0.85,1.36,0), xyz(0.98,0.88,0), xyz(0.5,0.75,0))
  write_rad_quad(out, "blue_plastic", "box", "3267", xyz(0.85,1.36,0), xyz(0.37,1.23,0), xyz(0.37,1.23,0.5), xyz(0.85,1.36,0.5))
  write_rad_quad(out, "blue_plastic", "box", "5137", xyz(0.98,0.88,0.5), xyz(0.98,0.88,0), xyz(0.85,1.36,0), xyz(0.85,1.36,0.5))
  write_rad_quad(out, "blue_plastic", "box", "6457", xyz(0.37,1.23,0.5), xyz(0.5,0.75,0.5), xyz(0.98,0.88,0.5), xyz(0.85,1.36,0.5))
  write_rad_material(out, radiance_metal_material("chrome", gray=.8, specularity=.9, roughness=0))
  write_rad_cylinder(out, "chrome", "fixture_support", xyz(2,1,1.5), .05, xyz(2,1,1.75))
end

radiance_rview(radiance_oconv(room_rad), xyz(2.25, .375, 1), xyz(2.25, .375, 1) + vxyz(-.25, .125, -.125), 50, (.5, .5, .5))

=#

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

mutable struct RadianceBackend{K,T,LOD} <: LazyBackend{K,T}
  shapes::Shapes
  shape_material::Dict{Shape,RadianceMaterial}
  materials::Dict
  sky::String
  buffer::LazyParameter{IOBuffer}
  camera::Loc
  target::Loc
  lens::Real
  sun_altitude::Real
  sun_azimuth::Real
  count::Integer
end

const Radiance{LOD} = RadianceBackend{RadianceKey, RadianceId, LOD}
# Traits
has_boolean_ops(::Type{Radiance}) = HasBooleanOps{false}()
#eager_realize(::Type{Radiance}) = EagerRealize{false}()

save_shape!(b::Radiance, s::Shape) =
  begin
    push!(b.shapes, s)
    b.shape_material[s] = default_radiance_material()
    s
  end

#=
The Radiance backend cannot realize shapes immediately, only when requested.
=#

void_ref(b::Radiance) = RadianceNativeRef(-1)

const radiance =
  Radiance{500}(Shape[],
                Dict{Shape,RadianceMaterial}(),
                Dict(),
                radiance_utah_sky_string(DateTime(2020, 9, 21, 10, 0, 0), 39, 9, 0, 5, true),
                LazyParameter(IOBuffer, IOBuffer),
                xyz(10,10,10),
                xyz(0,0,0),
                35,
                90,
                0,
                0)

buffer(b::Radiance) = b.buffer()
get_material(b::Radiance, key) = get!(b.materials, key, key.name)
get_material(b::Radiance, s::Shape) = get_material(b, get(b.shape_material, s, default_radiance_material()))
next_id(b::Radiance) =
  begin
      b.count += 1
      b.count - 1
  end

save_rad(path::String, b::Radiance=current_backend()) =
  let buf = b.buffer()
    take!(buf)
    for s in b.shapes
      realize(b, s)
    end
    open(path, "w") do out
      write(out, String(take!(buf)))
    end
  end

macro save_rad()
  :(save_rad(splitext($(string(__source__.file)))[1]*".rad"))
end

#

set_view(camera::Loc, target::Loc, lens::Real, b::Radiance) =
  begin
#    set_view(camera, target, lens, autocad)
    b.camera = camera
    b.target = target
    b.lens = lens
  end

get_view(b::Radiance) =
  b.camera, b.target, b.lens

backend_realistic_sky(b::Radiance, date, latitude, longitude, meridian, turbidity, withsun) =
  b.sky = radiance_utah_sky_string(date, latitude, longitude, meridian, turbidity, withsun)

backend_realistic_sky(b::Radiance, altitude, azimuth, turbidity, withsun) =
  b.sky = radiance_cie_sky_string(altitude, azimuth, turbidity, withsun)

#
delete_all_shapes(b::Radiance) =
  begin
#    delete_all_shapes(autocad)
    (empty!(b.shapes); empty!(b.materials); empty!(b.shape_material); b.count = 0; nothing)
  end

realize(b::Radiance, s::Sphere) =
  let id = next_id(b),
      mod = get_material(b, s),
      kind = "sphere",
      buf = buffer(b)
    write_rad_sphere(buf, mod, id, in_world(s.center), s.radius)
    id
  end

#=
backend(radiance)
delete_all_shapes()
sphere()
@save_rad()
=#

  #=
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
=#
realize(b::Radiance, s::Box) =
  let id = next_id(b),
      mod = get_material(b, s),
      kind = "box",
      buf = buffer(b)
    write_rad_box(buf, #="mat$(kind)$(mod)"=# mod, id, in_world(s.c), s.dx, s.dy, s.dz)
    id
  end

realize(b::Radiance, s::SurfaceGrid) =
  let buf = buffer(b),
      mod = get_material(b, s),
      id = next_id(b),
      n = 0
    quad_grid(s.points, s.closed_u, s.closed_v) do p1, p2, p3, p4
      write_rad_divided_quad(buf, mod, id, n, p1, p2, p3, p4)
      n += 1
    end
    void_ref(b)
  end
#=
open(room_rad, "w") do out
  write_rad_material(out, radiance_light_material("bright", gray=100))
  write_rad_material(out, radiance_plastic_material("red_plastic", red=.7, green=.05, blue=.05, specularity=.5, roughness=.5))
  write_rad_material(out, radiance_plastic_material("gray_paint", gray=.5))
  write_rad_sphere(out, "bright", "fixture", xyz(2, 1, 1.5), .125)
  write_rad_sphere(out, "red_plastic", "ball", xyz(.7, 1.125, .625), .125)
  write_rad_box(out, "gray_paint", "1", xyz(0,0,0), 3, 2, 1.75)
  write_rad_material(out, radiance_plastic_material("blue_plastic", red=.1, green=.1, blue=.6, specularity=.05, roughness=.1))
  write_rad_quad(out, "blue_plastic", "box", "1540", xyz(0.98,0.88,0), xyz(0.98,0.88,0.5), xyz(0.5,0.75,0.5), xyz(0.5,0.75,0))
  write_rad_quad(out, "blue_plastic", "box", "4620", xyz(0.5,0.75,0.5), xyz(0.37,1.23,0.5), xyz(0.37,1.23,0), xyz(0.5,0.75,0))
  write_rad_quad(out, "blue_plastic", "box", "2310", xyz(0.37,1.23,0), xyz(0.85,1.36,0), xyz(0.98,0.88,0), xyz(0.5,0.75,0))
  write_rad_quad(out, "blue_plastic", "box", "3267", xyz(0.85,1.36,0), xyz(0.37,1.23,0), xyz(0.37,1.23,0.5), xyz(0.85,1.36,0.5))
  write_rad_quad(out, "blue_plastic", "box", "5137", xyz(0.98,0.88,0.5), xyz(0.98,0.88,0), xyz(0.85,1.36,0), xyz(0.85,1.36,0.5))
  write_rad_quad(out, "blue_plastic", "box", "6457", xyz(0.37,1.23,0.5), xyz(0.5,0.75,0.5), xyz(0.98,0.88,0.5), xyz(0.85,1.36,0.5))
  write_rad_material(out, radiance_metal_material("chrome", gray=.8, specularity=.9, roughness=0))
  write_rad_cylinder(out, "chrome", "fixture_support", xyz(2,1,1.5), .05, xyz(2,1,1.75))
end

  radiance_rview(radiance_oconv(room_rad), xyz(2.25, .375, 1), xyz(2.25, .375, 1) + vxyz(-.25, .125, -.125), 50, (.5, .5, .5))
=#

#=
realize(b::ACAD, s::Cone) =
  ACADCone(connection(b), add_z(s.cb, s.h), s.r, s.cb)
realize(b::ACAD, s::ConeFrustum) =
  ACADConeFrustum(connection(b), s.cb, s.rb, s.cb + vz(s.h, s.cb.cs), s.rt)
=#
realize(b::Radiance, s::Cylinder) =
  let bot_id = next_id(b),
      top_id = next_id(b),
      side_id = next_id(b),
      mod = get_material(b, s),
      kind = "cylinder",
      buf = buffer(b),
      bot = in_world(s.cb),
      top = in_world(s.cb + vz(s.h, s.cb.cs)),
      normal = unitized(top-bot)
    write_rad_ring(buf, #="mat$(kind)bot$(mod)"=# mod, bot_id, bot, 0, s.r, -normal)
    write_rad_cylinder(buf, #="mat$(kind)side$(mod)"=# mod, side_id, bot, s.r, top)
    write_rad_ring(buf, #="mat$(kind)top$(mod)"=# mod, top_id, top, 0, s.r, normal)
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
# In Radiance, PathSet need to be converted to a single path
realize_prism(b::Radiance, top, bot, side, path::PathSet, h::Real) =
  realize_prism(b, top, bot, side, convert(ClosedPath, path), h)

realize_pyramid_frustum(b::Radiance, top, bot, side, bot_vs::Locs, top_vs::Locs, closed=true) =
  let buf = buffer(b)
    if closed
      write_rad_polygon(buf, bot, next_id(b), reverse(bot_vs))
      write_rad_polygon(buf, top, next_id(b), top_vs)
    end
    for vs in zip(bot_vs, circshift(bot_vs, 1), circshift(top_vs, 1), top_vs)
      write_rad_polygon(buf, side, next_id(b), vs)
    end
  end

realize_polygon(b::Radiance, mat, vs::Locs, acw=true) =
  let buf = buffer(b)
#    polygon(vs, backend=autocad)
    write_rad_polygon(buf, mat, next_id(b), acw ? vs : reverse(vs))
  end

#=

Upon daysim processing, we need to compute the useful daylight illumination, which is the fraction of time that a sensor is within a given range of illumination

To do that, we start by reading the .ill file generated

ill = CSV.read(path, delim=" ", datarow=1)

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

ill = CSV.read(path, delim=" ", datarow=1)

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

#=
Radiance families need to know the different kinds of materials
that go on each surface.
In some cases it might be the same material, but in others, such
as slabs, outside walls, etc, we will have different materials.
=#

abstract type RadianceFamily <: Family end

struct RadianceMaterialFamily <: RadianceFamily
  material::RadianceMaterial
end

radiance_material_family(mat::RadianceMaterial) =
  RadianceMaterialFamily(mat)

struct RadianceSlabFamily <: RadianceFamily
  top_material::RadianceMaterial
  bottom_material::RadianceMaterial
  side_material::RadianceMaterial
end

radiance_slab_family(top::RadianceMaterial, bot::RadianceMaterial=top, side::RadianceMaterial=bot) =
  RadianceSlabFamily(top, bot, side)

radiance_roof_family = radiance_slab_family


struct RadianceWallFamily <: RadianceFamily
  right_material::RadianceMaterial
  left_material::RadianceMaterial
end

radiance_wall_family(right::RadianceMaterial, left::RadianceMaterial=right) =
  RadianceWallFamily(right, left)

backend_get_family_ref(b::Radiance, f::Family, rf::RadianceFamily) = rf

export radiance_material_family,
       radiance_slab_family,
       radiance_roof_family,
       radiance_wall_family,
       default_radiance_material


set_backend_family(default_wall_family(), radiance,
  radiance_wall_family(radiance_generic_interior_wall_70))
set_backend_family(default_slab_family(), radiance,
  radiance_slab_family(radiance_generic_floor_20, radiance_generic_ceiling_80))
set_backend_family(default_roof_family(), radiance,
  radiance_roof_family(radiance_generic_floor_20, radiance_outside_facade_30))
set_backend_family(default_beam_family(), radiance,
  radiance_material_family(radiance_generic_metal))
set_backend_family(default_column_family(), radiance,
  radiance_material_family(radiance_generic_metal))
set_backend_family(default_door_family(), radiance,
  radiance_material_family(radiance_generic_furniture_50))
set_backend_family(default_panel_family(), radiance,
  radiance_material_family(radiance_glass_material("GenericGlass", gray=0.3)))

# Radiance-specific operations

default_radiance_ground_material = Parameter(radiance_generic_floor_20)

#=
create_ground_plane(shapes, material=default_radiance_ground_material()) =
  if shapes == []
    error("No shapes selected for analysis. Use add-radiance-shape!.")
  else
    let (p0, p1) = bounding_box(union(shapes)),
        (center, ratio) = (quad_center(p0, p1, p2, p3),
                  distance(p0, p4)/distance(p0, p2));
     ratio == 0 ?
      error("Couldn"t compute height. Use add-radiance-shape!.") :
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
realize(b::Radiance, s::Beam) =
    ref(right_cuboid(s.p0, 0.2, 0.2, s.p1, 0))

=#

realize(b::Radiance, w::Window) = nothing
realize(b::Radiance, w::Door) = nothing

#=

Sensors are need on the surfaces that are intended for analysis.
They can be placed in many different ways but that are a few standard ones,
namely those that follow an regular distribution.
=#

analysis_nodes_height = Parameter(0.5)
analysis_nodes_separation_u = Parameter{Union{Missing, Real}}(missing)
analysis_nodes_separation_v = Parameter{Union{Missing, Real}}(missing)
analysis_nodes_separation = Parameter(4.0)

analysis_nodes_u_separation()::Real =
  ismissing(analysis_nodes_separation_u()) ?
    analysis_nodes_separation() :
    analysis_nodes_separation_u()

analysis_nodes_v_separation()::Real =
  ismissing(analysis_nodes_separation_v()) ?
    analysis_nodes_separation() :
    analysis_nodes_separation_v()

#=
sensors_from_surfaces(surfaces=radiance_surfaces(),
 height::Real=analysis_nodes_height(),
 sep_u::Real=analysis_nodes_u_separation(),
 sep_v::Real=analysis_nodes_v_separation()) =
  (map(surface_from -> let_values(nu(nv)(surface_nu_nv(surface_from, sep_u, sep_v))(), let nodes = map_inner_surface_division(pt -> pt && pt+vz(height), surface_from, nu, nv); filter(identity, append*(nodes)) end),
                surfaces))

=#

sensors_locations(s::Slab, b::Radiance) =
  let out_pts = path_vertices(s.countour),
      in_ptss = map(path_vertices, s.openings),
      p_ns = Main.plane_surface_sensor_locations(out_pts, in_ptss)
  end

struct RadianceVisualization <: LightingAnalysis
end

export radiance_visualization
radiance_visualization(b::Radiance=radiance; light=(1,1,1)) =
  let path=radiance_simulation_path(),
      radpath = export_geometry(b, path),
      matpath = export_materials(b, path),
      skypath = export_sky(b, path),
      octpath = radiance_oconv(radpath, matpath, skypath)
      @info radpath
      @info matpath
      @info skypath
    radiance_rview(octpath, b.camera, b.target, b.lens, light)
  end

export radiance_render
render_view(path::String, b::Radiance) =
  let radpath = export_geometry(b, path),
      matpath = export_materials(b, path),
      skypath = export_sky(b, path),
      octpath = radiance_oconv(radpath, matpath, skypath)
      @info radpath
      @info matpath
      @info skypath
    #radiance_rpict(octpath, b.camera, b.target, b.lens)
    diva_render(octpath, b.camera, b.target, b.lens)
  end

struct DaysimAnalysis <: LightingAnalysis
end

struct RadianceMapAnalysis <: LightingAnalysis
  sensors::Locs
  path::AbstractString
  location::AbstractString
  datetime::DateTime
  sky::AbstractString
  simulation_parameters::AbstractString
  results::Parameter{Any} # HACK: FIX THIS
end

export radiance_map_analysis
radiance_map_analysis(sensors;
                      location="Lisbon",
                      datetime=Dates.now(UTC),
                      sky=sky_for(location, datetime),
                      simulation_parameters="-ab 2 -ad 1000 -as 20 -ar 300 -aa 0.1",
                      path=radiance_simulation_path()) =
  RadianceMapAnalysis(sensors, path, location, datetime, sky, simulation_parameters, Parameter{Any}(nothing))

@defop analyze(a::Analysis)

analyze(a::RadianceMapAnalysis, b::Radiance) =
  let radpath = export_geometry(b, a.path),
      matpath = export_materials(b, a.path),
      skypath = export_sky(b, a.path),
      octpath = path_replace_suffix(a.path, ".oct"),
      ptspath = export_sensors(b, a.path, a.sensors),
      datpath = path_replace_suffix(a.path, ".dat")
    radiance_cmd("""oconv "$(matpath)" "$(skypath)" "$(radpath)" > "$(octpath)" """)
    radiance_cmd("""rtrace -I -h -dp 2048 -ms 0.063 -ds .2 -dt .05 -dc .75 -dr 3 -st .01 -lr 12 -lw .0005 $(a.simulation_parameters) "$(octpath)" < "$(ptspath)" > "$(datpath)" """)
    a.results(read_rtrace_results(a.path))
  end

export_geometry(b::Radiance, path::AbstractString) =
  let radpath = path_replace_suffix(path, ".rad"),
      buf = b.buffer()
    take!(buf)
    for s in b.shapes
      realize(b, s)
    end
    add_ground_plane(b)
    open(radpath, "w") do out
      write(out, String(take!(buf)))
    end
    radpath
  end

add_ground_plane(b::Radiance) =
  @warn "Not generating ground plane"

used_materials(b::Radiance) =
  let materials = unique(map(f -> realize(s.family, b), b.shapes))
      materials
    end

export_materials(b::Radiance, path::AbstractString) =
  let matpath = path_replace_suffix(path, "_materials.rad")
    open(matpath, "w") do out
      for (k,v) in b.materials
        write_rad_material(out, k)
      end
    end
    matpath
  end

export_sky(b::Radiance, path::AbstractString) =
  let skypath = path_replace_suffix(path, "_sky.rad")
    open(skypath, "w") do out
      write(out, b.sky)
    end
    skypath
  end

radiance_sensors = Parameter(Loc[])

export_sensors(path::Path, sensors=radiance_sensors()) =
  let ptspath = path_replace_suffix(path, ".pts")
    open(ptspath, "w") do out
      with(current_layer, sensor_layer()) do
        for p in sensors
          let (wp, wv) = (in_world(p), in_world(vz(1, p)))
            #;Just to see the sensor location
            point(p)
            #    (line p (+z p 3))
            write(out, "$(wp.x),$(wp.y),$(wp.z),$(wv.x),$(wv.y),$(wv.z)\n")
          end
        end
      end
    end
    ptspath
  end
