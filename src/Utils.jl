export random, random_range, set_random_seed, RGB, rgb,
       division, map_division

#previous_random::Int = 12345
previous_random = 12345

set_random_seed(v::Int) =
  global previous_random = v

next_random(previous_random::Int) =
  let test = 16807*rem(previous_random,127773) - 2836*div(previous_random,127773)
    if test > 0
      if test > 2147483647
        test - 2147483647
      else
        test
      end
    else
      test + 2147483647
    end
  end

next_random!() =
  begin
    global previous_random = next_random(previous_random)
    previous_random
  end

random(x::Int) = rem(next_random!(), x)

random(x::Real) = x*next_random!()/2147483647.0

random_range(x0, x1) =
  if x0 == x1
    x0
  else
    x0 + random(x1 - x0)
  end

struct RGB
  r::Real
  g::Real
  b::Real
end

rgb(r::Real=0, g::Real=0, b::Real=0) =
    RGB(convert(UInt8, r), convert(UInt8, g), convert(UInt8, b))

required() = error("Required parameter")

division(t0, t1, n::Real, include_last::Bool=true) =
  let n = convert(Int, n), iter = range(t0, stop=t1, length=n + 1)
    collect(include_last ? iter : take(iter, n))
  end

# Generic objects are processed using map_division
division(obj::Any, n::Real) = map_division(identity, obj, n)


map_division(f, t0, t1, n::Real, include_last::Bool=true) =
  let n = convert(Int, n), iter = range(t0, stop=t1, length=n + 1)
    map(f, include_last ? iter : take(iter, n))
  end

map_division(f, u0, u1, nu::Real, include_last_u::Bool, v0, v1, nv::Real) =
  map_division(u -> map_division(v -> f(u, v), v0, v1, nv),
               u0, u1, nu, include_last_u)

map_division(f, u0, u1, nu::Real, v0, v1, nv::Real, include_last_v::Bool=true) =
  map_division(u -> map_division(v -> f(u, v), v0, v1, nv, include_last_v),
               u0, u1, nu)

map_division(f, u0, u1, nu::Real, include_last_u::Bool, v0, v1, nv::Real, include_last_v::Bool) =
  map_division(u -> map_division(v -> f(u, v), v0, v1, nv, include_last_v),
               u0, u1, nu, include_last_u)


# Renders and Films
export render_dir,
       render_user_dir,
       render_backend_dir,
       render_kind_dir,
       render_color_dir,
       render_ext,
       render_width,
       render_height,
       render_quality,
       render_exposure,
       set_render_dir,
       render_size,
       prepare_for_saving_file,
       render_pathname,
       render_view,
       rendering_with

# There is a render directory
const render_dir = Parameter(homedir())
# with a user-specific subdirectory
const render_user_dir = Parameter(".")
# with a backend-specific subdirectory
const render_backend_dir = Parameter(".")
# and with subdirectories for static images, movies, etc
const render_kind_dir = Parameter("Render")
# and with subdirectories for white, black, and colored renders
const render_color_dir = Parameter(".")
# containing files with different extensions
const render_ext = Parameter("png")

render_pathname(name::String) =
    realpath(
        joinpath(
            render_dir(),
            render_user_dir(),
            render_backend_dir(),
            render_kind_dir(),
            render_color_dir(),
            "$(name).$(render_ext())"))

const render_width = Parameter(1024)
const render_height = Parameter(768)
const render_quality = Parameter{Real}(0) # [-1, 1]
const render_exposure = Parameter{Real}(0)  # [-3, +3]
const render_floor_width = Parameter(1000)
const render_floor_height = Parameter(1000)

set_render_dir(val::String) = render_dir(realpath(val))

render_size(width::Integer, heigth::Integer) =
    (render_width(width), render_height(heigth))

prepare_for_saving_file(path::String) =
    let p = abspath(path)
        mkpath(dirname(path))
        rm(p, force=true)
        p
    end

export film_active, film_filename, film_frame, start_film, frame_filename, save_film_frame

const film_active = Parameter(false)
const film_filename = Parameter("")
const film_frame = Parameter(0)

start_film(name::String) =
    begin
        film_active(true)
        film_filename(name)
        film_frame(0)
    end

frame_filename(filename::String, i::Integer) =
    "$(filename)-frame-$(lpad(i,3,'0'))"

save_film_frame(obj::Any=true) =
    with(render_kind_dir, "Film") do
        render_view(frame_filename(film_filename(), film_frame()))
        film_frame(film_frame() + 1)
        obj
    end

rendering_with(f;
    dir=render_dir(),
    user_dir=render_user_dir(),
    backend_dir=render_backend_dir(),
    kind_dir=render_kind_dir(),
    color_dir=render_color_dir(),
    ext=render_ext(),
    width=render_width(),
    height=render_height(),
    quality=render_quality(),
    exposure=render_exposure(),
    floor_width=render_floor_width(),
    floor_height=render_floor_height()) =
    with(render_dir,dir) do
        with(render_user_dir, user_dir) do
            with(render_backend_dir, backend_dir) do
                with(render_kind_dir,kind_dir) do
                    with(render_color_dir,color_dir) do
                        with(render_ext,ext) do
                            with(render_width,width) do
                                with(render_height,height) do
                                    with(render_quality,quality) do
                                        with(render_exposure,exposure) do
                                           with(render_floor_width,floor_width) do
                                               with(render_floor_height,floor_height) do
                                                   f()
                                               end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
