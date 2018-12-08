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
    render_view(prepare_for_saving_file(frame_filename(film_filename(), film_frame())))
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
    with(render_dir, dir) do
        with(render_user_dir, user_dir) do
            with(render_backend_dir, backend_dir) do
                with(render_kind_dir, kind_dir) do
                    with(render_color_dir, color_dir) do
                        with(render_ext, ext) do
                            with(render_width, width) do
                                with(render_height, height) do
                                    with(render_quality, quality) do
                                        with(render_exposure, exposure) do
                                           with(render_floor_width, floor_width) do
                                               with(render_floor_height, floor_height) do
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

export default_lens,
       track_still_target,
       track_moving_target,
       walkthrough,
       panning,
       lens_zoom,
       dolly_effect_back,
       dolly_effect_forth


default_lens = Parameter(50)

set_view_save_frame(camera, target, lens=default_lens()) =
  begin
    view(camera, target, lens)
    save_film_frame()
  end

# Tracking

# Rotation around target

# target object remains still and the camera moves around it
# parameters: camera path (list of locations) + fixed target (1 location)
# number of locations in the path defines the number of frames

track_still_target(camera_path, target, lens=default_lens()) =
  for camera in camera_path
    set_view_save_frame(camera, target, lens)
  end

# Camera paths
rotation_path(center, r, start_phi, stop_phi, h, frames) =
  map(fi -> center+vcyl(r, fi, h),
      division(start_phi, stop_phi, frames))

spiral_path(center, r, r_delta, start_phi, stop_phi, h, h_delta, frames) =
  map((fi, r, h) -> center+vcyl(r, fi, h),
      division(start_phi, stop_phi, frames),
      division(r, r*r_delta, frames),
      division(h, h*h_delta, frames))

# Moving camera and target

# target object moves and the camera makes a similar movement
# parameters: camera (initial location) + target path (list of locations)
# number of locations in the path defines the number of frames

set_view_save_frames(camera_path, target_path, lens) =
  for (camera, target) in zip(camera_path, target_path)
    set_view_save_frame(camera, target, lens)
  end

twin_path(path, v) = map(p -> p+v, path)

track_moving_target(camera, target_path, lens=default_lens()) =
  let v = camera - target_path[1]
      camera_path = twin_path(target_path, v)
    set_view_save_frames(camera_path, target_path, lens)
  end

linear_path(p1, p2, frames) =
  let v = p2 - p1
    map_division(t -> p1 + v*t, 0, 1, frames)
  end

# Walkthroughs

# targets object moves and the camera follows a few positions behind
# parameters: path (list of locations) + number of distance locations

walkthrough(path, camera_spread, lens=default_lens()) =
  let target_path = path[camera_spread:end]
      camera_path = path[1:end-camera_spread]
    set_view_save_frames(camera_path, target_path, lens)
  end

# Panning

# Moving target, still camera

# Camera rotation
# parameters: fixed camera + target path
# number of locations in the path defines the number of frames

panning(camera, target_path, lens=default_lens()) =
  for target in target_path
    set_view_save_frame(camera, target, lens)
  end

#;(define my-path (rotation-path (x 0) -50 0 2pi 0 100)) ;example
#;(panning (x 0) my-path) ;example

# Zooming

# Lens zoom

# Changing lens focal length while camera position remains the same
# paremeters: fixed camera + fixed target + zomming scale (lens degrees to change) + number of frames
# zomming scale: delta > 1 => zooms in, delta < 1 => zooms out - |delta| >> zoom mais violeto

lens_zoom(camera, target, delta, frames, lens=default_lens()) =
  let lenses = range(lens, step=delta, length=frames)
    for lens in lenses
      set_view_save_frame(camera, target, lens)
    end
  end

# Dolly zoom

# Moving camera backwards as the lens increases the focal length
# parameters: camera start + fixed target + zomming scale (lens degrees to change) + number of frames

#=
dolly_effect(cur_camera, target, cur_lens, new_camera) =
  let cur_dist = distance(cur_camera, target)
      new_dist = distance(new_camera, target)
      new_lens = (cur_lens*new_dist)/cur_dist
    set_view_save_frame(new_camera, target, new_lens)
  end

dolly_effect_pull_push(back_forth, delta) =
  camera, target, lens = view()
  let d = distance(camera, target)
    let (new_camera_back, new_camera_forth) = (target+(camera-target)*(d+delta)/d, target+(camera-target)*(d-delta)/d)
      back/forth == 0 ? dolly_effect(camera, target, lens, new_camera_back) : dolly_effect(camera, target, lens, new_camera_forth)
    end
  end

dolly_zoom(back/forth, camera_start, target, lens, delta, frames) =
  begin
    set_view_save_frame(camera_start, target, lens)
    for i in 0:frames-1
      dolly_effect_pull/push(back/forth, delta)
    end
  end

current_view_dolly_zoom(back/forth, delta, frames) =
  for i in 0:frames-1
    dolly_effect_pull/push(back/forth, delta)
  end

=#
#;(view (xyz 5 -0.8 0.5) (xyz 5 0 0.5) 10)
#;(dolly-zoom (xyz 5 -0.8 0.5) (xyz 5 0 0.5) 0.05 100) ;example

# When we don't have access to the parameters of the current view
# parameters: zomming scale + camera start + fixed target + initial lens value + number of frames

dolly_effect_back(delta, camera, target, lens, frames) =
  let d = distance(camera, target)
      new_camera = target+(camera-target)*(d+delta)/d
      new_d = distance(new_camera, target)
      new_lens = (lens*new_d)/d
    if frames <= 0
      nothing
    else
      set_view_save_frame(camera, target, lens)
      dolly_effect_back(delta, new_camera, target, new_lens, frames-1)
    end
  end

dolly_effect_forth(delta, camera, target, lens, frames) =
  let d = distance(camera, target)
      new_camera = target+(camera-target)*(d-delta)/d
      new_d = distance(new_camera, target)
      new_lens = (lens*new_d)/d
    if frames <= 0
      nothing
    else
      set_view_save_frame(camera, target, lens)
      dolly_effect_forth(delta, new_camera, target, new_lens, frames-1)
    end
  end
