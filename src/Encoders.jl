export encode_String, decode_String,
       encode_string, decode_string,
       encode_string_array,
       Guid, Guids,
       encode_Guid, decode_Guid,
       encode_Guid_array, decode_Guid_array,
       is_empty_guid,
       encode_bool, decode_bool,
       encode_byte, decode_byte,
       encode_int, decode_int, decode_int_or_error, decode_int_or_nothing,
       encode_long, decode_long, decode_long_or_error,
       encode_int_array, decode_int_array, decode_int_or_error_array,
       encode_float, decode_float,
       encode_float_array, decode_float_array,
       encode_double, decode_double,
       encode_double_array, decode_double_array,
       encode_Point3d, decode_Point3d,
       encode_Point3d_array, decode_Point3d_array,
       encode_Point2d, decode_Point2d,
       encode_Point2d_array, decode_Point2d_array,
       encode_Vector3d, decode_Vector3d,
       encode_Frame3d, decode_Frame3d,
       decode_void,
       encode_object, decode_object,
       encode_object_array, decode_object_array

struct BackendError
  msg::String
  backtrace
end

show(io::IO, e::BackendError) =
  print(io, "Backend Error: $(e.msg)")

backend_error(c::IO) =
  throw(BackendError(decode_String(c), backtrace()))

## To present errors in the backends that call back to Julia
exception_backtrace(e) = backtrace()
exception_backtrace(e::BackendError) = e.backtrace

export errormsg
errormsg(e) =
  sprint((io, e) -> showerror(io, e, exception_backtrace(e), backtrace=true), e)

encode_String(c::IO, v::AbstractString) = begin
  str = string(v)
  size = length(str)
  array = UInt8[]
  while true
    byte = size & 0x7f
    size >>= 7
    if size > 0
      push!(array, byte | 0x80)
    else
      push!(array, byte)
      break
    end
  end
  write(c, array)
  write(c, str)
end

decode_String(c::IO) = begin
  loop(size::Int, shift::Int) = begin
    b = convert(Int, read(c, UInt8))
    size = size | ((b & 0x7f) << shift)
    if (b & 0x80) == 0
      String(read(c, size))
    else
      loop(size, shift + 7)
    end
  end
  loop(0, 0)
end

# C# uses two different names for strings
encode_string = encode_String
decode_string = decode_String


Guid = Vector{UInt8}
Guids = Vector{Guid}

is_empty_guid(g::Guid) = all(v -> v == 0, g)

encode_Guid(c::IO, v::Guid) = write(c, v)
decode_Guid(c::IO) =
  let guid = read(c, 16)
      if is_empty_guid(guid)
        backend_error(c)
      else
        guid
      end
  end


encode_Guid_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_Guid(c, e) end
end
decode_Guid_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{Guid}(undef, len)
  for i in 1:len
    r[i] = decode_Guid(c)
  end
  r
end



encode_string = encode_String
decode_string = decode_String

encode_bool(c::IO, v::Bool) = encode_byte(c, v ? UInt8(1) : UInt8(0))
decode_bool(c::IO) =
  let i = decode_byte(c)
    if i == 127
      backend_error(c)
    else
      i == 1
    end
  end

encode_byte(c::IO, v::UInt8) = write(c, convert(UInt8, v))
decode_byte(c::IO) = convert(UInt8, read(c, UInt8))

encode_int(c::IO, v::Int) = write(c, convert(Int32, v))
decode_int(c::IO) = convert(Int, read(c, Int32))

encode_long(c::IO, v::Int) = write(c, convert(Int64, v))
decode_long(c::IO) = convert(Int, read(c, Int64))

encode_double(c::IO, v::Real) = write(c, convert(Float64, v))
decode_double(c::IO) =
    let d = read(c, Float64)
        if isnan(d)
            backend_error(c)
        else
            d
        end
    end
encode_double3(c::IO, v0::Real, v1::Real, v2::Real) = begin
    encode_double(c, v0)
    encode_double(c, v1)
    encode_double(c, v2)
end

encode_float(c::IO, v::Real) = write(c, convert(Float32, v))
decode_float(c::IO) =
    let d = read(c, Float32)
        if isnan(d)
            backend_error(c)
        else
            convert(Float64, d)
        end
    end
encode_float3(c::IO, v0::Real, v1::Real, v2::Real) = begin
    encode_float(c, v0)
    encode_float(c, v1)
    encode_float(c, v2)
end

decode_int_or_error_numbered(err_num) = (c::IO) ->
  let i = decode_int(c)
    if i == err_num
      backend_error(c)
    else
      i
    end
  end

decode_int_or_error = decode_int_or_error_numbered(-1)

decode_int_or_nothing(c::IO) =
  let i = decode_int_or_error(c)
    if i == -2
      nothing
    else
      i
    end
  end

decode_long_or_error_numbered(err_num) = (c::IO) ->
  let i = decode_long(c)
    if i == err_num
      backend_error(c)
    else
      i
    end
  end

decode_long_or_error = decode_long_or_error_numbered(-1)

encode_string_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_string(c, e) end
end

encode_double_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_double(c, e) end
end
decode_double_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{Float64}(undef, len)
  for i in 1:len
    r[i] = decode_double(c)
  end
  r
end

encode_int_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_int(c, e) end
end
decode_int_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{Int}(undef, len)
  for i in 1:len
    r[i] = decode_int(c)
  end
  r
end

decode_int_or_error_numbered_array(err_num) =
    let decode_int_or_err_num = decode_int_or_error_numbered(err_num)
        (c::IO) -> begin
            len = decode_int_or_error(c)
            r = Vector{Int}(undef, len)
            for i in 1:len
                r[i] = decode_int_or_err_num(c)
            end
            r
        end
    end

decode_int_or_error_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{Int}(undef, len)
  for i in 1:len
    r[i] = decode_int_or_error(c)
  end
  r
end

# Must convert from local to world coordinates
encode_Point3d(c::IO, v::XYZ) = begin
  v = in_world(v)
  encode_double3(c, v.x, v.y, v.z)
end
decode_Point3d(c::IO) = xyz(decode_double(c), decode_double(c), decode_double(c), world_cs)
encode_Point3d_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_Point3d(c, e) end
end
decode_Point3d_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{XYZ}(undef, len)
  for i in 1:len
    r[i] = decode_Point3d(c)
  end
  r
end
encode_Point2d(c::IO, v::XYZ) = begin
  v = in_world(v)
  @assert v.z == zero(v.z) "Non zero Z coordinate";
  encode_double(c, v.x)
  encode_double(c, v.y)
end
decode_Point2d(c::IO) = xyz(decode_double(c), decode_double(c), 0.0, world_cs)
encode_Point2d_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v
    encode_Point2d(c, e)
  end
end
decode_Point2d_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{XYZ}(undef, len)
  for i in 1:len
    r[i] = decode_Point2d(c)
  end
  r
end

encode_Vector3d(c::IO, v::VXYZ) = begin
  v = in_world(v)
  encode_double3(c, v.x, v.y, v.z)
end
decode_Vector3d(c::IO) = vxyz(decode_double(c), decode_double(c), decode_double(c), world_cs)

encode_Vector3d_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_Vector3d(c, e) end
end
decode_Vector3d_array(c::IO) = begin
  len = decode_int_or_error(c)
  r = Vector{VXYZ}(undef, len)
  for i in 1:len
    r[i] = decode_Vector3d(c)
  end
  r
end

encode_Frame3d(c::IO, v::XYZ) = begin
  encode_Point3d(c, v)
  t = v.cs.transform
  encode_double3(c, t[1,1], t[2,1], t[3,1])
  encode_double3(c, t[1,2], t[2,2], t[3,2])
  encode_double3(c, t[1,3], t[2,3], t[3,3])
end

decode_Frame3d(c::IO) =
    u0(cs_from_o_vx_vy_vz(decode_Point3d(c), decode_Vector3d(c), decode_Vector3d(c), decode_Vector3d(c)))

decode_void(c::IO) = begin
  v = decode_byte(c)
  if v == 127
    backend_error(c)
  end
end

encode_object(c::IO, v::Bool) = (encode_byte(c, 0x0); encode_bool(c, v))
encode_object(c::IO, v::UInt8) = (encode_byte(c, 0x1); encode_byte(c, v))
encode_object(c::IO, v::Int32) = (encode_byte(c, 0x2); encode_int(c, v))
encode_object(c::IO, v::Int64) = (encode_byte(c, 0x3); encode_long(c, v))
encode_object(c::IO, v::Float32) = (encode_byte(c, 0x4); encode_float(c, v))
encode_object(c::IO, v::Float64) = (encode_byte(c, 0x5); encode_double(c, v))

encode_object_array(c::IO, v::Vector) = begin
  encode_int(c, length(v))
  for e in v encode_object(c, e) end
end

create_backend_connection(backend::AbstractString, port::Integer) =
  for i in 1:10
    try
      return connect(port)
    catch e
      @info("Please, start/restart $(backend).")
      sleep(5)
      if i == 9
        throw(e)
      end
    end
  end
