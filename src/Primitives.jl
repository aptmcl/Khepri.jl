decode_id(c::IO) =
  let id = decode_int(c)
    if id == -1
      error("Backend Error: $(decode_String(c))")
    else
      id
    end
  end

encode_BIMLevel = encode_int
decode_BIMLevel = decode_int_or_error
encode_FloorFamily = encode_int
decode_FloorFamily = decode_int_or_error

function request_operation(conn::IO, name)
  write(conn, Int32(0))
  encode_String(conn, name)
  op = read(conn, Int32)
  if op == -1
    error(name * " is not available")
  else
    op
  end
end

#=
create_op(name::ASCIIString, argtypes::Array{DataType}, rettype::DataType) = (
  op = request_operation(name);
  (args...) -> (conn = current_connection();
                write(conn, Int32(op));
                @assert length(args) == length(argtypes) "Incorrect number of args";
                for actual_arg in [argtype(arg) for (arg, argtype) in zip(args, argtypes)]
                  write(conn, actual_arg)
                end;
                read(conn, rettype)))

macro def_op(name, argtypes, rettype)
  :($(esc(name)) = create_op($(string(name)), $argtypes, $rettype))
end

@def_op(Circle, [XYZ, XYZ, Float64], Int32)
@def_op(Sphere, [XYZ, Float64], Int32)

circle = (op = request_operation("Circle");
          conn = current_connection();
          (c, n, r)-> (send_data(conn, Int8(op));
                       send_data(conn, c);
                       send_data(conn, n);
                       send_data(conn, Float64(r);
                       read(conn, Int32))))
=#

function parse_c_signature(sig)
  m = match(r"^ *(public|) *(\w+) *(\[\])? +(\w+) *\( *(.*) *\)", sig)
  ret = Symbol(m.captures[2])
  ret_is_array = m.captures[3]=="[]"
  name = Symbol(m.captures[4])
  params = split(m.captures[5], r" *, *", keepempty=false)
  function parse_c_decl(decl)
    m = match(r"^ *(\w+) *(\[\])? *(\w+)$", decl)
    (Symbol(m.captures[1]), m.captures[2]=="[]", Symbol(m.captures[3]))
  end
  (name, [parse_c_decl(decl) for decl in params], (ret, ret_is_array))
end

#=
const julia_type_for_c_type = Dict(
  :byte => :Int8,
  :double => :Float64,
  :float => :Float64,
  :int => :Int,
  :bool => :Bool,
  :Point3d => :XYZ,
  :Point2d => :XYZ,
  :Vector3d => :VXYZ,
  :string => :ASCIIString,
  :ObjectId => :Int32,
  :Entity => :Int32,
  :BIMLevel => :Int32,
  :FloorFamily => :Int32,
    ) #Either Int32 or Int64, depending on the architecture

julia_type(ctype, is_array) = is_array ? :(Vector{$(julia_type(ctype, false))}) : julia_type_for_c_type[ctype]
=#
export show_rpc, step_rpc
const show_rpc = Parameter(false)
const step_rpc = Parameter(false)
function initiate_rpc_call(conn, opcode, name)
    if step_rpc()
        print(stderr, "About to call $(name) [press ENTER]")
        readline()
    end
    if show_rpc()
        print(stderr, name)
    end
end
function complete_rpc_call(conn, opcode, result)
    if show_rpc()
        println(stderr, result == nothing ? "-> nothing" : "-> $(result)")
    end
    result
end

#=
This does not play well with reset_backend because the opcode are in global
scope and completely separated from the backend and, thus, the operation
still tries to use the same opcode even after a reset. Must think about
storing the generated functions on the backend itself, so that we can drop
it and restart.

Another option would be to have the table of generated methods separated
from the channel on the backend side, but I'm afraid that will prevent multiple
clients (not sure they are useful, though).
=#
function rpc(prefix, str)
    name, params, ret = parse_c_signature(str)
    func_name = Symbol(prefix, name)
    #Expr(:(=), Expr(:call, esc(name), [p[1] for p in params]...), Expr(:call, :+, params[1][1], 2))
    esc(quote
        $func_name =
          let opcode = -1
              buf = IOBuffer()
              (conn, $([:($(p[3])) for p in params]...)) -> begin
                opcode = Int32(request_operation(conn, $(string(name))))
                global $func_name = (conn, $([:($(p[3])) for p in params]...)) -> begin
                  initiate_rpc_call(conn, opcode, $(string(name)))
                  take!(buf) # Reset the buffer just in case there was an encoding error on a previous call
                  write(buf, opcode)
                  $([:($(Symbol("encode_", p[1], p[2] ? "_array" : ""))(buf, $(p[3]))) for p in params]...)
                  write(conn, take!(buf))
                  complete_rpc_call(conn, opcode, $(Symbol("decode_", ret[1], ret[2] ? "_array" : ""))(conn))
                end
                $func_name(conn, $([:($(p[3])) for p in params]...))
          end
        end
        #export $func_name
      end)
end
