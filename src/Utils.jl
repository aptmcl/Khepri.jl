# The List implementation is inspired by https://github.com/JuliaCollections/DataStructures.jl
export List, Nil, Cons, list, cons, nil, head, tail
abstract type List{T} end
struct Nil{T} <: List{T} end
struct Cons{T} <: List{T}
    head::T
    tail::List{T}
end

nil(T) = Nil{T}()
nil() = nil(Any)
cons(h, t::List{T}) where {T} = Cons{T}(h, t)
list() = nil()
Base.isempty(lst::Nil) = true
Base.isempty(lst::Cons) = false
Base.getindex(lst::Nil, i) = throw(BoundsError(lst, i))
Base.getindex(lst::Cons, i) = i == 1 ? lst.head : getindex(lst.tail, i-1)

list(elts...) =
    let l = nil()
        for i = length(elts):-1:1
            l = cons(elts[i], l)
        end
        l
    end

list(elts::T...) where {T} =
    let l = nil(T)
        for i = length(elts):-1:1
            l = cons(elts[i], l)
        end
        l
    end
head(x::Cons) = x.head
tail(x::Cons) = x.tail

import Base.==
==(x::Nil, y::Nil) = true
==(x::Cons, y::Cons) = (x.head == y.head) && (x.tail == y.tail)

Base.show(io::IO, lst::Nil) = print(io, "list()")
Base.show(io::IO, lst::Cons) =
    begin
        print(io, "list(")
        show(io, head(lst))
        for e in tail(lst)
            print(io, ", ")
            show(io, e)
        end
        print(io, ")")
    end

Base.length(l::Nil) = 0
Base.length(l::Cons) =
    let n = 0
        for i in l
            n += 1
        end
        n
    end

Base.map(f::Base.Callable, lst::List) = list((f(e) for e in lst)...)
Base.filter(f::Function, lst::List) = list((e for e in lst if f(e))...)

Base.cat() = list()
Base.cat(lst::List, lsts::List...) =
    let T = typeof(lst).parameters[1]
        n = length(lst)
        for l in lsts
            T2 = typeof(l).parameters[1]
            T = typejoin(T, T2)
            n += length(l)
        end
        elems = Vector{T}(undef, n)
        i = 1
        for e in lst
            elems[i] = e
            i += 1
        end
        for lst in lsts
            for e in lst
                elems[i] = e
                i += 1
            end
        end
        let l = nil(T)
            for i = i-1:-1:1
                l = cons(elems[i], l)
            end
            l
        end
    end

Base.iterate(l::List, ::Nil) = nothing
Base.iterate(l::List, state::Cons = l) = state.head, state.tail

# Lists can be converted to Arrays

Base.convert(::Type{Array{S,1}}, l::List{T}) where {S, T <: S} = collect(T, l)





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
