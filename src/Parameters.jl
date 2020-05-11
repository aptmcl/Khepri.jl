export with, Parameter, LazyParameter

mutable struct Parameter{T}
  value::T
end

# This no longer works in 1.0
#(p::Parameter{T}){T}()::T = p.value
#(p::Parameter{T}){T}(newvalue::T) = p.value = newvalue

(p::Parameter)() = p.value
(p::Parameter)(newvalue) = p.value = newvalue

#=
with{T}(f, p::Parameter{T}, newvalue::T) where T =
  begin
    oldvalue, p.value = p.value, newvalue
    try
      f()
    finally
      p.value = oldvalue
    end
  end
=#

# A more generic version (presumably, compatible with the previous one)
function with(f, p, newvalue)
  oldvalue = p()
  p(newvalue)
  try
    f()
  finally
    p(oldvalue)
  end
end

with(f, p, newvalue, others...) =
  with(p, newvalue) do
    with(f, others...)
  end

mutable struct LazyParameter{T}
  initializer::Union{DataType,Function} #This should be a more specific type: None->T
  value::Union{T, Nothing}
end

LazyParameter(T::DataType, initializer::Union{DataType,Function}) =
  LazyParameter{T}(initializer, nothing)

# This no longer works in 1.0
#(p::LazyParameter{T}){T}()::T = p.value == nothing ? (p.value = p.initializer()) : get(p.value)
#(p::LazyParameter{T}){T}(newvalue::T) = p.value = newvalue

(p::LazyParameter)() = isnothing(p.value) ? (p.value = p.initializer()) : p.value
(p::LazyParameter)(newvalue) = p.value = newvalue

import Base.reset
reset(p::LazyParameter{T}) where {T} = p.value = nothing
