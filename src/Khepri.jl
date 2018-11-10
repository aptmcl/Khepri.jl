module Khepri
using Base.Iterators
using Sockets
using PyCall

include("Parameters.jl")
include("Utils.jl")
include("Coords.jl")
include("Shapes.jl")
include("Geometry.jl")
include("Encoders.jl")
include("Primitives.jl")
include("Tikz.jl")
include("AutoCAD.jl")
include("Revit.jl")
#include("ArchiCAD.jl")
#include("Radiance.jl")
include("Unity.jl")
const com = PyNULL()

include("Robot.jl")
end
