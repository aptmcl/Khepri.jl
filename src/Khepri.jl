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
include("ArchiCAD.jl")
#include("Radiance.jl")

const com = PyNULL()

function __init__()
    copy!(com, pyimport("win32com.client"))
end

include("Robot.jl")
end
