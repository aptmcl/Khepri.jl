module Khepri
using Base.Iterators
using Sockets
using Dates
using ColorTypes
using Colors
using StaticArrays
using LinearAlgebra
using IntervalSets
using Dierckx
#using PyCall
using MeshCat
#using PlotlyJS

import Base.show, Base.zero, Base.iterate
import Base.convert
import Base.getindex, Base.firstindex, Base.lastindex, Base.broadcastable
import Base.union

include("Parameters.jl")
include("Utils.jl")
include("Coords.jl")
include("Paths.jl")
include("Shapes.jl")
include("BIM.jl")
include("Geometry.jl")
include("Encoders.jl")
include("Primitives.jl")
#include("Tikz.jl")
include("AutoCAD.jl")
include("Rhinoceros.jl")
#include("Revit.jl")
#include("ArchiCAD.jl")
include("Radiance.jl")
include("POVRay.jl")
#include("Notebook.jl")
#include("MeshCat.jl")
include("Grasshopper.jl")
include("Unity.jl")
#include("Unreal.jl")
include("Camera.jl")
#include("Frame3DD.jl")
#const com = PyNULL()
#include("Robot.jl")
include("BIMFamilies.jl")
end
