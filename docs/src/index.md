```@meta
CurrentModule = Khepri
```

# Khepri

Khepri is the main entry point for the Khepri architectural design framework.
It re-exports [KhepriAutoCAD](https://aptmcl.github.io/KhepriAutoCAD.jl/stable),
providing access to all Khepri geometry, BIM, and backend operations with
AutoCAD as the default backend.

## Quick Start

```julia
using Khepri

# All KhepriBase + KhepriAutoCAD functionality is available
backend(autocad)
box(xyz(0, 0, 0), 5, 5, 5)
```

## Switching Backends

To use a different backend, replace `using Khepri` with the specific backend
package:

```julia
using KhepriRevit    # For Autodesk Revit
using KhepriBlender  # For Blender
using KhepriRhino    # For Rhinoceros
using KhepriUnity    # For Unity
```

## Package Index

### Core
- [KhepriBase](https://aptmcl.github.io/KhepriBase.jl/stable) — Geometry, coordinates, backend abstraction, shape proxies

### CAD/BIM Backends
- [KhepriAutoCAD](https://aptmcl.github.io/KhepriAutoCAD.jl/stable) — Autodesk AutoCAD
- [KhepriRevit](https://aptmcl.github.io/KhepriRevit.jl/stable) — Autodesk Revit
- [KhepriBlender](https://aptmcl.github.io/KhepriBlender.jl/stable) — Blender
- [KhepriRhino](https://aptmcl.github.io/KhepriRhino.jl/stable) — Rhinoceros
- [KhepriFreeCAD](https://aptmcl.github.io/KhepriFreeCAD.jl/stable) — FreeCAD
- [KhepriUnity](https://aptmcl.github.io/KhepriUnity.jl/stable) — Unity
- [KhepriUnreal](https://aptmcl.github.io/KhepriUnreal.jl/stable) — Unreal Engine

### Visualization Backends
- [KhepriMakie](https://aptmcl.github.io/KhepriMakie.jl/stable) — Makie.jl
- [KhepriThreejs](https://aptmcl.github.io/KhepriThreejs.jl/stable) — Three.js
- [KhepriMeshCat](https://aptmcl.github.io/KhepriMeshCat.jl/stable) — MeshCat
- [KhepriXeokit](https://aptmcl.github.io/KhepriXeokit.jl/stable) — xeokit

### Rendering
- [KhepriPOVRay](https://aptmcl.github.io/KhepriPOVRay.jl/stable) — POV-Ray
- [KhepriRadiance](https://aptmcl.github.io/KhepriRadiance.jl/stable) — Radiance

```@index
```

```@autodocs
Modules = [Khepri]
```
