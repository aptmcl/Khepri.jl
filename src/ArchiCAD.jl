export
    archicad

#
abstract type ARCHKey end
const ARCHId = Int
const ARCHIds = Vector{ARCHId}
const ARCHRef = GenericRef{ARCHKey, ARCHId}
const ARCHRefs = Vector{ARCHRef}
const ARCHNativeRef = NativeRef{ARCHKey, ARCHId}
const ARCH = Socket_Backend{ARCHKey, ARCHId}

void_ref(b::ARCH) = ARCHNativeRef(-1)

create_archicad_connection() = create_backend_connection("ArchiCAD", 11002)

const archicad = ARCH(LazyParameter(TCPSocket, create_archicad_connection))
#current_backend(archicad)
