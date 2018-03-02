module DiskTracer

# package code goes here
export constants,
    model

# These statements just straight up dump the source code directly here, making DiskJockey.jl
# act as one giant file.
include("constants.jl")
include("model.jl")
include("geometry.jl")
include("trace.jl")

end # module
