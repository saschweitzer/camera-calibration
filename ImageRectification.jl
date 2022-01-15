###############################################################################
"""
Provides utilities to transform between physical (x,y,z) and pixel (u,v) coordinate systems.
"""
module ImageRectification
###############################################################################

###############################################################################
using StaticArrays              # provides SMatrix, SVector, for faster linear algebra
using ProgressMeter
###############################################################################

###############################################################################
"""
FindXYFromUVZ(u, v, z, DLT)

    Returns XY coordinates of a point.

    Calculated using the DLT formulation from eqn 3 in Holland et al. z is
    known a-priori (usually the water surface):
"""
function FindXYFromUVZ end
FindXYFromUVZ(u::Number, v::Number, z::Number, DLT::NamedTuple) = uvz2xy(u, v, z, DLT)
FindXYFromUVZ(uv::Array, z::Number, DLT::NamedTuple) = FindXYFromUVZ(view(uv,:,1), view(uv,:,2), z::Number, DLT)
FindXYFromUVZ(uv::Array, z::Array, DLT::NamedTuple) = FindXYFromUVZ(view(uv,:,1), view(uv,:,2), z::Array, DLT)
# 
function FindXYFromUVZ(u, v, z::Number, DLT::NamedTuple)
    x = Array{Float64}(undef, size(u))
    y = similar(x)
    for i in eachindex(u); x[i],y[i] = uvz2xy(u[i], v[i], z, DLT); end
    return x,y
end
function FindXYFromUVZ(u, v, z::Array, DLT::NamedTuple)
    x = Array{Float64}(undef, size(u))
    y = similar(x)
    for i in eachindex(u); x[i],y[i] = uvz2xy(u[i], v[i], z[i], DLT); end
    return x,y
end
# Helper functions for FindXYFromUVZ():
Amatrix(u, v, z, DLT::NamedTuple) = SMatrix{2,2}( u*DLT.L₉-DLT.L₁,  v*DLT.L₉-DLT.L₅, u*DLT.L₁₀-DLT.L₂, v*DLT.L₁₀-DLT.L₆ )
bvector(u, v, z, DLT::NamedTuple) = SVector( DLT.L₄-u + (DLT.L₃-u*DLT.L₁₁)*z, DLT.L₈-v + (DLT.L₇-v*DLT.L₁₁)*z )
uvz2xy(u, v, z, DLT) = Amatrix(u, v, z, DLT) \ bvector(u, v, z, DLT)
###############################################################################


###############################################################################    
""" Calculate the uv -> xyz coordinate conversion using each of the sets of DLTs: """
function CalculateModifiedxyCoordinates(u,v, WaterSurfaceElev, DLTs;
                                        wseNoiseFactor = 0.0)
    Perturb_wse = randn(size(DLTs,2)) .* wseNoiseFactor
    Modified_x, Modified_y = Rect.FindXYFromUVZ_DLT_z_sensitivity(u,v,
                                                                  # WaterSurfaceElev .+ randn(size(DLTs,2)) .* wseNoiseFactor,
                                                                  WaterSurfaceElev .+ Perturb_wse,
                                                                  DLTs)
    return Modified_x, Modified_y, Perturb_wse
end

###############################################################################    

###############################################################################
# This option is for sensitivity analysis of water surface elevation (and a single DLT set):
function FindXYFromUVZ_WSE_sensitivity(u::Array, v::Array, z::Array, DLT::NamedTuple)
    x = Array{Float64}(undef, size(u)..., length(z))
    y = Array{Float64}(undef, size(v)..., length(z))
    @showprogress 1 "Calculating $(length(u)) xy coordinates using $(size(DLT,2)) DLT sets:" for d in 1:length(z)
        for j in 1:size(u,2), i in 1:size(u,1)
            x[i,j,d],y[i,j,d] = FindXYFromUVZ(u[i,j], v[i,j], z[d], DLT)
        end
    end
    return x,y
end
# This option is for sensitivity analysis of multiple DLT sets:
function FindXYFromUVZ_DLT_sensitivity(u::Array, v::Array, z::Float64, DLTs::Array)
    x = Array{Float64}(undef, size(u)..., size(DLT,2))
    y = Array{Float64}(undef, size(v)..., size(DLT,2))
    @showprogress 1 "Calculating $(length(u)) xy coordinates using $(size(DLT,2)) DLT sets:" for d in 1:size(DLT,2)
        for j in 1:size(u,2), i in 1:size(u,1)
            x[i,j,d],y[i,j,d] = FindXYFromUVZ(u[i,j], v[i,j], z, DLT[d])
        end
    end
    return x,y
end
# This option is for sensitivity analysis of multiple DLT sets, along with
# multiple z values, for analysis of sensitivity to, e.g., water surface
# elevation:
function FindXYFromUVZ_DLT_z_sensitivity(u::Array, v::Array, z::Array, DLT::Array)
    (length(z) != size(DLT,2)) ? @error("lengths of z and DLT do not match!") : nothing
    x = Array{Float64}(undef, size(u)..., size(DLT,2))
    y = Array{Float64}(undef, size(v)..., size(DLT,2))
    @showprogress 1 "Calculating $(length(u)) xy coordinates using $(size(DLT,2)) DLT sets:" for d in 1:size(DLT,2) for j in 1:size(u,2), i in 1:size(u,1)
            x[i,j,d],y[i,j,d] = FindXYFromUVZ(u[i,j], v[i,j], z[d], DLT[d])
        end
    end
    return x,y
end
###############################################################################

###############################################################################
###############################################################################
###############################################################################
end                             # end of module
###############################################################################
###############################################################################
###############################################################################
