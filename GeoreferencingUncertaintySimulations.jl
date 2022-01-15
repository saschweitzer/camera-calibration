###############################################################################
"""

# General

This module contains Julia code for performing uncertainty analysis of camera
calbiration and georeferencing, as described in Schweitzer & Cowen, 2021.


# References
Schweitzer, S. A., & Cowen, E. A. (2021). "Instantaneous river-wide water surface
velocity field measurements at centimeter scales using infrared quantitative
image velocimetry". Water Resources Research, 57,
e2020WR029279. https://doi.org/10.1029/2020WR029279

# License
This software is released under the MIT license, however, the authors request to
be credited if you use the software, preferably by citing the above-mentioned
Schweitzer & Cowen (2021) paper.

MIT License

Copyright (c) 2022 Seth A. Schweitzer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
module GeoreferencingUncertaintySimulations
###############################################################################

###############################################################################
using Statistics
using ProgressMeter
# 
import ImageRectification; Rect = ImageRectification
import ExtrinsicCameraCalibration; Calib = ExtrinsicCameraCalibration
import RemoveImageDistortion; RadialDist = RemoveImageDistortion
###############################################################################

###############################################################################
export SimulateDLTFromGCPsWithAddedNoise
###############################################################################


###############################################################################
function GeorefGCPSimulationWrapper(TrueGCPxyz, TrueGCPuv, TrueWaterSurfaceElev, u, v, χinitial, ι;
                                     xyzNoiseFactor = 0.0,
                                     uvNoiseFactor = 0.0,
                                     wseNoiseFactor = 0.0,
                                     AddRadialDistortion = false,
                                     nSimulations = 500,
                                     MaxIterations = 50)
    InfoString = "n Simulations = $nSimulations\nAdding the following types of noise:"
    InfoString = (xyzNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nXYZ $xyzNoiseFactor"])
    InfoString = (uvNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nUV $uvNoiseFactor"])
    InfoString = (wseNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nWSE $wseNoiseFactor"])
    InfoString = join([InfoString, "\nMaxIterations $MaxIterations"])
    @info(InfoString)
    # 
    # Calculate the calibration without adding any noise:
    χcalculated = Calib.HollandEtAlExtrinsicCalibration(TrueGCPuv, TrueGCPxyz, χinitial, ι,
                                                            MaxIterations=MaxIterations)
    # Calculate the DLTs without adding any noise:
    DLTcalculated = Calib.CalculateDLTcoefficients(χcalculated, ι)
    # Calculate the x,y coordinates without adding any noise:
    Calculated_x, Calculated_y = Rect.FindXYFromUVZ(u,v, TrueWaterSurfaceElev, DLTcalculated)
    # Generate a set of DLTs that by running the calibraiotn with noise added to
    # the GCP uv and/or xyz coordinates:
    DLTsimulated, GCP_perturb_uv_xyz = SimulateDLTFromGCPsWithAddedNoise(TrueGCPxyz, TrueGCPuv, u, v, χinitial, ι;
                                                                         xyzNoiseFactor = xyzNoiseFactor,
                                                                         uvNoiseFactor = uvNoiseFactor,
                                                                         AddRadialDistortion = AddRadialDistortion,
                                                                         nSimulations = nSimulations,
                                                                         MaxIterations = MaxIterations)
    # Calculate the x,y coordinates of pixels that would be found using the
    # simulated DLTs (i.e., with noise added to the GCP coordinates). Optionally
    # add more noise related to uncertainty in determining the water surface
    # elevation at this step:
    Modified_x, Modified_y, Perturb_wse = CalculateModifiedxyCoordinates(u,v,
                                                                         TrueWaterSurfaceElev, DLTsimulated;
                                                                         wseNoiseFactor=wseNoiseFactor)
    Perturbations = (GCPxyz=GCP_perturb_uv_xyz[1], GCPuv=GCP_perturb_uv_xyz[2], wse=Perturb_wse)
    return (;Calculated_x, Calculated_y, Modified_x, Modified_y, InfoString, Perturbations)
end
###############################################################################




###############################################################################
"""
    function SimulateDLTFromGCPsWithAddedNoise(TrueGCPxyz, TrueGCPuv, u, v, z, χinitial, ι;
                                       xyzNoiseFactor = 1.0, uvNoiseFactor = 1.0,
                                       AddRadialDistortion = false,
                                       nSimulations = 500,
            GeorefSim.GeorefGCPSimulationWrapper2                           MaxIterations = 50)

    Calculates extrinsic calibration for a set of inputs with random noise.

    The function adds normally distributed noise with mean zero and input-dependent
    standard deviation to the pixel and/or physical coordinates of the GCPs.
    """
function SimulateDLTFromGCPsWithAddedNoise(TrueGCPxyz, TrueGCPuv, u, v, χinitial, ι;
                                   xyzNoiseFactor = 1.0, uvNoiseFactor = 1.0,
                                   AddRadialDistortion = false,
                                   nSimulations = 500,
                                   MaxIterations = 50)
    xyzNoiseFactor = isa(xyzNoiseFactor, Number) ? repeat([xyzNoiseFactor], 3) : xyzNoiseFactor
    uvNoiseFactor = isa(uvNoiseFactor, Number) ? repeat([uvNoiseFactor], 2) : uvNoiseFactor
    # DLTs = fill(NaN, 11, nSimulations)
    DLTs = Array{Any}(undef, nSimulations)
    GCP_xyz_perturb = [ hcat(randn(size(TrueGCPxyz,1)) .* xyzNoiseFactor[1],
                             randn(size(TrueGCPxyz,1)) .* xyzNoiseFactor[2],
                             randn(size(TrueGCPxyz,1)) .* xyzNoiseFactor[3]) for i in 1:nSimulations ]
    GCP_uv_perturb = [ hcat(randn(size(TrueGCPuv,1)) .* uvNoiseFactor[1],
                            randn(size(TrueGCPuv,1)) .* uvNoiseFactor[2]) for i in 1:nSimulations ]
    InputGCPuv = ( AddRadialDistortion ?
        hcat(RadialDist.SimulateRadialDistortion(TrueGCPuv[:,1], TrueGCPuv[:,2], ι)...) :
        TrueGCPuv )
    @showprogress 1 "simulating DLT uncertainty..." for i in 1:nSimulations
        χmodified = Calib.HollandEtAlExtrinsicCalibration(
            InputGCPuv .+ GCP_uv_perturb[i],
            TrueGCPxyz .+ GCP_xyz_perturb[i],
            χinitial,
            ι;
            MaxIterations=MaxIterations,
            PrintoutStats=false)
        DLTs[i] = Calib.CalculateDLTcoefficients(χmodified, ι)
    end
    return DLTs, (GCP_xyz_perturb=GCP_xyz_perturb, GCP_uv_perturb=GCP_uv_perturb)
end
###############################################################################    

###############################################################################    
""" Calculate the uv -> xyz coordinate conversion using each of the sets of DLTs """
function CalculateModifiedxyCoordinates(u,v, WaterSurfaceElev, DLTs;
                                        wseNoiseFactor = 0.0)
    Perturb_wse = randn(size(DLTs,2)) .* wseNoiseFactor
    Modified_x, Modified_y = Rect.FindXYFromUVZ_DLT_z_sensitivity(u,v,
                                                                  WaterSurfaceElev .+ Perturb_wse,
                                                                  DLTs)
    return Modified_x, Modified_y, Perturb_wse
end
###############################################################################    


###############################################################################    
"""
Calculate the distance b/w the simulated and "true" xy coordinates of each pixel
footprint:
"""
function CalculateDistanceFromOriginal(Modified_x, Calculated_x, Modified_y, Calculated_y)
    return hypot.((Modified_x .- Calculated_x), (Modified_y .- Calculated_y))
end
###############################################################################    


###############################################################################
function GeorefExtrinsicParamsSimulationWrapper(χinitial, ι,
                                                TrueWaterSurfaceElev,
                                                u, v;
                                                φNoiseFactor = 0.0,
                                                τNoiseFactor = 0.0,
                                                σNoiseFactor = 0.0,
                                                x_cNoiseFactor = 0.03,
                                                y_cNoiseFactor = 0.03,
                                                z_cNoiseFactor = 0.05,
                                                wseNoiseFactor = 0.03,
                                                AddRadialDistortion = false,
                                                nSimulations = 500,
                                                MaxIterations = 50)
    InfoString = "n Simulations = $nSimulations\nAdding the following types of noise:"
    roundedrad2deg(t) = round(rad2deg(t), digits=3)
    InfoString = (φNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nφ $(roundedrad2deg(φNoiseFactor)) [deg]"])
    InfoString = (τNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nτ $(roundedrad2deg(τNoiseFactor)) [deg]"])
    InfoString = (σNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nσ $(roundedrad2deg(σNoiseFactor)) [deg]"])
    InfoString = (z_cNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nz_c $z_cNoiseFactor"])
    InfoString = (wseNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nWSE $wseNoiseFactor"])
    InfoString = join([InfoString, "\nMaxIterations $MaxIterations"])
    @info(InfoString)
    # 
    DLTcalculated = Calib.CalculateDLTcoefficients(χinitial, ι)
    # Calculate the x,y coordinates without adding any noise:
    Calculated_x, Calculated_y = Rect.FindXYFromUVZ(u,v, TrueWaterSurfaceElev, DLTcalculated)
    # Generate a set of DLTs that by running the calibraiotn with noise added to
    # the GCP uv and/or xyz coordinates:
    DLTsimulated, ExtrinsicParamsPerturbs = SimulateDLTFromExtrinsicParamsWithAddedNoise(χinitial, ι, u, v;
                                                                                         φNoiseFactor = φNoiseFactor,
                                                                                         τNoiseFactor = τNoiseFactor,
                                                                                         σNoiseFactor = σNoiseFactor,
                                                                                         z_cNoiseFactor = z_cNoiseFactor,
                                                                                         nSimulations = nSimulations,
                                                                                         MaxIterations = MaxIterations)
    # Simulate radial distortion by changing the pixel coordinates of the input grid:
    u_sim, v_sim = AddRadialDistortion ? RadialDist.SimulateRadialDistortion(u, v, ι) : (u, v)
    # Calculate the x,y coordinates of pixels that would be found using the
    # simulated DLTs (i.e., with noise added to the extrinsic
    # calibration). Optionally add more noise related to uncertainty in
    # determining the water surface elevation at this step:
    Modified_x, Modified_y, Perturb_wse = CalculateModifiedxyCoordinates(u_sim, v_sim,
                                                                         TrueWaterSurfaceElev,
                                                                         DLTsimulated;
                                                                         wseNoiseFactor=wseNoiseFactor)
    # Calculate the horizontal difference in position between the set of xy
    # coordinates calculated with random noise added, and the set of coordinates
    # found from the original input (no additional noise added to GCP or water
    # surface elevation coordinates):
    DeltaDist = CalculateDistanceFromOriginal(Modified_x, Calculated_x, Modified_y, Calculated_y)
    # return DeltaDist, InfoString, (; ExtrinsicParamsPerturbs..., Perturb_wse)
    return (; Calculated_x, Calculated_y, Modified_x, Modified_y, InfoString, ExtrinsicParamsPerturbs..., Perturb_wse)
end
###############################################################################


###############################################################################    
function SimulateDLTFromExtrinsicParamsWithAddedNoise(
    χinitial, ι, u, v;
    φNoiseFactor = 0.0,
    τNoiseFactor = 0.0,
    σNoiseFactor = 0.0,
    x_cNoiseFactor = 0.05,
    y_cNoiseFactor = 0.05,
    z_cNoiseFactor = 0.1,                 
    nSimulations = 500,
    MaxIterations = 50)
    # 
    DLTs = Array{Any}(undef, nSimulations)
    φ_perturb = randn(nSimulations) .* φNoiseFactor
    τ_perturb = randn(nSimulations) .* τNoiseFactor
    σ_perturb = randn(nSimulations) .* σNoiseFactor
    x_c_perturb = randn(nSimulations) .* x_cNoiseFactor
    y_c_perturb = randn(nSimulations) .* y_cNoiseFactor
    z_c_perturb = randn(nSimulations) .* z_cNoiseFactor
    @showprogress 1 "simulating DLT uncertainty from extrinsic params..." for i in 1:nSimulations
        χmodified = ExtrinsicCameraCalibration.ExtrinsicCalibrationParameters(
            χinitial.φ + φ_perturb[i],
            χinitial.τ + τ_perturb[i],
            χinitial.σ + σ_perturb[i],
            χinitial.f,
            χinitial.x_c + x_c_perturb[i],
            χinitial.y_c + y_c_perturb[i],
            χinitial.z_c + z_c_perturb[i])
        #
        DLTs[i] = Calib.CalculateDLTcoefficients(χmodified, ι)
    end
    return DLTs, (; φ_perturb, τ_perturb, σ_perturb, x_c_perturb, y_c_perturb, z_c_perturb)
end
###############################################################################    


###############################################################################
function SimulateAndReturnXYcoordinates(TrueGCPxyz, TrueGCPuv, TrueWaterSurfaceElev, u, v, χinitial, ι;
                                        xyzNoiseFactor = 0.0,
                                        uvNoiseFactor = 0.0,
                                        wseNoiseFactor = 0.0,
                                        AddRadialDistortion = false,
                                        nSimulations = 500,
                                        MaxIterations = 50)
    InfoString = "n Simulations = $nSimulations\nAdding the following types of noise:"
    InfoString = (xyzNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nXYZ $xyzNoiseFactor"])
    InfoString = (uvNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nUV $uvNoiseFactor"])
    InfoString = (wseNoiseFactor == 0.0) ? InfoString : join([InfoString, "\nWSE $wseNoiseFactor"])
    InfoString = join([InfoString, "\nMaxIterations $MaxIterations"])
    @info(InfoString)
    # 
    # Calculate the calibration results without adding any noise:
    χcalculated = Calib.HollandEtAlExtrinsicCalibration(TrueGCPuv, TrueGCPxyz, χinitial, ι,
                                                            MaxIterations=MaxIterations)
    # Calculate the DLTs without adding any noise:
    # DLTcalculated = Calib.CalculateDLTcoefficients(χcalculated, ι)
    DLTcalculated = Calib.CalculateDLTcoefficients(χcalculated, ι)
    # Calculate the x,y coordinates without adding any noise:
    Calculated_x, Calculated_y = Rect.FindXYFromUVZ(u,v, TrueWaterSurfaceElev, DLTcalculated)
    # Generate a set of DLTs that by running the calibraiotn with noise added to
    # the GCP uv and/or xyz coordinates:
    DLTsimulated, GCP_perturb_uv_xyz = SimulateDLTFromGCPsWithAddedNoise(TrueGCPxyz, TrueGCPuv, u, v, χinitial, ι;
                                                                         xyzNoiseFactor = xyzNoiseFactor,
                                                                         uvNoiseFactor = uvNoiseFactor,
                                                                         AddRadialDistortion = AddRadialDistortion,
                                                                         nSimulations = nSimulations,
                                                                         MaxIterations = MaxIterations)
    # Calculate the x,y coordinates of pixels that would be found using the
    # simulated DLTs (i.e., with noise added to the GCP coordinates). Optionally
    # add more noise related to uncertainty in determining the water surface
    # elevation at this step:
    Modified_x, Modified_y, Perturb_wse = CalculateModifiedxyCoordinates(u,v,
                                                                         TrueWaterSurfaceElev, DLTsimulated;
                                                                         wseNoiseFactor=wseNoiseFactor)
    return Calculated_x, Calculated_y, Modified_x, Modified_y# , ρ, θ
end
###############################################################################


###############################################################################
randinrange(minval, maxval) = (minval + rand() * (maxval-minval))
randinrange(minval, maxval, n) = [ randinrange(minval, maxval) for _ in 1:n ]
blurintegers(x::Array) = x .+ rand(size(x)) .- 0.5 # adds a random value b/w (-0.5,0.5), to avoind integer locking
blurintegers(x::Number) = x .+ rand() .- 0.5 # adds a random value b/w (-0.5,0.5), to avoind integer locking
function GenerateGCPuvCoordinatesForSimulation(nGCPs,
                                               u,v,
                                               GCPzrange;
                                               EdgeWidth = 50, # how far in from the edges can GCPs be?
                                               IdealEdge = 10, # distance from edge for fixed, "ideal" GCPs
                                               GCPdistribution = "perimeter") # how to distribute the GCPs throughout the frame
    # 
    ### Define GCPs:
    # 
    # We'll randomly select u,v,z coordinates for the GCPs, possibly with some
    # constraints, and then calculate the true xyz coordinates from these
    # values. This will be the basis for the uncertainty calculations.
    #
    # GCPdistribution = "perimeter" # "left_and_right" # "corners", "perimeter", "anywhere", "left_and_right", "top_and_bottom"
    if GCPdistribution == "anywhere"
        # # ---------- Randomly distributed across image: ----------
        GCPu = randinrange(extrema(u)..., nGCPs)
        GCPv = randinrange(extrema(v)..., nGCPs)
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "corners" 
        # # ---------- Randomly in corners of image: ----------
        # First four GCPs, one in each corner:
        GCPu = [rand(min(u...):EdgeWidth),
                rand((max(u...)-EdgeWidth):max(u...)),
                rand(min(u...):EdgeWidth),
                rand((max(u...)-EdgeWidth):max(u...))]
        GCPv = [rand(min(v...):EdgeWidth),
                rand(min(v...):EdgeWidth),
                rand((max(v...)-EdgeWidth):max(v...)),
                rand((max(v...)-EdgeWidth):max(v...))]
        # Any additional GCPs get distributed among the corners randomly:
        if nGCPs > 4
            GCPu = vcat(GCPu, rand(union(min(u...):EdgeWidth, (max(u...)-EdgeWidth):max(u...)), nGCPs-4))
            GCPv = vcat(GCPv, rand(union(min(v...):EdgeWidth, (max(v...)-EdgeWidth):max(v...)), nGCPs-4))
        end
        GCPu =  blurintegers(GCPu)
        GCPv =  blurintegers(GCPv)
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        # # 
        # trueGCPuvz = hcat(Random.shuffle!(GCPu), Random.shuffle!(GCPv), Random.shuffle!(GCPz))
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "perimeter" 
        # ---------- Randomly along edges of image: ----------
        # set up half of the GCPs along the top/bottom edges, and half along the
        # left/right edges. If an odd number of GCPs is requested, generate
        # ceil(nGCPs/2) of each, leadaing to a total of nGCPs+1, then randomly
        # select nGCPs out of the set.
        # 
        # left/right
        GCPu_lr = rand(union(min(u...):EdgeWidth, (max(u...)-EdgeWidth):max(u...)), ceil(Int, nGCPs/2)) |> blurintegers
        GCPv_lr = randinrange(extrema(v)..., ceil(Int, nGCPs/2))
        # top/bottom:
        GCPu_tb = randinrange(extrema(u)..., ceil(Int, nGCPs/2))
        GCPv_tb = rand(union(min(v...):EdgeWidth, (max(v...)-EdgeWidth):max(v...)), ceil(Int, nGCPs/2)) |> blurintegers
        #
        GCPu = vcat(GCPu_lr, GCPu_tb) 
        GCPv = vcat(GCPv_lr, GCPv_tb)
        GCPz = randinrange(extrema(GCPzrange)..., length(GCPu))
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)[StatsBase.sample(1:size(GCPu,1), nGCPs; replace=false),:]
    elseif GCPdistribution == "perimeter_ideal"
        # Place one GCP in the center of each side of the image (IdealEdge pixels
        # from edge), then add additional GCPs up to nGCPs in the corners:
        GCPu = [round(Int,mean(u)), round(Int,mean(u)),   IdealEdge,          maximum(u)-IdealEdge]
        GCPv = [IdealEdge,          maximum(v)-IdealEdge, mean(v),            mean(v)             ]
        if nGCPs > 4
            GCPu = [GCPu..., IdealEdge]
            GCPv = [GCPv..., IdealEdge]
        end
        if nGCPs > 5
            GCPu = [GCPu..., maximum(u)-IdealEdge]
            GCPv = [GCPv..., maximum(v)-IdealEdge]
        end
        if nGCPs > 6
            GCPu = [GCPu..., IdealEdge]
            GCPv = [GCPv..., maximum(v)-IdealEdge]
        end
        if nGCPs > 7
            GCPu = [GCPu..., maximum(u)-IdealEdge]
            GCPv = [GCPv..., IdealEdge]
        end
        GCPu = round.(Int, GCPu)
        GCPv = round.(Int, GCPv)
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "left_and_right" 
        # ---------- Randomly along left and right edges of image: ----------
        # left/right
        # GCPu = rand(union(min(u...):EdgeWidth, (max(u...)-EdgeWidth):max(u...)), nGCPs) |> blurintegers
        GCPu = union( rand(min(u...):EdgeWidth, floor(Int,nGCPs/2)),
                      rand((max(u...)-EdgeWidth):max(u...), ceil(Int,nGCPs/2))) |> blurintegers
        GCPv = randinrange(extrema(v)..., nGCPs)
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "right" 
        # ---------- Randomly along right edge of image: ----------
        # left/right
        GCPu = rand( (max(u...)-EdgeWidth):max(u...), nGCPs) |> blurintegers
        GCPv = randinrange(extrema(v)..., nGCPs)
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "left" 
        # ---------- Randomly along left edge of image: ----------
        # left/right
        GCPu = rand(min(u...):EdgeWidth, nGCPs) |> blurintegers
        GCPv = randinrange(extrema(v)..., nGCPs)
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "top_and_bottom" 
        # ---------- Randomly along top and bottom edges of image: ----------
        GCPu = randinrange(extrema(u)..., nGCPs)
        GCPv = union(rand(min(v...):EdgeWidth, floor(Int,nGCPs/2)),
                     rand( (max(v...)-EdgeWidth):max(v...), ceil(Int, nGCPs/2))) |> blurintegers
        # GCPv = rand(union(min(v...):EdgeWidth, (max(v...)-EdgeWidth):max(v...)), nGCPs) |> blurintegers
        #
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "top" 
        # ---------- Randomly along top edge of image: ----------
        GCPu = randinrange(extrema(u)..., nGCPs)
        GCPv = rand((max(v...)-EdgeWidth):max(v...), nGCPs) |> blurintegers
        #
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif GCPdistribution == "bottom" 
        # ---------- Randomly along bottom edge of image: ----------
        GCPu = randinrange(extrema(u)..., nGCPs)
        GCPv = rand(min(v...):EdgeWidth, nGCPs) |> blurintegers
        #
        GCPz = randinrange(extrema(GCPzrange)..., nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    elseif lowercase(GCPdistribution) == "none"
        GCPu = fill(NaN, nGCPs)
        GCPv = fill(NaN, nGCPs)
        GCPz = fill(NaN, nGCPs)
        trueGCPuvz = hcat(GCPu, GCPv, GCPz)
    end
    return (trueGCPuv=trueGCPuvz[:,1:2], trueGCPz=trueGCPuvz[:,3])
end
###############################################################################    


###############################################################################    
end                             # end of module
###############################################################################    
