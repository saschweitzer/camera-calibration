###############################################################################
# 
# This script contains Julia code for performing uncertainty analysis of camera
# calbiration and georeferencing, as described in Schweitzer & Cowen, 2021; and
# Schweitzer & Cowen, 2022
# 
# # References
# Schweitzer, S. A., & Cowen, E. A. (2021). "Instantaneous river-wide water surface
# velocity field measurements at centimeter scales using infrared quantitative
# image velocimetry". Water Resources Research, 57,
# e2020WR029279. https://doi.org/10.1029/2020WR029279
# 
# Schweitzer, S. A., & Cowen, E. A. (2022, exptected) "A Method for Analysis of
# Spatial Uncertainty in Image Based Surface Velocimetry". Frontiers in Water,
# Water and Hydrocomplexity.
# 
# # License
# This software is released under the MIT license, however, the authors request to
# be credited if you use the software, preferably by citing the above-mentioned
# Schweitzer & Cowen papers.
# 
# MIT License
# 
# Copyright (c) 2022 Seth A. Schweitzer
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# 
###############################################################################

###############################################################################
using Statistics
#
import ExtrinsicCameraCalibration as Calib
import ImageRectification as Rect
import GeoreferencingUncertaintySimulations as GeorefSim
#
dropdim3(x) = dropdims(x;dims=3)
dim3quantiles(data; q=0.95) = mapslices(d->quantile(d, q), data; dims=3) |> d -> dropdims(d; dims=(findall(size(d).==1)[]))
###############################################################################

###############################################################################
# 
### Define camera position and orientation :
# 
## Camera optical center position:
x_c = 0.0
y_c = 0.0
z_c = 25.0
## Camera orientation:
# Generic setup is 45deg tilt for oblique view, 0.0deg tilt for UAV.
φ = deg2rad(0.0) # φ azimuth. clockwise from the north (y-axis). This is meaningless in a camera-centric coordinate system
τ = deg2rad(45.0) # τ tilt. angle away from the vertical of the camera centeral ray (0 is nadir view).
σ = deg2rad(0.0) # σ roll. rotation of the camera around the centeral ray. (0.0 is no roll, rotation is from the positive x axis)
###############################################################################
# Define actual water surface elevation:
trueWSE = 0.0
# Number of GCPs to use (may later add restrictions on their distribution):
nGCPs = 4 # 6 # 8
GCPzrange = (trueWSE, trueWSE+3.0) # we'll randomly generate a value in this range for each GCP
#
GCPdistribution =  "corners" # "left_and_right" # "corners", "perimeter", "perimeter_ideal", "anywhere", "left_and_right", "top_and_bottom"
###############################################################################
SimulationParametersInfo = join(["φ=$(round(rad2deg(φ),digits=1))°",
                                 "τ=$(round(rad2deg(τ),digits=1))°", 
                                 "σ=$(round(rad2deg(σ),digits=1))°",
                                 "z_c=$(z_c)",
                                 "$(nGCPs) GCPs",
                                 GCPdistribution], ", ")
###############################################################################

###############################################################################
# 
### GCP perturbation parameters:
#
# "NoiseFactors" represent the standard deviation of the distribution from which
# perturbations will be randomly drawn:
xyzNoiseFactor = [0.03,0.03,0.03]
uvNoiseFactor  = [0.5,0.5]
wseNoiseFactor =  0.03
#
φNoiseFactor = deg2rad(0.8)
τNoiseFactor = deg2rad(0.1)
σNoiseFactor = deg2rad(0.1)
x_cNoiseFactor =  0.02           # horizontal position (RTK)
y_cNoiseFactor =  0.02           # horizontal position (RTK)
z_cNoiseFactor = 0.05           # meters
# 
## Do we want to simulate radial distortion?
AddRadialDistortion = true     # false
# 
nSimulations = 700
###############################################################################

###############################################################################
# 
## Define camera intrinsic and extrinsic parameters, and DLT coefficients:
# Here we use the parameters of the camera used in Schweitzer & Cowen, 2021 WRR.
#
CameraPixelResolution = [1344 784]
LensFocalLength_m = 17e-3
# 
ιtrue = Calib.IntrinsicCalibrationParameters(mean([1 CameraPixelResolution[1]]), # u₀
                                             mean([1 CameraPixelResolution[2]]), # v₀
                                             1.0,   # λᵤ
                                             1.0,   # λᵥ
                                             14e-6, # nominal pixel pitch
                                             [1, 3, 5], # radial_distortion_poly_degrees
                                             [0.020522291524825344, -3.156396095419575e-8, -8.859823756718999e-14]) # radial_distortion_poly_coeffs
# 
χtrue = Calib.ExtrinsicCalibrationParameters(
    φ, τ, σ, 
    -LensFocalLength_m/ιtrue.PixelPitch,    # f; distance b/w FPA and optical center of camera. Focal length/pixel pitch. Note that this may be a negative number (since the FPA is behind the focal point)
    x_c, y_c, z_c,
)
#
trueDLT = Calib.CalculateDLTcoefficients(χtrue, ιtrue)
###############################################################################

###############################################################################
# Work on only a subset of pixels in the image, to improve performance:
row_skip_px = 6                 # 2, 4, 8, etc.
col_skip_px = 6                 # 2, 4, 8, etc.
###############################################################################

###############################################################################
# Set up the pixel grid over which we will run the calibration sensitivity analysis:
u = [u for v in 1:row_skip_px:CameraPixelResolution[2], u in 1:col_skip_px:CameraPixelResolution[1]];
v = [v for v in 1:row_skip_px:CameraPixelResolution[2], u in 1:col_skip_px:CameraPixelResolution[1]];
###############################################################################

###############################################################################
# 
### Generate the simulated GCP:
# 
# random uvz coordinates:
trueGCPuv, trueGCPz = GeorefSim.GenerateGCPuvCoordinatesForSimulation(
    nGCPs,
    u,v,
    GCPzrange;
    EdgeWidth = 50, # how far in from the edges can GCPs be?
    IdealEdge = 10, # distance from edge for fixed, "ideal" GCPs
    GCPdistribution = GCPdistribution) # how to distribute the GCPs throughout the frame. defined above
# -----------------------------------------------------------------------------
# Calculate GCP x,y coordinates from the uvz coordinates:
trueGCPx, trueGCPy = Rect.FindXYFromUVZ(trueGCPuv[:,1], trueGCPuv[:,2], trueGCPz, trueDLT)
trueGCPxyz = hcat(trueGCPx, trueGCPy, trueGCPz)
###############################################################################

###############################################################################
# 
# Calculate the "reference grid" of true pixel footprints:
# 
χcalculated = Calib.HollandEtAlExtrinsicCalibration(trueGCPuv, trueGCPxyz, χtrue, ιtrue; MaxIterations = 50)
DLTcalculated = Calib.CalculateDLTcoefficients(χcalculated, ιtrue)
XCalculated, YCalculated = Rect.FindXYFromUVZ(u,v, trueWSE, DLTcalculated)
###############################################################################

###############################################################################
# 
# We can now define the complete set of simulation input parameters:
# 
GCPSimulationInputs = (; trueGCPxyz, trueGCPuv, trueWSE, u, v, χcalculated, ιtrue)
ExtrinsicParamsSimulationInputs = (; χtrue, ιtrue, trueWSE, u, v)
# 
# And the set of uncertainty parameters to use for the simulations:
GCPSimulationUncertaintyParams = (; xyzNoiseFactor,
                                  uvNoiseFactor,
                                  wseNoiseFactor,
                                  AddRadialDistortion,
                                  nSimulations, # how many simulation scenarios to run
                                  MaxIterations = 50, # max iterations of calibration to run before giving up
                                  )
#
ExtrinsicParamsSimulationUncertaintyParams = (; φNoiseFactor, τNoiseFactor, σNoiseFactor,
                                              x_cNoiseFactor,
                                              y_cNoiseFactor,
                                              z_cNoiseFactor,
                                              wseNoiseFactor, #
                                              AddRadialDistortion,
                                              nSimulations, # how many simulation scenarios to run
                                              MaxIterations=50) # max iterations of calibration to run before giving up
###############################################################################

###############################################################################
# 
# Apply perturbations, and calculate results!
# 
###############################################################################

###############################################################################
# Perturb xyz and uv coordinates of GCPs, and WSE:
ObliqueView_Outcomes_xyz_uv_wse = GeorefSim.GeorefGCPSimulationWrapper(
    GCPSimulationInputs...;
    GCPSimulationUncertaintyParams...,
    nSimulations = nSimulations*2) # 
###############################################################################

###############################################################################
# Perturb only xyz coordinates of GCPs:
ObliqueView_Outcomes_xyz = GeorefSim.GeorefGCPSimulationWrapper(
    GCPSimulationInputs...;
    GCPSimulationUncertaintyParams...,
    uvNoiseFactor = 0.0,
    wseNoiseFactor = 0.0)
###############################################################################

###############################################################################
# Perturb only uv coordinates of GCPs:
ObliqueView_Outcomes_uv = GeorefSim.GeorefGCPSimulationWrapper(
    GCPSimulationInputs...;
    GCPSimulationUncertaintyParams...,
    xyzNoiseFactor = 0.0,
    wseNoiseFactor = 0.0)
###############################################################################

###############################################################################
# Perturb only WSE:
ObliqueView_Outcomes_wse = GeorefSim.GeorefGCPSimulationWrapper(
    GCPSimulationInputs...;
    GCPSimulationUncertaintyParams...,
    xyzNoiseFactor = 0.0,
    uvNoiseFactor = 0.0)
###############################################################################

###############################################################################
# Perturb Direct Georeferencing parameters:
ObliqueView_Outcomes_Extrinsic = GeorefSim.GeorefExtrinsicParamsSimulationWrapper(
    ExtrinsicParamsSimulationInputs...;
    ExtrinsicParamsSimulationUncertaintyParams...,
    nSimulations = nSimulations*2,
)
###############################################################################

###############################################################################
###############################################################################
#
# Gradients of uncertainty from here!!
#
###############################################################################
###############################################################################

###############################################################################
""" Calculate gradients in 1 or 2 dimensions, using a central differencing scheme. """
function central_diff(data_array::AbstractMatrix;
                      dims=nothing,
                      grid_spacing_1=1.0,
                      grid_spacing_2=1.0)
    if any(dims .== 1)
        d1    = diff(data_array, dims=1)/2
        _g1   = cat(d1[[1],:], d1; dims=1) # put the first row back in to recover the original shape
        _g1 .+= cat(d1, d1[[end],:]; dims=1)
        _g1 ./= grid_spacing_1
    end
    if any(dims .== 2)
        d2    = diff(permutedims(data_array), dims=1)/2
        _g2   = cat(d2[[1],:], d2; dims=1) # put the first row back in to recover the original shape
        _g2 .+= cat(d2, d2[[end],:]; dims=1)
        _g2 = permutedims(_g2)
        _g2 ./= grid_spacing_2
    end
    if dims == 1
        return _g1
    elseif dims == 2
        return _g2
    else
        return _g1 .+ _g2 .* im # return complex number, where _g1 is real part, _g2 is imaginary
    end
end
###############################################################################




###############################################################################
"""
Calculate the 2D gradient of difference between reference and simulated
coordinates (i.e., the georeferencing error), at each pixel footprint
location. Results are complex numbers at each node, in which the real part
represents dim1 (y direction) and the imaginary part represents dim2 (x
direction).
"""
function Calculate2DCoordinateGradient(SimulationOutcomes)
    # Calculate the (irregular) "grid" spacing along each coordinate:
    grid_spac1 = central_diff(SimulationOutcomes.Calculated_y[:,:]; dims=1) # distance between y coordinate of grid-adjacent points (y is dim 1)
    grid_spac2 = central_diff(SimulationOutcomes.Calculated_x[:,:]; dims=2) # distance between x coordinate of grid-adjacent points (x is dim 2)
    #
    # Calculate the gradients of uncertainty along each dimension, for each simulation:
    # y coordinate (varies along dim 1!)
    _y_coordinate_gradient = mapslices(y->central_diff(y; dims=(1,2),
                                                          grid_spacing_1=grid_spac1,
                                                          grid_spacing_2=grid_spac2),
                                       SimulationOutcomes.Modified_y .- SimulationOutcomes.Calculated_y;
                                       dims=(1,2))
    # Note: this is 2D, despite only considering the y-direction, because the
    # reference points are not laid out on a regular grid, and hence "grid-adjacent"
    # nodes have different coordinates in two dimensions. We calculate the slope
    # with which the difference between (1D) simulated and reference "y" coordinate
    # changes along each of the two directions (x,y), at each location on the
    # (reference) grid.
    # 
    # Similar to above, but for the "x" coordinate (varies along dim 2!):
    _x_coordinate_gradient = mapslices(x->central_diff(x; dims=(1,2),
                                                          grid_spacing_1=grid_spac1, 
                                                          grid_spacing_2=grid_spac2),
                                       SimulationOutcomes.Modified_x .- SimulationOutcomes.Calculated_x;
                                       dims=(1,2))
    # 
    # return an array of complex numbers, in which the real component is the
    # total gradient in the y direction, and the the imaginary component is the
    # gradient in the x direction:
    return (; _x_coordinate_gradient, _y_coordinate_gradient)
end
###############################################################################

###############################################################################
# 
# Calculate the gradients of uncertainty resulting from spatial uncertainty:
# 
ObliqueView_UncertGradients_xyz_uv_wse = Calculate2DCoordinateGradient(ObliqueView_Outcomes_xyz_uv_wse)
ObliqueView_UncertGradients_xyz = Calculate2DCoordinateGradient(ObliqueView_Outcomes_xyz)
ObliqueView_UncertGradients_uv = Calculate2DCoordinateGradient(ObliqueView_Outcomes_uv)
ObliqueView_UncertGradients_wse = Calculate2DCoordinateGradient(ObliqueView_Outcomes_wse)
ObliqueView_UncertGradients_Extrinsic = Calculate2DCoordinateGradient(ObliqueView_Outcomes_Extrinsic)
###############################################################################
