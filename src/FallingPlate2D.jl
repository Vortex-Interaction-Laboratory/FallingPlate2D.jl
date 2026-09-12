"""
FallingPlate2D — a thin layer over WaterLily.jl for a freely falling 2D plate
section in a quiescent fluid, solved in a frame translating with the plate.

Units: lengths in grid cells, fluid density 1, velocity scale `U = 1` equal to
the gravitational speed `sqrt(2(ρ − 1) f g)`, where `f` is the thickness.

Frames. Lab coordinates `(X, Y)`: `X` horizontal, `Y` vertical upwards, gravity
along `−Y`. Grid coordinates `(x, y) = (Y, −X)` (same handedness): the plate
falls towards `−x`, so the far-field flow in the plate frame is along `+x`,
enters at `x = 0` and leaves through the convective exit at `x = xmax`; the
lateral direction `y` is periodic. The body angle `θ` is anticlockwise about
`z`; the chord line makes the angle `φ = θ − π/2` with the horizontal,
positive anticlockwise, so `φ = 0` is broadside to the fall.
"""
module FallingPlate2D

using WaterLily, StaticArrays, LinearAlgebra, Printf

export Ellipse, RoundedPlate, area, inertia, added_mass, andersen_Istar
export RigidState, PlateDynamics, FallingPlateSim, falling_plate, step!, run!, rigid_velocity_fill!
export Gust, NOGUST, gust_velocity, gust_acceleration, gust_displacement, peak_acceleration
export summarise, to_lab, to_grid, apply_darcy!, prepare_flow!, screening_length

include("bodies.jl")
include("gust.jl")
include("rigid.jl")
include("darcy.jl")
include("simulation.jl")
include("diagnostics.jl")

end
