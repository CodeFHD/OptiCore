# OptiCore
Optical Elements Addon for Blender

This addon was motivated by the need to create accurate optical elements. Modelling them in Blender via a boolean operation on spheres can lead to artifacts for low f-number surfaces, if not very high u,v-numbers are used, making the computation very slow.
This addon creates surfaces with a fixed number of surface elements for all surface radii, creating accurate surfaces

Currently, it provides two functions, "Add Lens" and "Add Mirror", which can be found in Blender under "Add > Mesh"

All elements created with this addon should normally have an edge-split modifier and smooth shading applied.

## Add Lens

Creates an optical lens mesh with two spherical surfaces.
Planned to be updated in future with additional features such as aspheric surfaces.

The component origin is set at the on-axis intersection with surface 1.

#### Properties

Surface 1 Radius: Radius of first optical surface, can be 

Surface 2 Radius: Radius of second optical surface

N1: Number of vertices in radial direction

N2: Number of vertices in lateral direction

Lens Radius: Outer radius of the lens

Center Thickness: Distance between the on-axis intersection points between the two surfaces

## Add Mirror

Creates an optical mirror mesh with a spherical surface and a flat back.

(Currently) Equivalent to Add Lens with surface radius 2 = 0 and surface radius 1 multiplied by -1
Planned to be updated in future with additional features such as parabolic surfaces and off-axis surfaces