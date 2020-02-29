# OptiCore
Optical Elements Addon for Blender
Current version 2.0 for Blender 2.81

This addon was motivated by the need to create accurate optical elements. Modelling them in Blender via a boolean operation on spheres can lead to artifacts for low f-number surfaces unless very high (u,v)-numbers are used, making the computation very slow.
OptiCore creates surfaces with a fixed number of surface elements for all surface radii, creating accurate surfaces

Currently, it provides two functions, "Add Lens" and "Add Mirror", which can be found in Blender under "Add > Mesh"

For correct optical behaviour, all elements created with this addon should have smooth shading applied, together with edge-split to avoid artifacts at the lens edge. Both operators will be applied by default.

## Add Lens

Creates an optical lens mesh with two spherical surfaces.
Planned to be updated in future with additional features such as aspheric surfaces.

The component origin is currently fixed at the on-axis intersection with surface 1.

#### Properties

Surface 1 Radius: Radius of first optical surface, can be negative

Surface 2 Radius: Radius of second optical surface, can be negative

N1: Number of vertices in radial direction

N2: Number of vertices in lateral direction

Lens Radius: Outer radius of the lens

Center Thickness: Distance between the on-axis intersection points between the two surfaces. Note that the thickness increases if it becomes incompatible with the surface radii and lens radius

Material: Assign pre-defined Blender Material to allow live-view rendering while changing lens parameters. Note: This material must already exits on some other object in Blender, otherwise it despawns when changing any parameter

Smooth Shading: Activate smooth shading

Edge Split: Activate edge split

## Add Mirror

Creates a circular optical mirror mesh with a spherical surface and a flat back.

Surface Shape: Shape of the Mirror surface. Current Choice is between spherical and parabolic.

Origin Position: Location of the component origin. Default is the (on-axis) focal point. Alternatively the center of the  object and on the optical surface.

Surface Radius: Radius of the optical surface. Can be negative to create a convex mirror.

N1: Number of vertices in radial direction

N2: Number of vertices in lateral direction

Mirror Radius: Outer radius of the mirror

Back Thickness: Thickness at the thinnest point. Takes into account negative surface radius and offset angle

Offset Angle: Angle for an off-axis mirror.

Material: Assign pre-defined Blender Material to allow live-view rendering while changing mirror parameters. Note: This material must already exits on some other object in Blender, otherwise it despawns when changing any parameter.

Smooth Shading: Activate smooth shading

Edge Split: Activate edge split