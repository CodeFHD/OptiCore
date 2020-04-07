# OptiCore
Optical Elements Addon for Blender
Current version 2.0 for Blender 2.81

This addon was motivated by the need to create accurate optical elements. Modelling them in Blender via a boolean operation on spheres can lead to artifacts for low f-number surfaces unless very high (u,v)-numbers are used, making the computation very slow.
OptiCore creates surfaces with a fixed number of surface elements for all surface radii, creating accurate surfaces

Currently, it provides two functions, "Add Lens" and "Add Mirror", which can be found in Blender under "Add > Mesh"

For correct optical behaviour, all elements created with this addon should have smooth shading applied, together with edge-split to avoid artifacts at the lens edge. Both operators will be applied by default.

## Add Lens

Creates an optical lens mesh with two flat, spherical or apsheric surfaces.

The component origin is currently fixed at the on-axis intersection with surface 1.

#### Settings

- Surface 1 Type: Shape of first optical surface (spherical or aspheric)

- Surface 1 Radius: Radius of first optical surface, can be negative

- Surface 2 Type: Shape of second optical surface (spherical or aspheric)

- Surface 2 Radius: Radius of second optical surface, can be negative

- Lens Radius: Outer radius of the lens

- Center Thickness: Distance between the on-axis intersection points between the two surfaces. Note that the thickness increases if it becomes incompatible with the surface radii and lens radius

- N1: Number of vertices in radial direction

- N2: Number of vertices in lateral direction

- k: conical constant for surface 1 (only if aspheric is selected)

- A: Aspheric correction coefficients of orders [4,6,8] for surface 1 (only if aspheric is selected)

- k2: conical constant for surface 2 (only if aspheric is selected)

- A2: Aspheric correction coefficients of orders [4,6,8] for surface 2 (only if aspheric is selected)

- Material: Assign pre-defined Blender Material to allow live-view rendering while changing lens parameters

- Smooth Shading: Activate smooth shading

- Use Autosmooth (LuxCore v2.3): If Selected, uses autosmooth on regular mesh. If not selected, uses curstom split normals (only supported in LuxCore v2.4)

- D-shaped lens: Create a cross-section model of the lens

## Add Square Lens

Creates a optical lens mesh with a quadratic shape
Only flat and spherical sufaces are supported at the moment.

The component origin is currently fixed at the on-axis intersection with surface 1.

#### Settings

- Surface 1 Radius: Radius of first optical surface, can be negative

- Surface 2 Radius: Radius of second optical surface, can be negative

- Lens Width: Width of the lens

- Center Thickness: Distance between the on-axis intersection points between the two surfaces. Note that the thickness increases if it becomes incompatible with the surface radii and lens radius

- N: Number of vertices per axis

- Material: Assign pre-defined Blender Material to allow live-view rendering while changing lens parameters

- Smooth Shading: Activate smooth shading

- Use Autosmooth (LuxCore v2.3): If Selected, uses autosmooth on regular mesh. If not selected, uses curstom split normals (only supported in LuxCore v2.4)

## Add Mirror

Creates a circular optical mirror mesh with a spherical or parabolic surface and a flat back.
Parabolic Mirrors can be constructed with component origin at the focal point or mirror centre. For spherical mirrors, only mirror centre is available.
Parabolic mirrors can be 
A hole along the collimated beam axis can be included. 

#### Settings

- Surface Shape: Shape of the Mirror surface. Current Choice is between spherical and parabolic.

- Origin Position: Location of the component origin (only if parabolic is selected). Default is the (on-axis) focal point. Alternatively the center of the  object and on the optical surface.

- Surface Radius: Radius of the optical surface. Can be negative to create a convex mirror.

- Mirror Radius: Outer radius of the mirror

- Back Thickness: Thickness at the thinnest point. Takes into account negative surface radius and offset angle

- Offset Angle: Angle for an off-axis mirror (only if parabolic is selected).

- N1: Number of vertices in radial direction

- N2: Number of vertices in lateral direction

- Material: Assign pre-defined Blender Material to allow live-view rendering while changing mirror parameters. Note: This material must already exits on some other object in Blender, otherwise it despawns when changing any parameter.

- Smooth Shading: Activate smooth shading

- Use Autosmooth (LuxCore v2.3): If Selected, uses autosmooth on regular mesh. If not selected, uses curstom split normals (only supported in LuxCore v2.4)

- Central Hole: Include a hole along the collimated beam direction.

- Hole Radius: Radius of the central hole (only if central hole is selected).

## ToDo-List

### general programming
- Update Readme for Square Lens
- Square lens needs triangles as surface elements (Done). Parabolic mirror might also for off-axis, check!
- check for nonsense-geometries
- If faces touch, don't generate outer rim to avoid duplicate vertices and split-normals
- test split-normals vs. edge-split artefacts
- Cross-section models (D-shaped - done for regular lens), perhaps with variable cut width
- lens file import (seq, zmx, spd; open standards, patents, copyright?)
- lens system designer

### Optical Components to add
- flat annulus
- along-focus holes for OAPs
- elliptics (off-axis flat)
- aspheric surfaces for lenses and mirrors
- Prisms (right-angle, triangle, dove, roof, CC-RR, cats-eye-RR, penta, wedge[?])
- cylindrical lenses/mirrors
- Fresnel lens
- wedged window (application without interference?)
- Axicons/laser line generators
- Further types of aspheric surfaces
-- Schmidt corrector plate
-- free-form with input formula if possible
- Microlens array (test first if usable)
- doublets/triplets with well-defined interior surface