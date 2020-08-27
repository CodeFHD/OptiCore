


<p align="center">
<h1 align="center">OptiCore</h1>
</p>
<p align="center">
Optical Elements Addon for Blender
<br />

OptiCore is a Blender-addon to provide accurate and well-defined models of optical elements (lenses, mirrors, ...) as well as optimechanical components (optical bench, posts, lens mounts, ...). These models can not only help you to easily create models of optical laboratory experiments, but also offer a functional implementation of optical elements to create stunning images involving caustics.

</p>


The addon was motivated by the need to create accurate optical lenses. Modelling them in Blender via a boolean operation on spheres can lead to artifacts for low f-number surfaces unless very high (u,v)-numbers are used, making the computation very slow.
The name Opti"Core" is derived from the LuxCoreRender engine, which is a great tool for physically accurate rendering and also available as a powerful addon for Blender.

You can find a brief description of the available optical elements below. Please don't hesitate to contact me or open a GitHub-Issue if you have any requests.

For correct optical behaviour, all elements created with this addon should have smooth shading applied, together with the edge-split modifier, auto-smooth mesh or custom split normals. There appear to be differences between these methods, and I still have to figure out which offers the best behaviour. (Input and user feedback welcome!)

## Table of Contents
* [Optical Elements](#optical-elements)
  * [Lens](#lens)
  * [Square Lens](#square-lens)
  * [Mirror](#mirror)
* [Optomechanics](#optomechanics)
  * [Breadboard](#breadboard)
  * [Post](#post)
* [ToDo List](#todo-list)

## Optical Elements

### Lens

Creates an optical lens with two flat, spherical or apsheric surfaces. Set surface radius = 0 for a flat surface.

The component origin is placed at the on-axis intersection with surface 1.

Aspheric surfaces are defined with a conical constant and three coefficients for polynominal terms. See the respective Wikipedia-page for an explanation: <https://en.wikipedia.org/wiki/Aspheric_lens>

An option exists to create a cross-section model (D-shaped lens).

### Square Lens

Creates a optical lens mesh with a quadratic outline

### Mirror

Creates an optical mirror mesh with a circular outline, spherical or parabolic surface and a flat back.

Parabolic Mirrors can be constructed with component origin at the focal point or mirror centre. For spherical mirrors, only the mirror centre is available at the moment.

Parabolic mirrors can be constructed as off-axis parabolic mirrors.

A hole along the collimated beam axis can be included.

###Retroreflector

Creates an array of cubecorner retroreflectors.

Features two "Retroreflector types":
- Cubecorner, which is how reflectors are typically made, and
- "Trirectangular tetrahedron", which is basically only the lower half. (Uncommon, but why not have it, looks cool...)

Retroreflectors can be used in two ways:
- with a mirror material, in which case you look directly onto the cubes.
- with a glass/plastic material, in which case you look from the flat side. Retroreflection then occurs due to total internal reflection.

A "Tip offset" can be specified to apply an offset to the base bottom vertices of the cube-corners. The range of [-1,1] corresponds to 5% of the "CubeCorner Spacing". This allows to mimic by-design imperfections of retroreflectors for traffic use - if they were perfect, light would be reflected from headlights back into headlight and not reach the human eye.

Hint: slightly bevel all edges to give a realistic segemented look due to real imperfections of sharp corners.

## Optomechanics

### Breadboard

Creates a rectangular plate with beveled holes and rounded corners, i.e. the top-plate of an optical table.

Upon creation, the faces inside the holes are pre-selected. This allows you to easily assign a second material, e.g. to create a screw thread by shader.

### Post

Creates an 0.5-inch optical post (Thorlabs-style), with a 4mm hole at the top and 6mm hole at the bottom, as well as through hole.

Warning: Mesh geometry with the through-hole appears to be very complicated to properly shade smooth. Artifacts may occur with high glossiness.

Other post-sizes are to follow.

## ToDo List

### Optical Components To Add
- flat annulus around concave lens surfaces
- elliptics (off-axis flat mirrors)
- Prisms (right-angle, triangle, dove, roof, cats-eye-RR, penta, wedge[?])
- cylindrical lenses/mirrors
- Fresnel lenses
- wedged window
- Axicons/laser line generators
- Further types of aspheric surfaces
-- Schmidt corrector plate
-- free-form with input formula; if feasible
- Microlens array (test if at all usable with 32bit rendering)
- doublets/triplets with well-defined interior surface

### Optical Component Improvements
- All: Consistent default settings of size, i.e. mm vs m scale
- Lens: center thickness checks: small surface radii as well as aspheric
- Lens: If faces touch, don't generate one outer rim to avoid duplicate vertices and split-normals
- Lens: Origin at focus point (at least infinity focus)
- Mirror: along-focus holes for OAPs
- CC-RR: manufacturing-accurate closure at sides
- CC-RR: Include Bevel option
- CC-RR: Tip offset presets

### Other Programming
- Cross-section models (D-shaped - done for regular lens), perhaps with variable cut width
- lens file import (seq, zmx, spd - open standards? patents? copyright?)
- lens system designer GUI

### Performance Assessments

- split-normals vs. edge-split artefacts
- general limitations of 32bit ray tracing