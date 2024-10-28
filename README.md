


<p align="center">
<h1 align="center">OptiCore</h1>
</p>
<p align="center">
Optical Elements Addon for Blender
<br />

OptiCore is a Blender-addon to provide accurate and well-defined models of optical elements (lenses, mirrors, etc.) as well as a few selected optimechanical components (optical bench, posts). These models can not only help you to easily create renderings of optical laboratory experiments, but also stunning images involving caustics.

This addon will primarily be maintained for the latest Blender release. The last commit was tested on Blender version 4.2.1 LTS.

</p>

<p align="center">
<a href="http://www.youtube.com/watch?feature=player_embedded&v=D8rQBVI4lIg
" target="_blank"><img src="http://img.youtube.com/vi/D8rQBVI4lIg/0.jpg" 
alt="Introduction Video on YouTube" width="480" height="270" border="10" /></a>

</p>


The development of this addon was originally motivated by the need to create accurate visualisations of optical lenses. Modelling them in Blender via spheres and boolean operations can lead to artifacts for surfaces where only a small segemnt of the sphere is needed, unless very high (u,v)-number spheress are used. This, in turn, makes it very slow and RAM-intensive.
With this addon, lenses are created using parametrisations that are used in lesn design. All computations are internally handled with 64-bit precision, including surface normals with edge-splits for the best possible performance.
The name Opti"Core" is derived from the LuxCoreRender engine, which is a great rendering engine and inspired me to continue using Blender for lens system visualisations.
It should be noted however, that OptiCore is completely independent from any rendering engine (including LuxCore), as its main purpose is to generate the geometric objects. An internal sequential ray-tracing feature is included (see below for further explanation), which is however not connected to Blenders rendering engines. Handling of materials of a third-party rendering engine, or the rendering engines included in Blender by default, are out of the scope of this addon at the present time.

Below, you will find a brief description of the available optical elements and features. You can also find tutorial videos on my YouTube-channel [HowToPrint](https://www.youtube.com/@howtoprint6002). Please don't hesitate to open a GitHub-Issue or -Discussion if you have any questions or feature requests.

## Table of Contents
* [Optical Elements](#optical-elements)
  * [Lens](#lens)
  * [Square Lens](#square-lens)
  * [Mirror](#mirror)
  * [Retroreflector](#retroreflector)
  * [Retroreflector](#siemens-star)
  * [Retroreflector](#zemax-file-import)
* [Optomechanics](#optomechanics)
  * [Breadboard](#breadboard)
  * [Post](#post)

## Optical Elements

### Lens

Creates an optical lens with flat, spherical or apsheric surfaces. Set surface radius = 0 for a flat surface.

Singlet, Doublet and Triplet lenses are supported. (More than three elements cemented together are uncommon and are not covered at the moment due to the implications on code maintenance.)

The component origin is placed at the on-axis intersection with the first surface.

Aspheric surfaces are defined with a conical constant and three coefficients for polynominal terms. See the respective Wikipedia-page for an explanation: <https://en.wikipedia.org/wiki/Aspheric_lens>

The lens can also be created as a cross-section model.

A ray fan can be added to visualize the optical path of rays through the lens. This uses an internal, sequential ray tracing algorithm. Presently, flat, spherical and conic surfaces are supported; surfaces using the polynominal aspheric coefficients are not yet supported.

### Square Lens

Creates a optical lens mesh with a quadratic outline.

The feature set is currently reduces compared to regular lenses, i.e. spherical surfaces only and no option for a cross-section model.

### Mirror

Creates an optical mirror mesh with a circular outline, spherical, parabolic or aspheric surface and a flat back.

Parabolic Mirrors can be constructed with component origin at the focal point or mirror centre. For spherical mirrors, only the mirror centre is available at the moment.

Parabolic mirrors can be constructed as off-axis parabolic mirrors.

A hole along the collimated beam axis can be included.

### Retroreflector

Creates an array of cubecorner retroreflectors.

Features two "Retroreflector types":
- "Cubecorner", which is how reflectors are typically manufactured, and
- "Trirectangular tetrahedron", which is basically only the lower half. (Uncommon, but why not have it, looks cool...)

Retroreflectors can be used in two ways:
- with a mirror material, in which case you look directly onto the cubes.
- with a glass/plastic material, in which case you look from the flat side. Retroreflection then occurs due to total internal reflection.

A "Tip offset" can be specified to apply an offset to the base bottom vertices of the cube-corners. The range of [-1,1] corresponds to 5% of the "CubeCorner Spacing". This allows to mimic by-design imperfections of retroreflectors for traffic use - if they were perfect, light would be reflected from headlights back into headlight and not reach the human eye.

Hint: slightly bevel all edges to give a realistic segemented look due to real imperfections of sharp corners.

### Siemens Star

Creates a Siemens Star, a typical MTF or contrast test object.

The structure can be created open ended, or surrounded by a solid cylinder wall and the Siemens Star cut out as a negative Structure.

### Zemax file import

Zemax is a major professional optical design software. Its .zmx file format is common for exchanging optical designs, and many examples, e.g. from patent literature, are available in this format.
The description of the file format is not published. It is in ASCII-format, so human readable, but based heavily on abbreviations and therefore not easy to interpret fully. OptiCore therefore features an import routine that implements the basic features, such as surface curvatures, standard aspheric surface, or mirrors. However, there may be many aspects that lead to an incorrect import into Blender.

Besides the correct interpretation of the .zmx file, the ray tracing may fail because it contains obsolete optical glasses that are not in the included database. OptiCore currently includes a library that is generated from available, current data from glass manufacturers.

In combination of these aspects, it is suggested to first check the internal, sequential ray-tracing to determine if a good focus spot is rendered and the system looks plausible. If not, checking the Blender system console can tell you if there are glass library errors. In this case, transfer to Blender materials may still be successful if you can identify the mateiral properties yourself.

## Optomechanics

### Breadboard

Creates a rectangular plate with beveled holes and rounded corners, i.e. the top-plate of an optical table.

Upon creation, the faces inside the holes are pre-selected. This allows you to easily assign a second material, e.g. to create a screw thread by shader.

### Post

Creates an 0.5-inch optical post (Thorlabs-style), with a 4 mm hole at the top and 6 mm hole at the bottom, as well as a through hole.

Warning: This mesh geometry with a through-hole appears to be very complicated to properly shade smooth. Artifacts may occur with high glossiness. Autosmooth appears to work better than split normals.

Other post-sizes are to follow.