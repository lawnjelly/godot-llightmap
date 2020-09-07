# godot-llightmap
* Lightmap module for Godot Engine 3.2.2 or later
* Version 0.31 (September 7th, 2020)
* Lightmaps created can be used with standard Godot builds and templates, i.e. you only need the module for a preprocess
* (work in progress, there may be bugs, especially in the uv mapping but it is usable)

## Precompiled builds
* Go to https://github.com/lawnjelly/godot-titan
* Click the releases and get the latest version

<img src="/pics/llightmap_01screen.png" alt="screenshot" width="400"/>

## About
While JFons is finishing off the new official core lightmapper for 3.2, I've spent some time playing making this toy lightmapper. It works so well I'm now going to be using it to replace the lightmapping functionality in my other module, lportal.

LLightmap is designed to be easy to use and produce nice soft shadows, with high performance at runtime. The existing core lightmappers produce a different lightmap for each object. This can be slow due to texture swapping. LLightmap instead can uv map your scene with multiple objects, unwrapping so that all your objects share the same texture space. Thus a single lightmap can be used.

https://www.youtube.com/watch?v=m7fJTWoVYj0

https://www.youtube.com/watch?v=pBpF2raGA8A

### Features
* Forward and backward ray tracing options
* Volumetric lights for soft shadows
* Ambient occlusion
* Linear HDR exr intermediate files
* Linear HDR exr final texture or gamma corrected normalized png
* SSE2 (on x86 64 bit) / Multithread ray tracing
* Albedo taken into account for bounces
* Omnis, Spotlights, Directional lights
* Emissive materials
* Transparency
* Light Probes for realtime lighting of dynamic objects

### Still todo
* Sky panorama textures
* Roughness / metal from source textures for PBR reflections
* Option of multiple lightmaps

## Instructions
### Shader
Before we begin, we need to understand that LLightmap simply produces a texture. It doesn't automatically show it on objects. In order to show our lightmaps, we should use a custom shader on our static level meshes.

An example shader is as follows:
```
shader_type spatial;

// we are using a lightmap, we don't need realtime lighting
render_mode unshaded;

// these 2 are optional, and although unused in the shader,
// allow us to set materials to emit light in the lightmapping stage
uniform float emission;
uniform vec4 emission_color : hint_color;

// our input textures, a material texture, and the lightmap
uniform sampler2D texture_albedo : hint_albedo;
uniform sampler2D texture_lightmap : hint_albedo;

void fragment() {
	// lookup the colors at the uv location of our textures
	vec4 albedo_tex = texture(texture_albedo,UV);
	vec4 lightmap_tex = texture(texture_lightmap,UV2);
  
	// the overall albedo (color) will be the material texture TIMES the lightmap
	// (so it can be darkened).
	// you can optionally use a multiplier to allow lightening areas (the 2.0 here)
	ALBEDO = albedo_tex.rgb * lightmap_tex.rgb * 2.0;
}
```
You may end up using variations of this shader in practice, but it will get you started.

When you build your level geometry you should assign this shader, and in the shader parameters, set which texture to use for the mesh, and create a dummy png which will be our final lightmap, and assign this to the lightmap slot.

#### Initial setup
1) Make sure your scene to be lightmapped is under a spatial, name it e.g. 'level' (not the scene root).
2) Mark any static geometry inside the level to be lightmapped with the 'use in baked light' flag in the geometry section of the mesh.
3) Either place the lights as part of the level branch, or create another branch just for static lights.
4) Add the new 'LLightmap' node to the scene tree (preferably not in the level or lights branch).

If you click on the LLightmap inspector, we can now set things up for lightmapping.
1) Assign the 'meshes' to point to your level geometry branch.
2) Assign the 'light' to point to your lights branch.
3) Set the filenames you wish to use for lights intermediate, ambient occlusion intermediate, and the final combined lightmap.

#### Unwrapping
Before we can lightmap anything, we need to make sure the level geometry is uvmapped, usually in the 2nd uv channel. The first uv channel is typically used for the main texture, and we will use the 2nd uv for the lightmap.

You can uv map scenes in a third party modelling program such as blender, but this is inconvenient, and LLightmap can use xatlas to do this for us.

1) All the meshes to be unwrapped should have `Use_In_Baked_Light` set to on in their Geometry section. Meshes that aren't tagged will not be lightmapped (or affect the lightmap of the other geometry).
2) In order to add the new UVs, we will need to modify the level scene, so it makes sense to save it as a new tscn file. You should set the filename as 'uv_filename' in the uv unwrap section of the inspector.
3) Set the bake_mode to 'UVMap'.
4) Once this output filename is set, and the meshes is assigned correctly, and the bake mode is correct, hit the 'Bake Lightmap' button above the main 3d window.
5) This does a number of things by magic. First it merges all the marked geometry into a single mesh, then it unwraps the mesh, then it 'unmerges' the wrapped mesh back to the original objects. This is quite a complex process and can result in added vertices. Finally it saves the new uvmapped scene into the file we specified.
6) As the mesh data has been altered, the original level mesh is deleted, and you should load in its place the UVmapped level that was exported. It is *highly recommended* to restart Godot IDE before loading in the UVmapped level due to referencing bugs in Godot core.

A backup of the original branch is saved to `uvmap_backup.tscn` in your project folder, just in case. You are highly recommended to keep a backup of your original level before uvmapping, for further editing etc.

#### Baking
Once the scene is uvmapped, we can move onto the fun stage, baking some lightmaps.

1) Change the bake_mode to 'lightmap'.
2) Once we are sure we have some lights, the meshes and lights are assigned, and the filenames are chosen, hit the 'Bake Lightmap' button.
3) This may take a while to complete, depending on the settings you use. The defaults should bake quite quickly for a very rough version if you choose 'low' or 'medium' quality.
4) If all goes well the final lightmap will be produced in your specified directory, check this and try opening it in an image editor.

#### Showing the lightmap
If the shader on the level geometry is set up correctly, assigned to the correct lightmap file, then after baking the lightmap should automatically import and show up. If this doesn't happen, check the lightmap looks correct, then check the shaders on the meshes.

#### Ambient Occlusion
As well as lighting, LLightmap can bake ambient occlusion. You may decide to combine this with lights, or even use ambient occlusion on its own.

1) Change the bake_mode to 'AO'.
2) Hit the 'Bake Lightmap' button.
3) This may take some time to bake, once it is done, it should appear similarly to the lights.

#### Merging
Once you have made some intermediate lights and ambient occlusion exr files, you can merge them together.

1) Change the bake mode to 'Merge'.
2) Hit 'Bake Lightmap'.

Merging is much faster than baking lights, or AO. This is because baking requires many ray tests, whereas merging is simply mixing two textures together. When merging you can play with the overall brightness and gamma, and the balance between ambient occlusion and lights.

#### Final
* Once you are satisfied with your final lightmap (for a particular game level), make sure to save it somewhere.
* For your final game build, you don't need the LLightmap node, and should remove it (it won't be understood in standard builds of the engine anyway).
* You also don't need the intermediate lights and ao exr files. These should be deleted to reduce the export size of your game.
* Instead of using the final png, you can alternatively use an exr file for your lightmap (it will contain more accurate colors, with higher dynamic range). For mobile use png is still recommended.
* At the moment it is recommended to run noise reduction on the final lightmap with an image editing program (e.g. gimp, photoshop). Soon OpenImageDenoise will be in core and I will add a hook to use that as soon as available, which will make this more user friendly.

## Realtime lighting (Light probes)
See [LIGHT_PROBES.md](LIGHT_PROBES.md) for full information on using light probes

## Installation
* Get the latest godot engine 3.2 source, NOT 3.2 stable, it should be 3.2.2 or later as there has been a change to xatlas since then (commit 1bd5188 specifically, June 4th 2020):

Either:
1) go to https://github.com/godotengine/godot/tree/3.2 and click download zip for the bleeding edge latest
2) go to https://github.com/godotengine/godot/releases/tag/3.2.2-stable and click to download source code
3) A later official release of the source (e.g. 3.2.3 or 3.2.4 when available)

* Create a folder called 'llightmap' in the modules directory
* Clone / download this repository into the folder
* Compile as normal (instructions at godot website)
* The new node should be available called 'LLightmap'

As you only need it to create lightmaps on desktop as a preprocess, the lightmaps produced can be used with a standard build of the engine / standard templates.

I'm hoping to eventually make some builds for windows / linux x86_64 so users won't need to compile. The custom build is only needed to create the lightmaps / probe data. Once these are created they can be used in standard vanilla Godot engine (with a gdscript addon for probes, but that requires no compilation).

### Tips
* Set overall brightness with a combination of using a multiplier in the shader (see the example shader) and using `dynamic_range/normalize_bias` in the llightmap settings to make the lightmap overbright.
* For each light, you can control the volumetric effect (soft or hard shadows) by changing the x, y, and z scale in the node `Transform` properties. Note that due to a bug / feature in Godot, light scales get reset every time you move them in the editor, which can be annoying.
* Be sure to use bounces in order to get colors from textures. You can make a bounced scene darker by reducing the bounce power.
* For each light you can scale power with the `energy` parameter, and change color.
* For transparency, the material name should contain the string `_T_` (anywhere, at the end, start etc). And use a texture with alpha.
* In forward tracing you can scale the number of samples per light using the `indirect energy` light parameter. This is useful for directional lights which may need more samples.
* Spotlights have position, direction and spot angle, and volume with scale.
* When using spatial materials, the albedo texture will be found automatically. When using custom shaders, in order for LLightmap to find the texture colors for bouncing light, the uniform in the shader _must_ be called `texture_albedo`. Otherwise a plain white color will be used for bounces.
* Directional lights must point at least slightly downward. This isn't ideal but allows more consistent lighting, and the use case is mainly skies. For side lights, or lighting from below, area lights via omnis are a better bet.
