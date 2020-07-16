# godot-llightmap
Lightmap module for Godot Engine 3.x
Version 0.1
<img src="/pics/llightmap_01screen.png" alt="screenshot" width="400"/>

## About
While JFons is finishing off the new official core lightmapper for 3.2, I've spent some time playing making this toy lightmapper. It works so well I'm now going to be using it to replace the lightmapping functionality in my other module, lportal.

LLightmap is designed to be easy to use and produce nice soft shadows, with high performance at runtime. The existing core lightmappers produce a different lightmap for each object. This can be slow due to texture swapping. LLightmap instead can uv map your scene with multiple objects, unwrapping so that all your objects share the same texture space. Thus a single lightmap can be used.

## Features
* Forward and backward ray tracing options
* Ambient occlusion
* Linear HDR exr intermediate files
* Linear HDR exr final texture or gamma corrected normalized png

## Installation
* Get the latest godot engine 3.x source
* Create a folder called 'llightmap' in the modules directory
* Clone / download this repository into the folder
* Compile as normal (instructions at godot website)
* The new node should be available called 'LLightmap'

Due to using SSE2, this will probably only currently compile on x86_64, on linux and (hopefully) windows. I'll fix up the reference methods so it will compile on other platforms soon. However, as you only need it to create lightmaps on desktop as a preprocess, the lightmaps produced can be used with a standard build of the engine / standard templates.


