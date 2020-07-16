# godot-llightmap
Lightmap module for Godot Engine 3.x
Version 0.1

## Installation
* Get the latest godot engine 3.x source
* Create a folder called 'llightmap' in the modules directory
* Clone / download this repository into the folder
* Compile as normal (instructions at godot website)
* The new node should be available called 'LLightmap'

Due to using SSE2, this will probably only currently compile on x86_64, on linux and (hopefully) windows. I'll fix up the reference methods so it will compile on other platforms soon. However, as you only need it to create lightmaps on desktop as a preprocess, the lightmaps produced can be used with a standard build of the engine / standard templates.

<img src="/pics/llightmap_01screen.png" alt="screenshot" width="400"/>

