# Light Probes

## About
Although lightmaps are great for producing high quality static lighting, they won't provide realtime lighting on dynamic objects. There are various approaches to combine lightmaps with realtime lighting. We could simply include our lights and have them only affect our dynamic objects. This is fairly easy to do - the example lightmap shaders are marked as `unshaded`, so will not be affected by realtime lights.

There are a few problems with using realtime lights, mostly related to performance.

1) The lights won't cast shadows from the environment, unless we allow the background to render into the shadow map. Otherwise for example lights might affect objects on the other side of walls.

2) Consider a level with 30 lights. How do we know which ones should be affecting our objects? We might have to render large numbers of lights, which is very expensive.

3) Realtime lights take no account of indirect light.

An established solution to deal with this problem is the use of 'light probes'.

We place many probes throughout our level, and use them to sample the local lighting conditions in a similar way to generating the lightmap. Then at runtime, we load this probe data, find the nearest probes to our sample points (e.g. players), and interpolate the probe data to give an approximation of local lighting.

This approach can work very well, although it doesn't deal with dynamic shadows.

### Pros
* Can be a high / best performing solution for realtime lighting
* Can use indirect light (bounces, emission)

### Cons
* Data files can potentially be huge
* Probes are static, like the lightmaps
* No dynamic shadows (but can help with schemes for dynamic shadows)
* Require some thought to shaders to make use of them

## Priorities
There are a number of tradeoffs involved in the size of the lightprobe data, the information stored, the speed at runtime and the final quality of the result.

Unlike much of Godot, LLightmap is designed from the get go for high performance on low powered hardware, especially mobiles. As such the primary lightprobes are designed to take up very little data (for small download sizes), and to run as fast as possible. It is possible that I may add an optional higher quality lightprobe system at a later date.

My aim is to produce level probe data files that when compressed come to around 50kb or less. This should allow a typical 20 level game to use 1 meg for light probe data.

## Runtime
At runtime there are two main areas that need to be dealt with:

1) Loading and sampling the light data
2) Shaders that turn these samples into lighting that looks good

Ideally the light data would be sampled at runtime using fast c++ code, however this requires either a specially compiled version of the engine (and all templates for export platforms), or gdnative. For ease of distribution and compatibility I have therefore written an addon to load and sample the lightprobe data in gdscript. Later I will also make available a similar c++ module for those needing ultimate speed and prepared to compile templates.

### Shaders
Although the shaders for the lightmaps are usually quite simple, the shaders for dynamic objects are a bit more complex. This is because they need to bypass the existing Godot lighting, and we do the lighting ourselves in the shader.

The process is as follows, for an example player:

1) Sample the lightdata at a point roughly in the centre of the player. This will provide us:

* Ambient indirect light color
* Primary light position
* Primary light color

These need to be passed to the material shader as uniforms:

```
	var sample = m_Probes.sample(sample_pos)

	var mat : Material = mesh_instance.get_surface_material(0)

	mat.set_shader_param("light_pos", sample.light_pos)
	mat.set_shader_param("light_color", sample.light_color)
	mat.set_shader_param("light_indirect", sample.light_indirect)
```

2) The shader performs the magic, by calculating diffuse, ambient and specular lighting from this sample.

I will provide some example shaders, which can be used 'as is', but feel free to modify these to get the effect you are after in your particular game.





