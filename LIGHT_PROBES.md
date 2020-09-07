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


## Appendix

### Medium Quality
* Vertex lighting, using faux roughness term from blue channel in fragment shader
* Suitable for skinned or non-animated models

```
shader_type spatial;

// we will be doing all the lighting ourself
render_mode unshaded;

// this is your opportunity to set the model texture
uniform sampler2D texture_albedo : hint_albedo;

// for lighting to work, you must update these uniforms each frame
// see the light probes documentation
uniform vec3 light_pos;
uniform vec4 light_color;
uniform vec4 light_indirect;

// specular is interpolated separately to diffuse + ambient in this shader
varying vec3 v_specular;

void vertex() {
	// Get the camera position in world space
	// by transforming a point at the origin by the camera matrix
	vec4 view_pos = vec4(0, 0, 0, 1);
	view_pos = CAMERA_MATRIX * view_pos;
	
	// absolute pain but there's a bug in the core shaders .. if skinning is applied
	// then using world_vertex_coords breaks the skinning. So we have to do the world
	// transform TWICE!! Once here, and once in the core shader.
	// Note that this step can potentially be avoided in non-skinned models
	// for better performance.
	
	// get the vertex in world space
	vec4 vert_world = WORLD_MATRIX * vec4(VERTEX, 1.0);
	
	// normal also needs to be transformed to world space
	vec3 normal = normalize((WORLD_MATRIX * vec4(NORMAL, 0.0)).xyz);
	
	// view direction
	vec3 vdir = vert_world.xyz - view_pos.xyz;
	vdir = normalize(vdir);
	
	// light direction
	vec3 ldir = (light_pos) - vert_world.xyz;

	// normalize light direction and get distance to light at same time
	float dist = length(ldir);

	// could divide by zero but I don't think this should cause GPU to freak out
	// GLES spec says divide by zero gives unspecified result, but mustn't interrupt shader
	ldir *= (1.0 / dist);
	
	// the surface normal dot the light direction gives us our diffuse lighting
	float d = dot(normal, ldir);
	
	// use dot product to calculate reflection direction of the light ray on the surface
	vec3 rdir = ldir - (d * 2.0 * normal);
	
	// specular light blob depends on dot reflection direction to view dir
	float dot_refl = dot(vdir, rdir);
	
	// don't want negative specular values
	float specular = max(0.0, dot_refl);

	// optional, apply a power operator to the specular highlight to get a finer blob	
	//spec = pow(spec, 4);
	
	// diffuse light, should be none from behind the surface
	d = max(0.0, d);
	
	// magic falloff, prevents divide by zero and more linear that just inverse square
	dist += 10.0;
	float falloff = 1.0 / (0.01 * (dist * dist));
	falloff = max (0.0, falloff);
	
	d *= falloff;
	specular *= falloff;
	
	// scale diffuse
	d *= 0.8;
	
	// to test just the indirect lighting, set these both to zero
	// d = 0.0;
	// specular = 0.0f;
	
	// the diffuse + ambient color.
	// we are letting the indirect light have some influence on the diffuse light here -
	// this isn't necessary but is artistic licence.
	COLOR = vec4(((light_color.rgb + light_indirect.rgb) * d) + light_indirect.rgb, 1);
	
	// the specular depends on the angles and the light color
	v_specular = specular * light_color.rgb;
}

void fragment() {
	// albedo color
	vec4 albedo_tex = texture(texture_albedo,UV);
	
	// whole point of sending specular separately is we can apply some
	// 'faux' roughness. We are just using part of the albedo for roughness,
	// however this could be sent as a separate channel (e.g. alpha), however
	// alpha is not present in some compressed formats.
	vec3 spec = v_specular * (1.0 - albedo_tex.b);
	
	// the separate specular could be avoided on low power mobile, as the shader
	// will run faster without roughness.
	
	// (diffuse + ambient + specular) * the albedo	
	ALBEDO = (COLOR.rgb + spec) * albedo_tex.rgb;
}
```
