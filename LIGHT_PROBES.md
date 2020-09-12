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

## Baking Probes
* Simply create your lightmap as usual, making sure the output exr file is present for the lightmap (the AO is not used).
* Once you are happy with the lightmap, change the bake mode to `Light Probes`.
* Click `bake`.
* The data file (needed at runtime) will be output to the exact same folder / filename as your combined filename, except the extension will be changed to `.probe`. You might want to keep your combined lightmap (e.g. PNG) and the probe file together for each level.
* When baking probes, there are some settings in the `Light Probes` rollout in the inspector for the LLightmap node.
* Probe density determines how many probes are used. More probes gives more resolution, but larger data file sizes. You may need higher density on larger maps, but note that the number of probes is roughly the cube of this value, so data sizes / bake times may increase exponentially. 
* Probe samples. This determines how many samples are taken for indirect lighting, and against each light. The default should be fine, however you can reduce it for rough versions.

## Runtime
At runtime there are two main areas that need to be dealt with:

1) Loading and sampling the light data
2) Shaders that turn these samples into lighting that looks good

Ideally the light data would be sampled at runtime using fast c++ code, however this requires either a specially compiled version of the engine (and all templates for export platforms), or gdnative. For ease of distribution and compatibility I have therefore written an addon to load and sample the lightprobe data in gdscript. Later I will also make available a similar c++ module for those needing ultimate speed and prepared to compile templates.

### Using the gdscript LightProbes addon
The light probe functionality is contained within a single file, [lightprobes.gd](runtime/lightprobes.gd). You don't need to install this, just include it somewhere within your project.

Once `lightprobes.gd` is part of the project, you will be able to create a LightProbes object to handle the runtime tasks. You can consider `lightprobes.gd` to be a black box, and you only need to concern yourself with calling two functions:

* LightProbes::load_file(var filename)
* SampleResult LightProbes::sample(var position : Vector3)

You can see the SampleResult structure by looking in the `lightprobes.gd` file if you are interested. It is quite simple and contains light position, color, and indirect light color.

In order to call and make use of the SampleResult we need to consider some aspects of your game. Your game object has a position (xyz Vector3) and also should have a MeshInstance node which is used to display the object. These could be one and the same (a MeshInstance is derived from Spatial and also contains a position).

We need a position in order to know where on the map the lighting should be sampled. And we need the mesh instance in order to gain access to the `Material`. Materials are used in Godot to determine the color, texture and lighting of objects. We need to pass the lighting information from the SampleResult to the material being used by that particular object. We do this by using `shader parameters`. See the Godot docs for more info, but it should be quite easy to understand.

https://docs.godotengine.org/en/stable/classes/class_shadermaterial.html#class-shadermaterial-method-set-shader-param

You should pass the lighting info to the material every frame, using the `_process` function, then after the `_process` functions have run, the engine renders the objects. Now they have the correct lighting information available to the shader, the lighting should look correct.

e.g. 

Usually a player scene would take the form:
```
Spatial (player_node)
    MeshInstance (player_mesh)
```
or something similar. You might be using a skeleton to animate the mesh, but there should be a mesh instance somewhere (or possibly multiple, although this is not so recommended for performance reasons).

```
# your lightprobes object could e.g. be a global like this
var m_Probes : LightProbes = LightProbes.new()

var m_Player : Spatial
var m_PlayerMesh : MeshInstance

func _ready():
        # load a new probe file whenever you load a level
	# (doesn't have to be in a _ready function)
	m_Probes.load_file("res://Lightmaps/LightMap.probe")
	
	# the relative paths depend on your project and scene tree
	m_Player = $player_node
	m_PlayerMesh = $player_node/player_mesh

func _process(delta):
	# where the node holding my player is called player_node
	# and the MeshInstance of the player is called player_mesh 

	# where are we going to sample? roughly the middle of the dynamic object
	# so I add a bit of height for the player feet
	var sample_pos = m_Player.translation + Vector3(0, 0.5, 0)
	
	# make the sample - this returns a SampleResult structure
	var sample = m_Probes.sample(sample_pos)
	
	# get the material from our mesh, we need this to set the shader uniforms
	var mat : Material = m_PlayerMesh.get_surface_material(0)

	mat.set_shader_param("light_pos", sample.light_pos)
	mat.set_shader_param("light_color", sample.light_color)
	mat.set_shader_param("light_indirect", sample.light_color_indirect)
```

The above basically shows the workflow. You create a LightProbes object over the lifetime of your game, load the corresponding `probe` file when loading a level, then for each dynamic object you make a sample, then pass the information in the sample to the shader so it can render the lighting correctly.

You may need to make the material unique for each object so that the uniforms are unique for each object (I'm not sure on this yet, this is untested).

`player_node` and `player_mesh` refer to e.g. a spatial holding the object and the mesh instance. Or they could be the same. Make the best choice over whether to use global_transform or translation (which is a local translate) for your sample point - this will depend on your game.

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
