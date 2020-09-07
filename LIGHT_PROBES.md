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

