extends Reference

class_name LightProbes

###############################################################
# CLIENT INTERFACE

# the user can set this if desired for a different look
var m_Directional_Distance : float = 15.0

# samples are returned with this structure
class SampleResult:
	var light_pos : Vector3
	var light_color : Color
	var light_color_indirect : Color

func load_file(var szFilename):
	return _load_file(szFilename) 

# main client function to sample the light data
# returns a SampleResult
func sample(var pos : Vector3):
	return _sample(pos)

###############################################################
# INTERNAL

var _m_File : File
var _m_ProbeMap : _ProbeMap
var _m_bError : bool = false
var _m_bDebugFrame : bool = false


func _sample(var pos : Vector3):
	var orig_pos = pos
	
	# get pos in voxel space
	pos -= _m_ProbeMap.ptMin
	pos /= _m_ProbeMap.voxel_size

	# within?
	var pt : _Vec3i = _Vec3i.new()
	pt.from(pos)
	
	# fractions through the voxel
	var fx = pos.x - pt.x
	var fy = pos.y - pt.y
	var fz = pos.z - pt.z
	var inv_x = 1.0 - fx
	var inv_y = 1.0 - fy
	var inv_z = 1.0 - fz
	
	pt = _clamp_to_map(pt)

	var sr : SampleResult = SampleResult.new()
	
	# octaprobe
	var opr : _Octaprobe = _get_octaprobe(pt)

	# indirect light
	sr.light_color_indirect.r = _trilinear(fx, fy, fz, opr.indirect_r)
	sr.light_color_indirect.g = _trilinear(fx, fy, fz, opr.indirect_g)
	sr.light_color_indirect.b = _trilinear(fx, fy, fz, opr.indirect_b)
	sr.light_color_indirect.a = 1.0

	var total_influence = 0.0
	var total_pos = Vector3()
	var total_power = 0.0
	var total_col = Color()
	var max_power = 0.0

	# trilinear interpolation
	for l in range (opr.lights.size()):
		var ol : _Octalight = opr.lights[l]
		
		# find interpolated power for this light
		var power = _trilinear(fx, fy, fz, ol.powers)
			
		# get the light info
		var pl : _ProbeLight = _m_ProbeMap.lights[ol.light_id]
		
		# special cases for types of light
		# directional
		if pl.type == 2:
			# always make the source pos of a directional light parallel to the rays
			pl.pos = orig_pos - (pl.dir * m_Directional_Distance)
		
		# ray offset
		var offset : Vector3 = pl.pos - orig_pos

		# +1 prevents divide by zero
		var dist = offset.length_squared() + 1.0
		
		# the influence of each light determines how they mix
		var influence = (1.0 / dist) * power
		
		
		total_influence += influence;
		total_pos += pl.pos * influence
		total_col += pl.color * influence
		total_power += power

		# use the maximum light power rather than adding the lights
		# to prevent overbright
		if (power > max_power):
			max_power = power
		
#		if m_bDebugFrame:
#			print ("light pos " + str(pl.pos) + " our pos " + str(orig_pos) + " dist " + str(dist) + " power " + str(power))
	
	if (total_influence > 0.0):
		# divide by zero? - maybe influence can never be zero
		sr.light_color = total_col / total_influence
		sr.light_pos = total_pos / total_influence
		
		# apply power to color
		sr.light_color *= max_power

	# sr should be auto blanked if no winner
		
	return sr


class _Vec3i:
	func Set(var xx, var yy, var zz):
		x = xx
		y = yy
		z = zz
	func from(var pt : Vector3):
		x = pt.x
		y = pt.y
		z = pt.z
	func to_string()->String:
		return str(x) + ", "+ str(y) + ", " + str(z)
		
	var x : int
	var y : int
	var z : int

class _Octalight:
	var light_id : int
	var powers = [] # 8 powers per octalight

class _Octaprobe:
	var indirect_r = []
	var indirect_g = []
	var indirect_b = []
	var lights = []

class _Probe:
	var col_indirect : Color = Color()
	var contributions = []

class _Contribution:
	var light_id : int
	var power : float

class _ProbeMap:
	var voxel_size : Vector3
	var ptMin
	var dims : _Vec3i = _Vec3i.new()
	var XTimesY : int
	var probes = []
	var lights = []
	var octaprobes = []

class _ProbeLight:
	var type : int
	var pos : Vector3
	var dir : Vector3
	var energy : float
	var rang : float
	var color : Color
	var spot_angle_radians : float


func _create_octaprobe(var x, var y, var z):
	var prs = []
	prs.push_back(_get_probe_xyz(x, y, z))
	prs.push_back(_get_probe_xyz(x, y, z+1))
	prs.push_back(_get_probe_xyz(x+1, y, z+1))
	prs.push_back(_get_probe_xyz(x+1, y, z))
	prs.push_back(_get_probe_xyz(x, y+1, z))
	prs.push_back(_get_probe_xyz(x, y+1, z+1))
	prs.push_back(_get_probe_xyz(x+1, y+1, z+1))
	prs.push_back(_get_probe_xyz(x+1, y+1, z))
	
	var octaprobe : _Octaprobe = _Octaprobe.new()
	
	# first step fill the lights
	for p in range (8):
		var pr : _Probe = prs[p]

		# indirect light
		octaprobe.indirect_r.push_back(pr.col_indirect.r)
		octaprobe.indirect_g.push_back(pr.col_indirect.g)
		octaprobe.indirect_b.push_back(pr.col_indirect.b)
		
		# go through contributions
		for n in range (pr.contributions.size()):
			var cont : _Contribution = pr.contributions[n]
			
			# is the light id already in the list? if not add it
			var exists = false
	
			for ol in range (octaprobe.lights.size()):
				if (octaprobe.lights[ol].light_id == cont.light_id):
					exists = true
					break
			
			if exists == false:
				var octalight = _Octalight.new()
				octalight.light_id = cont.light_id
				octaprobe.lights.push_back(octalight)

	# second step add the values for each corner
	for ol_id in range (octaprobe.lights.size()):
		var ol : _Octalight = octaprobe.lights[ol_id]
		var light_id = ol.light_id
		
		for pr_count in range (8):
			var pr : _Probe = prs[pr_count]
			
			var power : float = 0.0
			
			# go through contributions
			for n in range (pr.contributions.size()):
				var cont : _Contribution = pr.contributions[n]
				
				if cont.light_id == light_id:
					power = 1.0# cont.power
					break
		
			# add the power to that octalight (either specified in the probe, or 0.0 if not specified)
			ol.powers.push_back(power)
			
	# finally add the octaprobe
	_m_ProbeMap.octaprobes.push_back(octaprobe)

	pass

func _create_octaprobes():
	for z in range (_m_ProbeMap.dims.z):
		for y in range (_m_ProbeMap.dims.y):
			for x in range (_m_ProbeMap.dims.x):
				_create_octaprobe(x, y, z)

func _test_interpolate():
	var r
	r = _bilinear(0, 0, 0, 1, 0, 0)
	r = _bilinear(1, 0, 0, 1, 0, 0)
	r = _bilinear(0, 1, 0, 1, 0, 0)
	r = _bilinear(1, 1, 0, 1, 0, 0)


func _trilinear(var fx, var fy, var fz, var samples)->float:
	var c000 = samples[0]
	var c001 = samples[1]
	var c101 = samples[2]
	var c100 = samples[3]
	var c010 = samples[4]
	var c011 = samples[5]
	var c111 = samples[6]
	var c110 = samples[7]
	var e = _bilinear(fx, fy, c000, c100, c010, c110)
	var f = _bilinear(fx, fy, c001, c101, c011, c111)
	return e * (1 - fz) + (f * fz)

func _bilinear(var fx : float, var fy : float, var c00 : float, var c10 : float, var c01 : float, var c11 : float)->float:
	var a : float = c00 * (1.0 - fx) + c10 * fx
	var b : float = c01 * (1.0 - fx) + c11 * fx
	
	return a * (1.0 - fy) + (b * fy)


func _get_probe_xyz(var x, var y, var z)->_Probe:
	var pt : _Vec3i = _Vec3i.new()
	pt.Set(x, y, z)
	return _get_probe(pt)

func _get_octaprobe_xyz(var x, var y, var z)->_Octaprobe:
	var pt : _Vec3i = _Vec3i.new()
	pt.Set(x, y, z)
	pt = _limit_to_map(pt)

	return _get_octaprobe(pt)
	
func _get_octaprobe(var pt : _Vec3i)->_Octaprobe:
	var i : int = pt.z * _m_ProbeMap.XTimesY
	i += pt.y * _m_ProbeMap.dims.x
	i += pt.x
	
	assert (i < _m_ProbeMap.probes.size())
	return _m_ProbeMap.octaprobes[i]
	

func _limit_to_map(var pt : _Vec3i)->_Vec3i:
	# cap to map
	if (pt.x < 0):
		pt.x = 0
	if (pt.x >= _m_ProbeMap.dims.x):
		pt.x = _m_ProbeMap.dims.x-1
	if (pt.y < 0):
		pt.y = 0
	if (pt.y >= _m_ProbeMap.dims.y):
		pt.y = _m_ProbeMap.dims.y-1
	if (pt.z < 0):
		pt.z = 0
	if (pt.z >= _m_ProbeMap.dims.z):
		pt.z = _m_ProbeMap.dims.z-1
		
	return pt
	

func _get_probe(var pt : _Vec3i)->_Probe:
	pt = _limit_to_map(pt)
	
	var i : int = pt.z * _m_ProbeMap.XTimesY
	i += pt.y * _m_ProbeMap.dims.x
	i += pt.x
	
	assert (i < _m_ProbeMap.probes.size())
	return _m_ProbeMap.probes[i]
	

func _clampi(var i : int, var mn : int, var mx : int)->int:
	if (i < mn):
		i = mn
	if (i > mx):
		i = mx
	return i

func _clamp_to_map(var pt : _Vec3i)->_Vec3i:
	pt.x = _clampi(pt.x, 0, _m_ProbeMap.dims.x-1)
	pt.y = _clampi(pt.y, 0, _m_ProbeMap.dims.y-1)
	pt.z = _clampi(pt.z, 0, _m_ProbeMap.dims.z-1)
	return pt

func _read_vec3()->Vector3:
	var v = Vector3()
	v.x = _m_File.get_float()
	v.y = _m_File.get_float()
	v.z = _m_File.get_float()
	return v


func _load_lights():
	var nLights = _m_File.get_16()
	
	for n in range (nLights):
		var l : _ProbeLight = _ProbeLight.new()
		
		l.type = _m_File.get_8()
		
		l.pos = _read_vec3()
		l.dir = _read_vec3()
		l.energy = _m_File.get_float()
		l.rang = _m_File.get_float()
		
		l.color.r = _m_File.get_float()
		l.color.g = _m_File.get_float()
		l.color.b = _m_File.get_float()
		
		l.spot_angle_radians = _m_File.get_float()
		
		_m_ProbeMap.lights.push_back(l)
	
	pass



func _load_probeA(var x, var y, var z):
	var p : _Probe = _Probe.new()
	
	var nContribs = _m_File.get_8()
	
	for n in range (nContribs):
		var c : _Contribution = _Contribution.new()
		c.light_id = _m_File.get_8()
		
		p.contributions.push_back(c)


	_m_ProbeMap.probes.push_back(p)
	pass

func _load_probeB(var count : int):
	var p : _Probe = _m_ProbeMap.probes[count]
	
	# color
	var r : float = float (_m_File.get_8())
	r /= 127.0
	var g : float = float (_m_File.get_8())
	g /= 127.0
	var b : float = float (_m_File.get_8())
	b /= 127.0
	
	p.col_indirect = Color(r, g, b, 1.0)
	
	var nContribs = p.contributions.size()
	
	for n in range (nContribs):
		var power : float = float (_m_File.get_8())
		power /= 127.0
		p.contributions[n].power = power
	
	# restore
	_m_ProbeMap.probes[count] = p
	pass

func _load_probes():
	for z in range (_m_ProbeMap.dims.z):
		for y in range (_m_ProbeMap.dims.y):
			for x in range (_m_ProbeMap.dims.x):
				_load_probeA(x, y, z)
				
	var count : int = 0
	for z in range (_m_ProbeMap.dims.z):
		for y in range (_m_ProbeMap.dims.y):
			for x in range (_m_ProbeMap.dims.x):
				_load_probeB(count)
				count += 1
		
	pass

func _load_file(var szFilename):
	
	# create new probemap
	_m_ProbeMap = null
	_m_ProbeMap = _ProbeMap.new()
	
	_m_File = File.new()
	
	var err = _m_File.open(szFilename, File.READ)
	if err != OK:
		return false
	
	# load fourcc
	var fourcc_matches = 0
	
	# must be Prob (in ascii)
	var c0 = _m_File.get_8()
	if (c0 == 80):
		fourcc_matches += 1
	var c1 = _m_File.get_8()
	if (c1 == 114):
		fourcc_matches += 1
	var c2 = _m_File.get_8()
	if (c2 == 111):
		fourcc_matches += 1
	var c3 = _m_File.get_8()
	if (c3 == 98):
		fourcc_matches += 1
	
	if (fourcc_matches != 4):
		OS.alert("Error", "Not a probe file")
		return false
	
	# version
	var version = _m_File.get_16()
	if (version != 100):
		OS.alert("Probe file wrong version", "Re-export with matching version")
		return false
	
	_m_ProbeMap.dims.x = _m_File.get_16()
	_m_ProbeMap.dims.y = _m_File.get_16()
	_m_ProbeMap.dims.z = _m_File.get_16()
	
	_m_ProbeMap.XTimesY = _m_ProbeMap.dims.x * _m_ProbeMap.dims.y
	
	_m_ProbeMap.ptMin = _read_vec3()
	_m_ProbeMap.voxel_size = _read_vec3()
	
	_load_lights()
	
	_load_probes()
	
	_m_File.close()
	
	_create_octaprobes()
	return true
