#include "lmerger.h"
#include "llighttypes.h"

#ifdef TOOLS_ENABLED
#include "thirdparty/xatlas/xatlas.h"

extern bool (*array_mesh_lightmap_unwrap_callback)(float p_texel_size, const float *p_vertices, const float *p_normals, int p_vertex_count, const int *p_indices, const int *p_face_materials, int p_index_count, float **r_uv, int **r_vertex, int *r_vertex_count, int **r_index, int *r_index_count, int *r_size_hint_x, int *r_size_hint_y);
#endif

namespace LM {


int Merger::m_iUVPadding = 4;

Node * Merger::FindSceneRoot(Node * pNode) const
{
	if (pNode->get_parent())
	{
		return FindSceneRoot(pNode->get_parent());
	}

	return pNode;
}

MeshInstance * Merger::Merge(Spatial * pRoot, int padding)
{
	m_iUVPadding = padding;

	FindMeshes(pRoot);

	MeshInstance * pMerged = memnew(MeshInstance);
	pMerged->set_name("lightmap_proxy");

	Spatial * pParent = pRoot->get_parent_spatial();
	if (!pParent)
		return 0;

//	Node * pSceneRoot = FindSceneRoot(pParent);
//	pMerged->set_owner(pSceneRoot);

	pParent->add_child(pMerged);


//	print_line(pParent->get_name());
//	int nChildren = pParent->get_child_count();
//	for (int n=0; n<nChildren; n++)
//	{
//		Node * pChild = pParent->get_child(n);
//		print_line("\t" + pChild->get_name());
//	}

	if (!MergeMeshes(*pMerged))
	{
		pMerged->queue_delete();
		return 0;
	}
//	return 0;

	return pMerged;
}


bool Merger::MergeMeshes(MeshInstance &merged)
{
	PoolVector<Vector3> verts;
	PoolVector<Vector3> normals;
	PoolVector<int> inds;


	for (int n=0; n<m_Meshes.size(); n++)
	{
		Merge_MeshInstance(*m_Meshes[n], verts, normals, inds);
	}

	print_line("Merging, num verts is " + itos(verts.size()));

	// if there are no verts, there has been an error, abort
	if (!verts.size())
	{
		WARN_PRINT_ONCE("Merger::MergeMeshes : No vertices, aborting");
		OS::get_singleton()->alert("Error: No vertices, aborting. Check output log for details.", "MergeMeshes");
		return false;
	}

	// lightmap unwrap
	//LightmapUnwrap(verts, normals, inds, uv2s);

	// unmerge the original meshes (write the uvs2 back)


	Ref<ArrayMesh> am;
	am.instance();
	Array arr;
	arr.resize(Mesh::ARRAY_MAX);
	arr[Mesh::ARRAY_VERTEX] = verts;
	arr[Mesh::ARRAY_NORMAL] = normals;
	arr[Mesh::ARRAY_INDEX] = inds;
//	arr[Mesh::ARRAY_TEX_UV2] = uv2s;

	// bug fix. The lightmapping code removes duplicate vertices, which we DONT want
	// as it makes the merged mesh get out of sync with the original meshes.
	// To prevent this we will create a dummy set of UV1s that are unique.
//	PoolVector<Vector2> dummy_uvs;
//	for (int n=0; n<verts.size(); n++)
//	{
//		dummy_uvs.append(Vector2(n, n));
//	}
//	arr[Mesh::ARRAY_TEX_UV] = dummy_uvs;

	// sanity check on the indices
	for (int n=0; n<inds.size(); n++)
	{
		int i = inds[n];
		assert (i < verts.size());
	}


//	void add_surface_from_arrays(PrimitiveType p_primitive, const Array &p_arrays, const Array &p_blend_shapes = Array(), uint32_t p_flags = ARRAY_COMPRESS_DEFAULT);

	am->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, arr, Array(), Mesh::ARRAY_COMPRESS_DEFAULT);

	//if (bLightmapUnwrap)
	LightmapUnwrap(am, merged.get_global_transform());

	// duplicate the UV2 to uv1 just in case they are needed
	arr[Mesh::ARRAY_TEX_UV] = arr[Mesh::ARRAY_TEX_UV2];

	merged.set_mesh(am);

	// set mesh to use in baked lighting
	merged.set_flag(GeometryInstance::FLAG_USE_BAKED_LIGHT, true);

	return true;
}

void Merger::Merge_MeshInstance(const MeshInstance &mi, PoolVector<Vector3> &verts, PoolVector<Vector3> &norms, PoolVector<int> &inds)
{
	// some godot jiggery pokery to get the mesh verts in local space
	Ref<Mesh> rmesh = mi.get_mesh();

	if (rmesh->get_surface_count() == 0)
		return;

	Array arrays = rmesh->surface_get_arrays(0);
	PoolVector<Vector3> p_vertices = arrays[VS::ARRAY_VERTEX];
	PoolVector<Vector3> p_normals = arrays[VS::ARRAY_NORMAL];
	PoolVector<int> p_indices = arrays[VS::ARRAY_INDEX];
	//PoolVector<int>::Read ir = mesh_indices.read();

	// NEW .. the checking for valid triangles should be on WORLD SPACE vertices,
	// NOT model space

	// special case, if no indices, create some
	int num_indices_before = p_indices.size();
	if (!EnsureIndicesValid(p_indices, p_vertices))
	{
		print_line("\tignoring INVALID TRIANGLES (duplicate indices or zero area triangle) detected in " + mi.get_name() + ", num inds before / after " + itos (num_indices_before) + " / " + itos(p_indices.size()));
	}

	// no normals, can't do
	if (!p_normals.size())
	{
		//assert (0);
		WARN_PRINT_ONCE("Merge_MI : mesh with no normals, ignoring");
		return;
	}


	// the first index of this mesh is offset from the verts we already have stored in the merged mesh
	int first_index = verts.size();

//	print_line("Merge MI : " + mi.get_name() + "\tFirstVert : " + itos(first_index) + "\tNumUVs : " + itos(p_vertices.size()));

	// transform verts to world space
	Transform trans = mi.get_global_transform();

	for (int n=0; n<p_vertices.size(); n++)
	{
		Vector3 ptWorld = trans.xform(p_vertices[n]);

		// hacky way for normals, we should use transpose of inverse matrix, dunno if godot supports this
		Vector3 ptNormA = Vector3(0, 0, 0);
		Vector3 ptNormWorldA = trans.xform(ptNormA);
		Vector3 ptNormWorldB = trans.xform(p_normals[n]);

		Vector3 ptNorm = ptNormWorldB - ptNormWorldA;

		verts.push_back(ptWorld);
		norms.push_back(ptNorm);
	}

	// indices
	for (int n=0; n<p_indices.size(); n++)
		inds.push_back(p_indices[n] + first_index);


	// new .. test for zero size tris.
//	int nTris = p_indices.size() / 3;
//	int indCount = 0;
//	for (int t=0; t<nTris; t++)
//	{
//		const Vector3 &v0 = verts[p_indices[indCount++]];
//		const Vector3 &v1 = verts[p_indices[indCount++]];
//		const Vector3 &v2 = verts[p_indices[indCount++]];

//		// edge lengths
//		float l0 = (v1 - v0).length();
//		float l1 = (v2 - v1).length();
//		float l2 = (v0 - v2).length();

//		const float epsilon = 0.0001f;
//		assert (l0 > epsilon);
//		assert (l1 > epsilon);
//		assert (l2 > epsilon);
//	}

}

void Merger::FindMeshes(Spatial * pNode)
{
	// mesh instance?
	MeshInstance * pMI = Object::cast_to<MeshInstance>(pNode);
	if (pMI)
	{
		if (IsMeshInstanceSuitable(*pMI))
		{
			print_line("found mesh : " + pMI->get_name());
			m_Meshes.push_back(pMI);
		}
	}

	for (int n=0; n<pNode->get_child_count(); n++)
	{
		Spatial * pChild = Object::cast_to<Spatial>(pNode->get_child(n));
		if (pChild)
		{
			FindMeshes(pChild);
		}
	}
}

bool Merger::LightmapUnwrap(Ref<ArrayMesh> am, const Transform &trans)
{
#ifdef TOOLS_ENABLED
	array_mesh_lightmap_unwrap_callback = xatlas_unwrap;

	// we can add the UV2 coords from here
	Error err = am->lightmap_unwrap(trans);
	if (err != OK) {
		print_line ("UV Unwrap failed, mesh may not be manifold?");
		return false;
	}
#endif
	return true;
}

#ifdef TOOLS_ENABLED
bool Merger::xatlas_unwrap(float p_texel_size, const float *p_vertices, const float *p_normals, int p_vertex_count, const int *p_indices, const int *p_face_materials, int p_index_count, float **r_uv, int **r_vertex, int *r_vertex_count, int **r_index, int *r_index_count, int *r_size_hint_x, int *r_size_hint_y)
{

	//set up input mesh
	xatlas::MeshDecl input_mesh;
	input_mesh.indexData = p_indices;
	input_mesh.indexCount = p_index_count;
	input_mesh.indexFormat = xatlas::IndexFormat::UInt32;

	input_mesh.vertexCount = p_vertex_count;
	input_mesh.vertexPositionData = p_vertices;
	input_mesh.vertexPositionStride = sizeof(float) * 3;
	input_mesh.vertexNormalData = p_normals;
	input_mesh.vertexNormalStride = sizeof(uint32_t) * 3;
	input_mesh.vertexUvData = NULL;
	input_mesh.vertexUvStride = 0;

	xatlas::ChartOptions chart_options;
	xatlas::PackOptions pack_options;

	pack_options.maxChartSize = 4096;
	pack_options.blockAlign = true;
	pack_options.texelsPerUnit = 1.0 / p_texel_size;
	pack_options.padding = m_iUVPadding; // 4

	xatlas::Atlas *atlas = xatlas::Create();
	printf("Adding mesh..\n");
	xatlas::AddMeshError::Enum err = xatlas::AddMesh(atlas, input_mesh, 1);
	ERR_FAIL_COND_V_MSG(err != xatlas::AddMeshError::Enum::Success, false, xatlas::StringForEnum(err));

	printf("Generate..\n");
	xatlas::Generate(atlas, chart_options, xatlas::ParameterizeOptions(), pack_options);

	*r_size_hint_x = atlas->width;
	*r_size_hint_y = atlas->height;

	float w = *r_size_hint_x;
	float h = *r_size_hint_y;

	if (w == 0 || h == 0) {
		return false; //could not bake because there is no area
	}

	const xatlas::Mesh &output = atlas->meshes[0];

	*r_vertex = (int *)malloc(sizeof(int) * output.vertexCount);
	*r_uv = (float *)malloc(sizeof(float) * output.vertexCount * 2);
	*r_index = (int *)malloc(sizeof(int) * output.indexCount);

	float max_x = 0;
	float max_y = 0;
	for (uint32_t i = 0; i < output.vertexCount; i++) {
		(*r_vertex)[i] = output.vertexArray[i].xref;
		(*r_uv)[i * 2 + 0] = output.vertexArray[i].uv[0] / w;
		(*r_uv)[i * 2 + 1] = output.vertexArray[i].uv[1] / h;
		max_x = MAX(max_x, output.vertexArray[i].uv[0]);
		max_y = MAX(max_y, output.vertexArray[i].uv[1]);
	}

	printf("Final texture size: %f,%f - max %f,%f\n", w, h, max_x, max_y);
	*r_vertex_count = output.vertexCount;

	for (uint32_t i = 0; i < output.indexCount; i++) {
		(*r_index)[i] = output.indexArray[i];
	}

	*r_index_count = output.indexCount;

	xatlas::Destroy(atlas);
	printf("Done\n");
	return true;
}
#endif // TOOLS_ENABLED

} // namespace
