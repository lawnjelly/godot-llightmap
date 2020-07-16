#include "llightmapper_editor_plugin.h"

//void LLightmapEditorPlugin::_uvmap()
//{
//	if (lightmap)
//	{
//		lightmap->uvmap();
//	}
//}

void LLightmapEditorPlugin::_bake() {

	if (lightmap)
	{
		lightmap->lightmap_bake();
	}
	/*
	if (lightmap) {
		BakedLightmap::BakeError err;
		if (get_tree()->get_edited_scene_root() && get_tree()->get_edited_scene_root() == lightmap) {
			err = lightmap->bake(lightmap);
		} else {
			err = lightmap->bake(lightmap->get_parent());
		}

		switch (err) {
			case BakedLightmap::BAKE_ERROR_NO_SAVE_PATH:
				EditorNode::get_singleton()->show_warning(TTR("Can't determine a save path for lightmap images.\nSave your scene (for images to be saved in the same dir), or pick a save path from the BakedLightmap properties."));
				break;
			case BakedLightmap::BAKE_ERROR_NO_MESHES:
				EditorNode::get_singleton()->show_warning(TTR("No meshes to bake. Make sure they contain an UV2 channel and that the 'Bake Light' flag is on."));
				break;
			case BakedLightmap::BAKE_ERROR_CANT_CREATE_IMAGE:
				EditorNode::get_singleton()->show_warning(TTR("Failed creating lightmap images, make sure path is writable."));
				break;
			default: {
			}
		}
	}
	*/
}

void LLightmapEditorPlugin::edit(Object *p_object) {

	LLightmap *s = Object::cast_to<LLightmap>(p_object);
	if (!s)
		return;

	lightmap = s;
}

bool LLightmapEditorPlugin::handles(Object *p_object) const {

	return p_object->is_class("LLightmap");
}

void LLightmapEditorPlugin::make_visible(bool p_visible) {

	if (p_visible) {
		bake->show();
		//uvmap->show();
	} else {

		bake->hide();
		//uvmap->hide();
	}
}

EditorProgress *LLightmapEditorPlugin::tmp_progress = NULL;

void LLightmapEditorPlugin::bake_func_begin(int p_steps) {

	ERR_FAIL_COND(tmp_progress != NULL);

	tmp_progress = memnew(EditorProgress("bake_lightmaps", TTR("Bake Lightmaps"), p_steps, true));
}

bool LLightmapEditorPlugin::bake_func_step(int p_step, const String &p_description) {

	ERR_FAIL_COND_V(tmp_progress == NULL, false);
	return tmp_progress->step(p_description, p_step, false);
}

void LLightmapEditorPlugin::bake_func_end() {
	ERR_FAIL_COND(tmp_progress == NULL);
	memdelete(tmp_progress);
	tmp_progress = NULL;
}

void LLightmapEditorPlugin::_bind_methods() {

	ClassDB::bind_method("_bake", &LLightmapEditorPlugin::_bake);
	//ClassDB::bind_method("_uvmap", &LLightmapEditorPlugin::_uvmap);
}

LLightmapEditorPlugin::LLightmapEditorPlugin(EditorNode *p_node) {

	editor = p_node;
	bake = memnew(ToolButton);
	bake->set_icon(editor->get_gui_base()->get_icon("Bake", "EditorIcons"));
	bake->set_text(TTR("Bake Lightmap"));
	bake->hide();
	bake->connect("pressed", this, "_bake");
	add_control_to_container(CONTAINER_SPATIAL_EDITOR_MENU, bake);
	lightmap = NULL;

//	uvmap = memnew(ToolButton);
//	uvmap->set_icon(editor->get_gui_base()->get_icon("Bake", "EditorIcons"));
//	uvmap->set_text(TTR("UV Map"));
//	uvmap->hide();
//	uvmap->connect("pressed", this, "_uvmap");
//	add_control_to_container(CONTAINER_SPATIAL_EDITOR_MENU, uvmap);


	LightMapper::bake_begin_function = bake_func_begin;
	LightMapper::bake_step_function = bake_func_step;
	LightMapper::bake_end_function = bake_func_end;
}

LLightmapEditorPlugin::~LLightmapEditorPlugin() {
}
