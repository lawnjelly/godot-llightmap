/* register_types.cpp */

#include "register_types.h"

#include "core/class_db.h"
#include "gdlightmapper.h"
#include "llightmapper_editor_plugin.h"


void register_llightmap_types() {

	ClassDB::register_class<LLightmap>();

#ifdef TOOLS_ENABLED
    EditorPlugins::add_by_type<LLightmapEditorPlugin>();
#endif
}

void unregister_llightmap_types() {
   //nothing to do here
}
