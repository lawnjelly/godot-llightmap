#pragma once

namespace LM {


class SceneSaver
{
public:
	bool SaveScene(Node * pNode, String szFilename, bool reset_filenames = false);
	void SetOwnerRecursive(Node * pNode, Node * pOwner);
	void SetFilenameRecursive(Node * pNode);
};

} // namespace
