#pragma once

namespace LM {


class SceneSaver
{
public:
	bool SaveScene(Node * pNode, String szFilename);
	void SetOwnerRecursive(Node * pNode, Node * pOwner);
};

} // namespace
