#include "HuffTree.h"

namespace mtc{

TreeNode::TreeNode()
{
	leafvalue = -1; left = NULL; right=NULL;
}
TreeNode::~TreeNode()
{
}
HuffTree::HuffTree()
{
	root = new TreeNode();
}
HuffTree::~HuffTree()
{

}
void HuffTree::Insert(const string& binary, int value)
{
	Insert(root,binary,value);
}

void HuffTree::Insert(TreeNode* pos, const string& binary, int value)
{
	if(binary.empty())
	{
		pos->leafvalue = value;
		return;
	}
	else
	{
		if (binary[0] == ' ') //remove space
			Insert(pos,binary.substr(1,binary.length()-1),value);
		else if(binary[0]=='0') // to left
		{
			TreeNode* temp = new TreeNode();
			pos->left = temp;
			pos=temp;
			Insert(pos,binary.substr(1,binary.length()-1),value);
		}
		else //to right
		{
			TreeNode* temp = new TreeNode();
			pos->right = temp;
			pos=temp;
			Insert(pos,binary.substr(1,binary.length()-1),value);
		}
	}
}

int HuffTree::Search(TreeNode* pos, const string& binary)
{
	if(binary.empty())
	{
		return pos->leafvalue;
	}
	else
	{
		if(binary[0]=='0') // to left
		{
			if (pos->left!=NULL)
				return	Search(pos->left,binary.substr(1,binary.length()-1));
			else
				return -2;
		}
		else //to right
		{
			if (pos->right!=NULL)
				return Search(pos->right,binary.substr(1,binary.length()-1));
			else
				return -2;
		}
	}
}

TreeNode* HuffTree::Search(TreeNode* pos, const char bit)
{
	if(bit==0)
		return NULL;
	else
	{
		if(bit=='0') // to left
		{
			return pos->left;
		}
		else
		{
			return pos->right;
		}
	}
}
}
