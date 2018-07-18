#include "kdtree.hpp"

#include <set>

TreeNode* TreeNode::build(std::vector<const Triangle*> &members_in,
		const BoundingBox &bbox_in,
		const int depth,
		const int max_depth)
{
	TreeNode *treenode = new TreeNode;
	treenode->left = treenode->right = NULL;
	treenode->bbox = bbox_in;

	/*
	* take care of terminating cases
	*/
	if(members_in.size() <= 1)
	{
		treenode->members = members_in;
		//printf("ANN\n");
		return treenode;
	}
	else if(depth >= max_depth)
	{
		treenode->members = members_in;
		//printf("KHAR\n");
		return treenode;
	}

	/*
	 * Break the bounding box into two pieces.
	 */
	BoundingBox bb1, bb2;
	treenode->bbox.divide_longest_axis(bb1, bb2);

	/*
	 * Put each triangle in the corresponding node.
	 */
	std::vector<const Triangle*> mem1, mem2;
	std::set<const Triangle*> set1;
	for (auto it = members_in.begin() ; it != members_in.end() ; ++it)
	{
		const Triangle *tri = *it;
		if(bb1.is_containing(tri->bbox))
		{
			mem1.push_back(tri);
			set1.insert(tri);
		}
		if(bb2.is_containing(tri->bbox)) mem2.push_back(tri);
	}

	/*
	 * Check to stop if more than 50% of triangles are the same.
	 */
	int n_intersec = 0;
	for (auto it = mem2.begin() ; it != mem2.end() ; ++it)
	{
		if(set1.count(*it)) n_intersec++;
	}

	//printf("built node! mem l c: %d mem r c: %d. \n", (int)mem1.size(), (int)mem2.size());

	if( n_intersec < 0.5 * (int)members_in.size())
	{
		treenode->left = build(mem1, bb1,depth+1, max_depth);
		treenode->right = build(mem2, bb2,depth+1, max_depth);
	}
	else
	{
		treenode->members = members_in;
	}

	return treenode;
}

void TreeNode::free_children()
{
	if(left) left->free_children();
	if(right) right->free_children();
	delete left;  left  = NULL;
	delete right; right = NULL;
}

void TreeNode::get_members_bbox(std::vector<const BoundingBox*>& bboxvec)
{
	if(left) left->get_members_bbox(bboxvec);
	if(right) right->get_members_bbox(bboxvec);
	bboxvec.push_back(&this->bbox);
}
