/*
 * kdtree.hpp
 *
 *  Created on: Nov 27, 2016
 *      Author: shayan
 */

#ifndef INCLUDE_KDTREE_HPP_
#define INCLUDE_KDTREE_HPP_

#include "object.hpp"

class TreeNode
{
public:
	BoundingBox bbox;
	TreeNode *left;
	TreeNode *right;
	std::vector<const Triangle*> members;

	TreeNode(){}
	static TreeNode* build(std::vector<const Triangle*> &members_in,
			const BoundingBox &bbox_in,
			const int depth,
			const int max_depth);

	// Writing to a file
	void write_vtk(Vector * = NULL);
	void get_members_bbox(std::vector<const BoundingBox*>& bboxvec);

	// Free memory
	void free_children();

};



#endif /* INCLUDE_KDTREE_HPP_ */
