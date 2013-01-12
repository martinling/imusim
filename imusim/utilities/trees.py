"""
Support for basic tree structures.
"""
# Copyright (C) 2009-2011 University of Edinburgh
#
# This file is part of IMUSim.
#
# IMUSim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# IMUSim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with IMUSim.  If not, see <http://www.gnu.org/licenses/>.

class TreeNode(object):
    """
    A node in a tree structure.

    @ivar parent: The parent node in the tree structure.
    @ivar children: List of child nodes.
    """
    def __init__(self, parent):
        """
        Construct a new L{TreeNode}

        If parent is not `None`, then the new node will be automatically
        added to the tree structure.

        @param parent: The parent of this node, or `None` for the
            root of a tree.
        """
        self.parent = parent
        self.children = []

    @property
    def parent(self):
        """
        The parent node of this node.

        Setting this property will automatically insert the node into
        the tree structure by adding the node to its parent's list of
        children. If the node had a previous parent it is first removed
        from the old parent's children.
        """
        return self._parent
    @parent.setter
    def parent(self, newParent):
        if getattr(self, '_parent', None) is not None:
            self._parent.children.remove(self)
        if newParent is not None:
            newParent.children.append(self)
        self._parent = newParent

    @property
    def hasParent(self):
        """
        Whether this node has a parent.
        """
        return self.parent is not None

    @property
    def hasChildren(self):
        """
        Whether this node has any children.
        """
        return len(self.children) > 0

    @property
    def root(self):
        """
        The root node of the tree structure this node is part of.
        """
        return self.parent.root if self.hasParent else self

    @property
    def depth(self):
        return 0 if not self.hasParent else (self.parent.depth + 1)

    def ascendTree(self):
        """
        Generator returning the nodes of the tree encountered while ascending
        from the current node to the root node.

        E.g. for the tree::

                A
               / \
              B   C
             / \   \
            D   E   F

        D.ascendTree() would return a generator yielding D, B, A.

        @return: A generator that returns the nodes of the tree encountered
            while ascending from the current node to the root.
        """
        yield self
        if self.hasParent:
            for point in self.parent.ascendTree():
                yield point

    def preorderTraversal(self, preFunc=None, postFunc=None,
            condition=lambda p: True):
        """
        Perform a pre-order depth-first traversal of the tree structure.

        The traversal will start from this node. At each node traversed,
        callback functions will be called with the current node as the
        only argument.

        @param preFunc: Function to call before traversing to the next node.
        @param postFunc: Function to call after traversing the child nodes.
        @param condition: Condition that nodes must meet to be included.

        @return: A generator that returns the nodes of the (sub-)tree rooted
            at this joint, that meet the given condition, in depth first
            pre-order.
        """

        cond = condition(self)
        if cond:
            if preFunc is not None: preFunc(self)
            yield self
        for child in self.children:
            for grandchild in child.preorderTraversal(preFunc,
                    postFunc, condition):
                yield grandchild
        if cond:
            if postFunc is not None: postFunc(self)

    def __iter__(self):
        return self.preorderTraversal()
