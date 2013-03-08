// modified_preorder_tree_traversal.cpp
// Author: Ruchira S. Datta
// Copyright (c) 2008, Regents of the University of California
// All rights reserved.
//
// Redistiribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// o Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// o Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// o Neither the name of the University of California, Berkeley nor the names
// of its contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR 
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "modified_preorder_tree_traversal.h"
#include <iterator>
using namespace std;

int findLeftRightIds(const bpp::Tree &tree, 
                        LeftRightIdsOfNodeIdMap &idMap,
                        map<int, int> &nodeIdOfLeftIdMap,
                        map<int, int> &levelOfNodeIdMap,
                        int nodeId, int leftId, int level) {
  levelOfNodeIdMap[nodeId] = level;
  int rightId = leftId + 1;
  const vector<int> &children = tree.getSonsId(nodeId);
  vector<int>::const_iterator node_iter;
  for( node_iter = children.begin(); node_iter != children.end();
      ++node_iter ) {
    rightId = findLeftRightIds(tree, idMap, nodeIdOfLeftIdMap, 
                                levelOfNodeIdMap, *node_iter, 
                                rightId, level + 1);
  }

  idMap[nodeId] = make_pair(leftId, rightId);
  nodeIdOfLeftIdMap[leftId] = nodeId;
  return (rightId + 1);
}
