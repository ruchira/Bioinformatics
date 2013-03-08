// get_taxon_of_leaf_map_from_file.h
// Author: Ruchira S. Datta
// Copyright (c) 2011, Regents of the University of California
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

#ifndef GET_TAXON_OF_LEAF_MAP_FROM_FILE_H
#define GET_TAXON_OF_LEAF_MAP_FROM_FILE_H

#include "orthologs_common.h"
#include <string>
#include <map>
#include <vector>

using namespace std;

extern void getTaxonOfLeafMapFromFile(const string &taxa_file, 
        const map<DBIdentifierT, int> &leaf_of_sequence_header,
        const map<string, DBIdentifierT> &sequence_header_of_leaf_name,
        map<int, DBIdentifierT> &species_of_node,
        vector<string> &taxon_names,
        map<string, DBIdentifierT> &taxon_of_taxon_name);

#endif
