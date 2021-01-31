/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag
/// Based on
/// http://www.biorxiv.org/content/biorxiv/early/2017/03/24/103812.full.pdf

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval_set.hpp>

#include "align.h"
#include "hash.h"
#include "hit.h"

/******************************************************************************/

typedef boost::icl::discrete_interval<int> Interval;
typedef boost::icl::interval_map<int, std::set<std::pair<Interval, Interval>>>
    Subtree;
typedef boost::icl::interval_map<int, Subtree> Tree;

/******************************************************************************/

std::vector<Hit> search(int query_winnow_start,
                        std::shared_ptr<Index> query_hash,
                        std::shared_ptr<Index> ref_hash, Tree &tree,
                        const bool same_genome, const int init_len,
                        const bool allow_extend, const bool report_fails);
