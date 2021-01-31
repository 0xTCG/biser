/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Authors: alimg, inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>

#include "common.h"
#include "hit.h"

/******************************************************************************/

std::vector<Hit> merge(std::vector<Hit> &hits, const int merge_dist);

void merge_main(int argc, char **argv);
