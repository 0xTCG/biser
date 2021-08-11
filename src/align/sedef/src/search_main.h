/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <string>

#include "common.h"
#include "hash.h"
#include "search.h"

/******************************************************************************/

std::vector<std::vector<std::string>>
generate_translation(const std::string &ref_path, bool print = false);
void search_main(int argc, char **argv);
void trans_main(int argc, char **argv);
