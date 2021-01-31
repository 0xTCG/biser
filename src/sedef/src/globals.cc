/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag

/******************************************************************************/

#include "globals.h"

/******************************************************************************/

using namespace Globals;

int Search::KMER_SIZE = 12;
int Search::WINDOW_SIZE = 16;
int Search::MIN_UPPERCASE = Search::KMER_SIZE;

double Search::MAX_ERROR = 0.30;
double Search::MAX_EDIT_ERROR = 0.15;
double Search::GAP_FREQUENCY = 0.005;
int Search::MIN_READ_SIZE = KB * (1 - Search::MAX_ERROR); // 700 by default

int Align::MATCH = 5;
int Align::MISMATCH = -4;
int Align::GAP_OPEN = -40;
int Align::GAP_EXTEND = -1;

int Chain::MAX_CHAIN_GAP = Search::MAX_ERROR * Search::MIN_READ_SIZE;

double Extend::RATIO = 5;
int Extend::MAX_EXTEND = 15 * KB;
int Extend::MERGE_DIST = 250;

int Stats::MAX_OK_GAP = -1;
int Stats::MIN_SPLIT_SIZE = KB;
int Stats::MIN_UPPERCASE = 100;
double Stats::MAX_SCALED_ERROR = 0.5;
