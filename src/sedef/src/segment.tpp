/// 786

/// This file is subject to the terms and conditions defined in
/// file 'LICENSE', which is part of this source code package.

/// Author: inumanag
/// Range segment tree: 
/// https://pdfs.semanticscholar.org/6a87/0c8b438174c2f08fc69f3883b73e507b2dea.pdf

/******************************************************************************/

template<typename T>
SegmentTree<T>::SegmentTree(std::vector<T> &a): 
	anchors(a), activated(0)
{
	std::sort(anchors.begin(), anchors.end());

	int size = (1 << (32 - __builtin_clz(anchors.size() - 1)));
	tree.resize(size << 1);

	int tree_i = 0;
	int m = initialize(0, 0, anchors.size(), tree_i);
	m++;
	assert(tree_i == anchors.size());
	assert(m <= tree.size());
}

template<typename T>
int SegmentTree<T>::rmq(const SegmentTree<T>::Tp &p, const SegmentTree<T>::Tp &q, int i) const // [p, q] are both inclusive
{
	if (i >= tree.size()) {
		return -1;
	} else if (tree[i].a != -1) { // leaf
		if (p <= anchors[tree[i].a].x && anchors[tree[i].a].x <= q) {
			return i;
		} else {
			return -1;
		}
	} else {
		int pv = tree[i].p;
		if (pv == -1) {
			return -1; // nothing in [0, q]
		}
		assert(tree[pv].a != -1);
		if (p <= anchors[tree[pv].a].x && anchors[tree[pv].a].x <= q) {
			return pv;
		} else {
			assert(2 * i + 1 < tree.size());
			if (q <= tree[2 * i + 1].h) { // h is inclusive
				return rmq(p, q, 2 * i + 1);
			} else if (p > tree[2 * i + 1].h) {
				return rmq(p, q, 2 * i + 2);
			} else {
				int m1 = rmq(p, q, 2 * i + 1);
				int m2 = rmq(p, q, 2 * i + 2); 
				if (m1 == -1) 
					return m2;
				if (m2 == -1) 
					return m1;
				assert(tree[m1].a != -1);
				assert(tree[m2].a != -1);
				return (anchors[tree[m1].a].score >= anchors[tree[m2].a].score) ? m1 : m2;
			}
		}
	}
}

template<typename T>
int SegmentTree<T>::rmq(const SegmentTree<T>::Tp &p, const SegmentTree<T>::Tp &q) const
{
	int i = rmq(p, q, 0);
	return (i == -1 ? -1 : tree[i].a);
}

template<typename T>
void SegmentTree<T>::activate(const SegmentTree<T>::Tp &q, int score)
{
	int leaf = 0;
	for (leaf = 0; leaf < tree.size() && (tree[leaf].a == -1 || q != anchors[tree[leaf].a].x); ) {
		leaf = 2 * leaf + 1 + (q > tree[2 * leaf + 1].h);
	}
	assert(leaf < tree.size());
	assert(q == tree[leaf].h);
	assert(tree[leaf].a != -1); // leaf
	anchors[tree[leaf].a].score = score;

	for (int i = 0; i < tree.size(); ) {
		assert(tree[leaf].a != -1);
		if (tree[i].p == -1 || anchors[tree[leaf].a].score >= anchors[tree[tree[i].p].a].score)
			std::swap(tree[i].p, leaf);
		assert(tree[i].p != -1);
		if (leaf == -1) 
			break;

		assert(tree[leaf].a != -1);
		assert(2 * i + 1 < tree.size());
		i = 2 * i + 1 + (anchors[tree[leaf].a].x > tree[2 * i + 1].h);
	}

	activated++;
	assert(activated <= anchors.size());
}

template<typename T>
void SegmentTree<T>::deactivate(const SegmentTree<T>::Tp &q)
{
	int leaf = 0;
	for (leaf = 0; leaf < tree.size() && (tree[leaf].a == -1 || q != anchors[tree[leaf].a].x); ) {
		leaf = 2 * leaf + 1 + (q > tree[2 * leaf + 1].h);
	}
	assert(leaf < tree.size());
	assert(q == tree[leaf].h);
	assert(tree[leaf].a != -1); // leaf
	anchors[tree[leaf].a].score = SegmentTree<T>::MIN;;

	for (int i = 0; i < tree.size(); ) {
		if (tree[i].p == -1) {
			break;
		} else if (tree[i].p == leaf) {
			if (tree[i].a != -1) { // leaf
				tree[i].p = -1;
			} else {
				// ^_^
				assert(2 * i + 1 < tree.size());
				if (2 * i + 2 < tree.size() && 
					tree[2 * i + 2].p != -1 &&
					(tree[2 * i + 1].p == -1 ||
					 anchors[tree[tree[2 * i + 2].p].a].score > anchors[tree[tree[2 * i + 1].p].a].score))
				{
					tree[i].p = leaf = tree[2 * i + 2].p;
					i = 2 * i + 2;
				} else {
					tree[i].p = leaf = tree[2 * i + 1].p;
					i = 2 * i + 1;
				}
			}
		} else {
			i = 2 * i + 1 + (q > tree[2 * i + 1].h);
		}
	}

	activated--;
	assert(activated >= 0);
}

template<typename T>
int SegmentTree<T>::initialize(int i, int s, int e, int &tree_i)
{
	// assert(i < tree.size());
	if (i >= tree.size()) {
		return -1;
	} else if (s + 1 == e) {
		assert(tree_i < anchors.size());
		tree[i] = Point(-1, tree_i, anchors[tree_i].x);
		anchors[tree_i].score = SegmentTree<T>::MIN;
		tree_i++;
		return i;
	} else {
		int bnd = (s + e + 1) / 2;
		int a = initialize(2 * i + 1, s, bnd, tree_i);
		int b = initialize(2 * i + 2, bnd, e, tree_i);
		// assert(2 * i + 1 < tree.size());
		tree[i] = Point(-1, -1, tree[2 * i + 1 + (2 * i + 2 < tree.size())].h);
		return std::max(a, std::max(i, b));
	}
}

template<typename T>
bool SegmentTree<T>::empty() const
{
	return (activated == 0);
}

/******************************************************************************/

template<typename T>
void SegmentTree<T>::plot(int w, int l, int i, int s, int e, std::vector<std::vector<std::string>> &PLOT)
{
	if (i >= tree.size()) return;
	int bnd = (s + e + 1) / 2;
	if (tree[i].a == -1) plot(w/2, l+1, 2*i+1, s, bnd, PLOT);
	PLOT[0][l] += fmt::format(
		fmt::format("{{:^{}}}", w),
		fmt::format("{}/{}{}", tree[i].h.first, tree[i].h.second,
			tree[i].a == -1 ? "" : "*")
	);
	PLOT[1][l] += fmt::format(
		fmt::format("{{:^{}}}", w),
		fmt::format("{}", 
			tree[i].p != -1 
				? fmt::format("{}/{}", 
					anchors[tree[tree[i].p].a].x.first, anchors[tree[tree[i].p].a].x.second) 
				: tree[i].a != -1 ? fmt::format("({})", 
					anchors[tree[i].a].score == SegmentTree<T>::MIN ? -1 : anchors[tree[i].a].score) : "")
	);
	if (tree[i].a == -1) plot(w/2, l+1, 2*i+2, bnd, e, PLOT);
}

template<typename T>
std::string SegmentTree<T>::plot() const
{
	std::vector<std::vector<std::string>> PLOT(2, std::vector<std::string>(50));

	int w = 6 * pow(2, ceil(log(tree.size()) / log(2))-1);
	plot(w, 0, 0, 0, anchors.size(),PLOT);

	int cw=w/4, ll=1;
	std::string o = "";
	for (int si=0;si<PLOT[0].size()&&PLOT[0][si]!="";si++) {
		o += fmt::format("{}\n", PLOT[0][si]);
		o += fmt::format("{}\n", PLOT[1][si]);
		if (PLOT[0][si+1]!="") { 
			for (int i = 0; i < ll; i++) {
				std::string h; 
				for (int c=1;c<cw;c++) h+="─";
				o += fmt::format("{0}┌{1}┴{1}┐{0} ",
					std::string(cw-1, ' '), h);
			}
			o += "\n";
			ll*=2;
			cw/=2;
		}
	}
	return o;
}
