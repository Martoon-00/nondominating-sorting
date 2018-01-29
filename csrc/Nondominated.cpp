#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <memory>

using namespace std;

//
// Range
//
template <class V>
struct range_t {
    V left;
    V right;

    range_t(V l, V r): left(l), right(r) {}

    static range_t<V> thin(V value) {
        return range_t<V>(value, value);
    }

    bool contains(V value) {
        return left <= value && value <= right;
    }

    bool operator==(range_t o){
        return left == o.left && right == o.right;
    }

};

template <class V>
ostream& operator<<(ostream& os, range_t<V> r) {
    return os << "[" << r.left << ", " << r.right << "]";
}

//
// Utilities
//

template <class T>
function<T(T)> identity() {
    return [](T x){ return x; };
}

template <class V>
ostream& operator<<(ostream& os, vector<V> vs) {
    os << "{";
    for (size_t i = 0; i < vs.size(); i++) {
        os << vs[i];
        if (i != vs.size() - 1)
            os << ", ";
    }
    return os << "}";
}

template <class T>
T findMedian(vector<T>& list) {
    if (list.size() < 25) {
        sort(list.begin(), list.end());
        return list[list.size() / 2];
    }

    auto result = vector<T>();
    for (size_t i = 0; i + 5 < list.size(); i += 5) {
        sort(list.begin() + i, list.begin() + i + 5);
        result.push_back(list[i + 2]);
    }
    return findMedian(result);
}

template <class V>
vector<V> removeAdjacentDuplicates(vector<V>& list) {
    vector<V> unique = vector<V>();
    if (list.size() > 0) {
        unique.push_back(list[0]);
        for (size_t i = 1; i < list.size(); i++) {
            if (list[i - 1] != list[i]) {
                unique.push_back(list[i]);
            }
        }
    }
    return unique;
}

template <class V>
struct split_t {
    vector<V> L;
    vector<V> M;
    vector<V> R;

    split_t(vector<V>& l, vector<V>& m, vector<V>& r): L(l), M(m), R(r) {}
    split_t() {};
};

template <class V, class C>
split_t<V> split(vector<V> list, C med, function<C(V)> convert) {
    split_t<V> result = split_t<V>();
    for (V value : list) {
        auto coord = convert(value);
        if (coord < med) {
            result.L.push_back(value);
        } else if (coord > med) {
            result.R.push_back(value);
        } else {
            result.M.push_back(value);
        }
    }
    return result;
}

template <class V>
std::ostream &operator<<(std::ostream &os, split_t<V> const &split) {
    os << "L: ";
    for (auto it = split.L.begin(); it < split.L.end(); it++)
        os << *it << ", ";
    os << endl << "M: ";
    for (auto it = split.M.begin(); it < split.M.end(); it++)
        os << *it << ", ";
    os << endl << "R: ";
    for (auto it = split.R.begin(); it < split.R.end(); it++)
        os << *it << ", ";
    os << endl;

    return os;
}

template <class V>
struct cached {
    bool evaluated = false;
    V cache;
    function<V()> evaluate;

    cached(const function<V()>& eval): evaluate(eval) {};

    V operator()() {
        if (evaluated) {
            return cache;
        } else {
            auto res = evaluate();
            cache = res;
            evaluated = true;
            return res;
        }
    }
};


struct timer {
    clock_t start;

    timer(): start(clock()) {}

    double sec() {
        return (clock() - start) / (double) CLOCKS_PER_SEC;
    }

};
ostream& operator<<(ostream &os, timer t) {
    os << t.sec() << "s";
    return os;
}

//
// point
//

typedef int coord_t;

struct point {
    static int COORD_X;
    static int COORD_Y;

    vector<coord_t> coords;

    point(vector<coord_t>& cs): coords(cs) {}
    point(vector<coord_t>&& cs): coords(move(cs)) {}

    coord_t operator[](int i) {
        return coords[i];
    }

    size_t dimensions() {
        return coords.size();
    }

    bool dominates(point& o) {
        bool meet_strict = false;
        for (size_t i = 0; i < coords.size(); i++) {
            if (coords[i] > o.coords[i])
                return false;
            if (coords[i] < o.coords[i])
                meet_strict = true;
        }
        return meet_strict;
    }

    bool operator==(point& o) {
        return coords == o.coords;
    }
};

int point::COORD_X = 0;
int point::COORD_Y = 1;
ostream& operator<<(ostream &os, point p) {
    os << "(";
    for (size_t i = 0; i < p.dimensions(); i++) {
        os << p[i];
        if (i + 1 != p.dimensions())
            os << ", ";
    }
    os << ")";
    return os;
}

struct point_id {
    int i;
    explicit point_id(int i_): i(i_) {};
    point_id(const point_id& p): i(p.i) {};
    point_id(): i(-1) {};
    point_id& operator=(const point_id& o) {
        i = o.i;
        return *this;
    }
    explicit operator int() {
        return i;
    }
    bool operator<(const point_id& o) const {
        return i < o.i;
    }
    bool operator==(const point_id& o) const {
        return i == o.i;
    }
    bool operator!=(const point_id& o) const {
        return i != o.i;
    }
};
ostream& operator<<(ostream &os, point_id id) {
    return os << "#" << id.i;
}

// In algorithm we numerate ranks stating from 1, 0 indicates unknown value.
// It's convenient in initialization.
typedef size_t rang;

//
// Segment tree
//
template <class K, class V, class M>
struct seg_tree {
    virtual K get_key(V) = 0;
    virtual M get_sum(V) = 0;
    virtual M add_sum(M, M) = 0;

    struct tree {
        seg_tree* outer;
        range_t<K> range;
        M summary;

        tree(seg_tree* o, range_t<K>& r): outer(o), range(r), summary() {
        }

        virtual M request(::range_t<K> reqRange) = 0;

        virtual void modify(K key, function<void(set<V>&)> modifier) = 0;

        virtual void forEach(function<void(V)> consumer) = 0;
    };

    struct leaf: tree {
        set<V> values;

        leaf(seg_tree* o, range_t<K> range): tree(o, range) {};

        M request(range_t<K> reqRange) override {
            if (reqRange.contains(this -> range.left)) {
                return this -> summary;
            } else {
                return M();
            }
        }

        void modify(K key, function<void(set<V>&)> modifier) override {
            if (!this -> range.contains(key))
                return;
            modifier(this -> values);
            updateSummary();
        }

        void updateSummary() {
            M summary = M();
            for (V value : values) {
                summary = this -> outer -> add_sum(summary, this -> outer -> get_sum(value));
            }
            this -> summary = summary;
        }

        void forEach(function<void(V)> consumer) override {
            for_each(values.begin(), values.end(), consumer);
        }
    };

    struct node: tree {
        tree *left;
        tree *right;
        function<bool(K)> isLeft;

        node(seg_tree* o, tree* l, tree* r, function<bool(K)> il, range_t<K> rng)
            : tree(o, rng), left(l), right(r), isLeft(il)
        {}

        M request(range_t<K> reqRange) override {
            if (reqRange == this -> range) {  // totally matches
                return this -> summary;
            } else if (reqRange.left > this->range.right
                       || this->range.left > reqRange.right) {  // totally in nowhere
                return M();
            } else if (reqRange.right < this->right->range.left) {  // go left
                return left -> request(reqRange);
            } else if (reqRange.left > this->left->range.right) {  // go right
                return right -> request(reqRange);
            } else {  // go in both directions
                M leftRes = left -> request(range_t<K>(reqRange.left, this -> range.right));
                M rightRes = right -> request(range_t<K>(this -> range.left, reqRange.right));
                return this -> outer -> add_sum(leftRes, rightRes);
            }
        }

        void modify(K key, function<void(set<V>&)> modifier) override {
            if (this -> left -> range.contains(key)) {
                this -> left -> modify(key, modifier);
            } else {
                this -> right -> modify(key, modifier);
            }
            this -> summary = this -> outer -> add_sum(this -> left -> summary, this -> right -> summary);
        }

        void forEach(function<void(V)> consumer) override {
            this -> left -> forEach(consumer);
            this -> right -> forEach(consumer);
        }
    };

    tree* str;

    void form(vector<V>& values) {
        if (!values.size())
            throw runtime_error("Can't build empty segment tree");

        vector<K> keys = vector<K>();
        for (V value : values) {
            keys.push_back(get_key(value));
        }
        sort(keys.begin(), keys.end());
        vector<K> uniqKeys = removeAdjacentDuplicates(keys);
        this -> buildTree(uniqKeys);
    }

    seg_tree<K, V, M>& buildTree(vector<K>& sortedKeys) {
        this -> str = buildTreeImpl(sortedKeys.begin(), sortedKeys.end());
        return *this;
    }

    tree* buildTreeImpl(typename vector<K>::iterator const begin, typename vector<K>::iterator const end) {
        if (begin == end)
            throw runtime_error("Can't build empty seg tree");

        if (begin + 1 == end) {
            return new leaf(this, range_t<K>::thin(*begin));
        } else {
            auto range = range_t<K>(*begin, *(end - 1));
            auto med = begin + (end - begin) / 2;
            auto med_value = *med;
            return new node(
                this,
                buildTreeImpl(begin, med),
                buildTreeImpl(med, end),
                [&med_value](K k){ return k < med_value; },
                range
                );
        }
    }

    M request(range_t<K> range) {
        return str -> request(range);
    }

    void insert(V value) {
        str -> modify(get_key(value), [&value](set<V>& s){ s.insert(value); });
    }

    void remove(V value) {
        str -> modify(get_key(value), [&value](set<V>& s){ s.erase(value); });
    }

};

template <class K, class V, class M>
ostream& operator<<(ostream &os, seg_tree<K, V, M>& tree) {
    auto values = vector<V>();
    tree.str -> forEach([&os](V v){ os << v << " "; });
    return os << endl;
}


//
// Solution
//

struct solver {
    // id -> point
    vector<point>& points;
    // id -> rank
    vector<rang>& ranks;

    solver(vector<point>& ps): points(ps), ranks(*(new vector<rang>(points.size()))) {
    }

    solver(vector<point>& ps, vector<rang>& rs): points(ps), ranks(rs) {
        if (!points.size())
            throw runtime_error("Can't solve for no points");
    }

    rang get_rank(point_id i) {
        return ranks[int(i)];
    }

    void update_rank(point_id i, rang new_rank) {
        if (ranks[int(i)] < new_rank)
            ranks[int(i)] = new_rank;
    }

    point& getPoint(point_id i) {
        return points[int(i)];
    }

    virtual void solve_for_ids(vector<point_id>& ids) = 0;

    vector<rang> solve() {
        vector<point_id> ids = vector<point_id>();
        for (size_t i = 0; i < points.size(); i++) {
            ids.push_back(point_id(i));
        }
        solve_for_ids(ids);

        for (auto id: ids) {
            if (get_rank(id) == 0)
                throw runtime_error("Rank remained unknown");
            ranks[int(id)]--;
        }
        return ranks;
    }

};

struct dumb_solver: solver {
    dumb_solver(vector<point>& points): solver(points) {
        ranks.resize(points.size());
    }

    dumb_solver(vector<point>& points, vector<rang>& ranks): solver(points, ranks) { }

    void solve_for_ids(vector<point_id> &ids) override {
        sort_dumb(ids);
    }

    // reusable vector with rank evaluators.
    static vector<unique_ptr<cached<rang>>> rank_evals;

    void sort_dumb(vector<point_id>& ids) {
        if (ids.size() == 0)
            return;
        if (ids.size() == 1) {
            update_rank(ids[0], 1);
            return;
        }

        if (rank_evals.size() != points.size()) 
            rank_evals.resize(points.size());

        for (auto id : ids) {
            rank_evals[int(id)] = unique_ptr<cached<rang>>(
              new cached<rang>([id, this, &ids](){
                rang max_rank = 0;
                for (auto dep : ids) {
                    if (getPoint(dep).dominates(getPoint(id))) {
                        max_rank = max(max_rank, (*rank_evals[int(dep)])());
                    }
                }
                return max_rank + 1;
            }));
        }
        for (auto id : ids) {
            update_rank(id, (*rank_evals[int(id)])());
        }

    }

    void update_dumb(vector<point_id>& known, vector<point_id>& request) {
        if (!known.size())
            return;

        for (auto req : request) {
            for (auto kn : known) {
                if (getPoint(kn).dominates(getPoint(req))){
                    update_rank(req, get_rank(kn) + 1);
                }
            }
        }
    }
};
vector<unique_ptr<cached<rang>>> dumb_solver::rank_evals = vector<unique_ptr<cached<rang>>>();

struct sac_solver: solver {

    sac_solver(vector<point>& points): solver(points) {}

    void solve_for_ids(vector<point_id>& ids) override {
        sort(ids.begin(), ids.end(), 
             [&](const point_id& id1, const point_id& id2){ return lex_sort_cmp(id1, id2); });
        sortSAC(points[0].dimensions(), ids);
    }

    coord_t getX(point_id id) {
        return getPoint(id)[point::COORD_X];
    }
    coord_t getY(point_id id) {
        return getPoint(id)[point::COORD_Y];
    }

    bool lex_sort_cmp(const point_id& id1, const point_id& id2) {
        point p1 = getPoint(id1);
        point p2 = getPoint(id2);
        coord_t x1 = p1[point::COORD_X];
        coord_t y1 = p1[point::COORD_Y];
        coord_t x2 = p2[point::COORD_X];
        coord_t y2 = p2[point::COORD_Y];
        return x1 < x2 || (x1 == x2 && y1 < y2);
    };


    void sortSAC(int dimension, vector<point_id>& ids) {
        // cerr << "sortSAC " << dimension << " " << ids << endl;

        if (dimension <= 1) {
            throw runtime_error("Can't sort for so few dimensions");
        } else if (ids.size() == 0){
        } else if (ids.size() == 1) {
            update_rank(ids[0], 1);
        } else if (ids.size() < 500) {
            dumb_solver(points, ranks).solve_for_ids(ids);
        } else if (dimension == 2) {
            sortSweep(ids);
        } else {
            // Helpers
            function<coord_t(point_id)> getPointCoord = [&](point_id id) {
                return getPoint(id)[dimension - 1];
            };

            // Get median
            auto zs = vector<coord_t>();
            for (point_id id : ids) {
                zs.push_back(getPointCoord(id));
            }
            coord_t median = findMedian(zs);

            // Split
            split_t<point_id> cut = split(ids, median, getPointCoord);

            // Call recursively
            sortSAC(dimension, cut.L);
            updateSAC(dimension, cut.L, cut.M);
            sortSAC(dimension - 1, cut.M);
            updateSAC(dimension, cut.L, cut.R);
            updateSAC(dimension, cut.M, cut.R);
            sortSAC(dimension, cut.R);
        }

    }

    void updateSAC(int dimension, vector<point_id>& known, vector<point_id>& request) {
        // cerr << "updateSAC " << dimension << " " << known << " " << request << endl;

        if (dimension <= 1) {
            throw runtime_error("Can't update for so few dimensions");
        } else if (request.size() == 0 || known.size() == 0) {
        } else if (known.size() + request.size() <= 500) {
            dumb_solver(points, ranks).update_dumb(known, request);
        } else if (dimension == 2) {
            updateSweep(known, request);
        } else {
            // Helpers
            function<coord_t(point_id)> getPointCoord = [&](point_id id) {
                return getPoint(id)[dimension - 1];
            };

            // Get median
            auto zs = vector<coord_t>();
            for (point_id id : request) {
                zs.push_back(getPointCoord(id));
            }
            coord_t median = findMedian(zs);

            // Split
            split_t<point_id> knownCut = split(known, median, getPointCoord);
            split_t<point_id> requestCut = split(request, median, getPointCoord);

            // Recursive calls
            updateSAC(dimension, knownCut.L, requestCut.L);
            updateSAC(dimension - 1, knownCut.L, requestCut.M);
            updateSAC(dimension - 1, knownCut.L, requestCut.R);
            updateSAC(dimension - 1, knownCut.M, requestCut.R);
            updateSAC(dimension, knownCut.R, requestCut.R);
        }

    }

    typedef map<coord_t, point_id> line_t;

    struct ranks_tree_t: seg_tree<coord_t, point_id, rang> {
        sac_solver* sac;

        ranks_tree_t(sac_solver* sac_): sac(sac_) {};

        coord_t get_key(point_id id) { return sac -> getY(id); }
        rang get_sum(point_id id) { return sac -> get_rank(id); }
        rang add_sum(rang a, rang b) { return max(a, b); }
    };

    ranks_tree_t* makeRanksTree(vector<point_id>& ids) {
        auto* tree = new ranks_tree_t(this);
        tree -> form(ids);
        return tree;
    }

    int get_next_rank_after(ranks_tree_t* ranks_tree, point_id id) {
        rang maxRank = ranks_tree -> request(range_t<coord_t>(numeric_limits<int>::min(), getY(id)));
        return maxRank + 1;
    }

    void insert_point(line_t* line, ranks_tree_t* ranks_tree, vector<point_id>* front, point_id id) {
        if (line != nullptr) {
            line -> insert(pair<coord_t, point_id>(getY(id), id));
        }
        if (ranks_tree != nullptr) {
            ranks_tree -> insert(id);
        }
        if (front != nullptr) {
            (*front)[get_rank(id)] = id;
        }
    }

    void remove_point(line_t* line, ranks_tree_t* ranks_tree, vector<point_id>* front, point_id id) {
        if (line != nullptr) {
            line -> erase(getY(id));
        }
        if (ranks_tree != nullptr) {
            ranks_tree -> remove(id);
        }
        if (front != nullptr) {
            (*front)[get_rank(id)] = point_id();
        }
    }

   
    void sortSweep(vector<point_id>& ids) {
        auto line = line_t();
        auto ranks_tree = makeRanksTree(ids);

        point_id prevId = point_id();
        for (point_id id : ids) {
            if (prevId != point_id() && getPoint(prevId) == getPoint(id)) {
                update_rank(id, get_rank(prevId));
            } else {
                coord_t key = getY(id);
                auto higher_point = line.upper_bound(key);

                rang new_rank = get_next_rank_after(ranks_tree, id);

                if (higher_point != line.end() && new_rank == get_rank(higher_point -> second)) {
                    remove_point(&line, ranks_tree, nullptr, higher_point -> second);
                }

                update_rank(id, new_rank);
                insert_point(&line, ranks_tree, nullptr, id);
            }
            prevId = id;
        }
    }

    // reusable front, clear after use.
    vector<point_id> front = vector<point_id>(points.size());

    void updateSweep(vector<point_id>& knowns, vector<point_id>& requests) {
        // cerr << "updateSweep " << knowns << " " << requests << endl;

        auto line = map<coord_t, point_id>();
        auto ranks_tree = makeRanksTree(knowns);

        size_t known_i = 0;
        for (point_id req : requests) {
            while (known_i < knowns.size() && !lex_sort_cmp(req, knowns[known_i])) {
                auto known = knowns[known_i];
                known_i++;

                // processing known point
                rang rank = get_rank(known);
                point_id old_known = front[rank];
                if (old_known == point_id()) {
                    insert_point(&line, ranks_tree, &front, known);
                } else if (getY(old_known) > getY(known)) {
                    remove_point(&line, ranks_tree, &front, old_known);
                    insert_point(&line, ranks_tree, &front, known);
                } else {
                }
            }

            // processing request point
            int new_rank = get_next_rank_after(ranks_tree, req);
            update_rank(req, new_rank);
        }

        for (auto kn: knowns) {
            front[get_rank(kn)] = point_id();
        }
    }

};

//
// Random
//
struct gen {
    static size_t seed;

    static void refresh_seed(size_t s) {
        seed = s * 237 * 237;
    }

    static int int_g(size_t k) {
        size_t a = (seed + 99139) * 99139;
        seed = a / 95279;
        size_t n = a % 95279;
        return (int) (n % (2 * k + 1)) - k;
    }

    static point point_g(int d) {
        auto v = vector<coord_t>();
        for (int i = 0; i < d; i++) {
            v.push_back(int_g(100));
        }
        return point(v);
    }

    static vector<point> points_g(int n, int d) {
        auto points = vector<point>();
        for (int i = 0; i < n; i++)
            points.push_back(point_g(d));
        return points;
    }
};
size_t gen::seed;

bool test_correctness() {
    for (int i = 1; i <= 10000; i++) {
        cout << "Iteration #" << i << endl;
        gen::refresh_seed(i);
        auto points = gen::points_g(100, 4);
        auto ranks_sac = sac_solver(points).solve();
        auto ranks_dumb = dumb_solver(points).solve();
        if (ranks_sac != ranks_dumb) {
            cout << "For " << points << endl;
            cout << "expected " << ranks_dumb << ", got" << ranks_sac << endl;
            return false;
        }
        cout << "OK!" << endl;
    }
    return true;
}


int main(){
    timer timer;

    gen::refresh_seed(3234);
    auto points = gen::points_g(100000, 4);
    sac_solver(points).solve();

    cout << timer << endl;


    // test_correctness();

    return 0;
}

