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

//
// Utilities
//

template <class T>
function<T(T)> identity() {
    return [](T x){ return x; };
}

template <class T>
T findMedian(vector<T>& list) {
    if (list.size() < 25) {
        sort(list.begin(), list.end());
        return list.get(list.size() / 2);
    }

    auto result = vector<T>();
    for (int i = 0; i + 5 < list.size(); i += 5) {
        list.sort(list.begin() + i, list.begin() + i + 5);
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

template <class V>
split_t<V> split(vector<V> list, V med, bool less(V, V)) {
    split_t<V> result = split_t<V>();
    for (V value : list) {
        if (less(value, med)) {
            result.L.push_back(value);
        } else if (less(med, value)) {
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
// Segment tree
//
template <class K, class V, class M>
struct seg_tree {
    function<K(V)> get_key;
    function<M(V)> get_sum;
    function<M(M,M)> add_sum;

    struct tree {
        seg_tree* outer;
        range_t<K> range;
        M summary;

        tree(seg_tree* o, range_t<K> r): outer(o), range(r), summary() {
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
                throw logic_error("Modification visited leaf in nowhere");
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

    static seg_tree* make(function<K(V)> get_key, function<M(V)> get_sum, function<M(M,M)> add_sum, vector<V>& values) {
        if (!values.size())
            throw runtime_error("Can't build empty segment tree");

        vector<K> keys = vector<K>();
        for (V value : values) {
            keys.push_back(get_key(value));
        }
        sort(keys.begin(), keys.end());
        vector<K> uniqKeys = removeAdjacentDuplicates(keys);
        auto tree = new seg_tree<K, V, M>(get_key, get_sum, add_sum);
        tree -> buildTree(uniqKeys);
        return tree;
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

    seg_tree(function<K(V)> get_key_, function<M(V)> get_sum_, function<M(M, M)> add_sum_)
        : get_key(get_key_), get_sum(get_sum_), add_sum(add_sum_)
        {}

    M request(range_t<K> range) {
        return str -> request(range);
    }

    void insert(V value) {
        str -> modify(get_key(value), [&value](set<V>& s){ s.insert(value); });
    }

    void remove(V value) {
        str -> modify(get_key(value), [&value](set<V>& s){ s.remove(value); });
    }

};

template <class K, class V, class M>
ostream& operator<<(ostream &os, seg_tree<K, V, M>& tree) {
    auto values = vector<V>();
    tree.forEach([&os](V v){ os << v << " "; });
    return os << endl;
}






int main(){
    vector<int> values = vector<int>{1,2,3,4,5};
    seg_tree<int, int, int>* tree = seg_tree<int, int, int>::make(identity<int>(), identity<int>(), [](int a, int b){ return a + b; }, values);
    tree -> insert(3);
    cout << tree -> request(range_t<int>(2, 3)) << endl;
    return 0;
}

