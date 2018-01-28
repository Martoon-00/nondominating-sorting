#include <iostream>
#include <stdlib.h>

using namespace std;

template <class V>
struct range {
    V left;
    V right;

    range(V l, V r): left(l), right(r) {}

    template <V>
    static range<V> thin(V value) {
        return range(value, value);
    }

    bool contains(V value) {
        return left <= value && value <= right;
    }

    bool operator==(range o){
        return left == o.left && right == o.right;
    }

};

int main(){
    auto r = range<int>(1, 3);
    int k = 1;
    cout << r.contains(0) << " " << r.contains(k) << " " << r.contains(2) << endl;

    return 0;
}

