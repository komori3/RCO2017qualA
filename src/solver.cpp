#include <bits/stdc++.h>
#include <random>
#ifdef _MSC_VER
#include <ppl.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro_io **/

/* tuple */
// out
namespace aux {
    template<typename T, unsigned N, unsigned L>
    struct tp {
        static void output(std::ostream& os, const T& v) {
            os << std::get<N>(v) << ", ";
            tp<T, N + 1, L>::output(os, v);
        }
    };
    template<typename T, unsigned N>
    struct tp<T, N, N> {
        static void output(std::ostream& os, const T& v) { os << std::get<N>(v); }
    };
}
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << '[';
    aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t);
    return os << ']';
}

template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x);

/* pair */
// out
template<class S, class T>
std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) {
    return os << "[" << p.first << ", " << p.second << "]";
}
// in
template<class S, class T>
std::istream& operator>>(std::istream& is, const std::pair<S, T>& p) {
    return is >> p.first >> p.second;
}

/* container */
// out
template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) {
    bool f = true;
    os << "[";
    for (auto& y : x) {
        os << (f ? "" : ", ") << y;
        f = false;
    }
    return os << "]";
}
// in
template <
    class T,
    class = decltype(std::begin(std::declval<T&>())),
    class = typename std::enable_if<!std::is_same<T, std::string>::value>::type
>
std::istream& operator>>(std::istream& is, T& a) {
    for (auto& x : a) is >> x;
    return is;
}

/* struct */
template<typename T>
auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) {
    out << t.stringify();
    return out;
}

/* setup */
struct IOSetup {
    IOSetup(bool f) {
        if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); }
        std::cout << std::fixed << std::setprecision(15);
    }
} iosetup(true);

/** string formatter **/
template<typename... Ts>
std::string format(const std::string& f, Ts... t) {
    size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
    std::vector<char> b(l + 1);
    std::snprintf(&b[0], l + 1, f.c_str(), t...);
    return std::string(&b[0], &b[0] + l);
}

template<typename T>
std::string stringify(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

/* dump */
#define ENABLE_DUMP
#ifdef ENABLE_DUMP
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
#else
#define dump(...) void(0);
#endif

/* timer */
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 3.0e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 3.0e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() { return (time() - t - paused) * 1000.0; }
} timer;

/* rand */
struct Xorshift {
    uint64_t x = 88172645463325252LL;
    void set_seed(unsigned seed, int rep = 100) { x = uint64_t((seed + 1) * 10007); for (int i = 0; i < rep; i++) next_int(); }
    unsigned next_int() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    unsigned next_int(unsigned mod) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % mod; }
    unsigned next_int(unsigned l, unsigned r) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % (r - l + 1) + l; } // inclusive
    double next_double() { return double(next_int()) / UINT_MAX; }
} rnd;

/* shuffle */
template<typename T>
void shuffle_vector(std::vector<T>& v, Xorshift& rnd) {
    int n = v.size();
    for (int i = n - 1; i >= 1; i--) {
        int r = rnd.next_int(i);
        std::swap(v[i], v[r]);
    }
}

/* split */
std::vector<std::string> split(std::string str, const std::string& delim) {
    for (char& c : str) if (delim.find(c) != std::string::npos) c = ' ';
    std::istringstream iss(str);
    std::vector<std::string> parsed;
    std::string buf;
    while (iss >> buf) parsed.push_back(buf);
    return parsed;
}

template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) {
    std::fill((T*)array, (T*)(array + N), val);
}

template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }



using pii = std::pair<int, int>;
using ll = long long;

constexpr int N = 50;
constexpr int NN = N * N;
constexpr int K = 8;

int A[N][N];

void init(std::istream& in) {
    int buf; in >> buf >> buf >> buf;
    std::vector<std::string> S(N);
    in >> S;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = S[i][j] - '0';
        }
    }
}

// https://qiita.com/ref3000/items/af18a4532123c22a19a4
struct PolyominoEnumeratorNaive {

    struct Cell;
    using CellPtr = std::shared_ptr<Cell>;
    struct Cell {
        int x, y, num;
        bool is_border;
        Cell(int x, int y, int num, bool is_border) : x(x), y(y), num(num), is_border(is_border) {}
        CellPtr clone() const {
            return std::make_shared<Cell>(x, y, num, is_border);
        }
    };

    struct Node;
    using NodePtr = std::shared_ptr<Node>;
    struct Node {
        std::vector<CellPtr> cells;
        std::deque<int> available_numbers;
        NodePtr clone() const {
            auto node = std::make_shared<Node>();
            for (auto cell : cells) {
                node->cells.push_back(cell->clone());
            }
            node->available_numbers = available_numbers;
            return node;
        }
    };

    int N;
    std::stack<NodePtr> node_stack;
    std::vector<std::vector<std::pair<int, int>>> poly;

    PolyominoEnumeratorNaive(int N) : N(N) {}

    void confirm_border_cell_number(NodePtr node, int num) {
        static constexpr int dx[] = { 0, -1, 1, 0 };
        static constexpr int dy[] = { -1, 0, 0, 1 };

        auto& cells = node->cells;
        auto cell = *std::find_if(cells.begin(), cells.end(), [&](CellPtr cell) { return cell->num == num; });
        cell->is_border = true;
        for (int i = 0; i < 4; i++) {
            int x = cell->x + dx[i];
            int y = cell->y + dy[i];
            if (!(y > 0 || (y == 0 && x >= 0))) continue;
            if (std::find_if(cells.begin(), cells.end(), [&](CellPtr cell) { return cell->x == x && cell->y == y; }) != cells.end()) continue;
            CellPtr new_cell = std::make_shared<Cell>(x, y, cells.size() + 1, false);
            cells.push_back(new_cell);
            node->available_numbers.push_back(cells.size());
        }
    }

    void next_step() {
        auto node = node_stack.top();
        if (node_stack.size() > N || node->available_numbers.empty()) {
            node_stack.pop();
        }
        else {
            auto num = node->available_numbers.front(); node->available_numbers.pop_front();
            auto new_node = node->clone();
            confirm_border_cell_number(new_node, num);
            node_stack.push(new_node);
        }
    }

    void run() {
        {
            NodePtr root = std::make_shared<Node>();
            CellPtr cell = std::make_shared<Cell>(0, 0, 1, false);
            root->cells.push_back(cell);
            root->available_numbers.push_back(1);
            node_stack.push(root);
        }
        int found = 0;
        while (!node_stack.empty()) {
            next_step();
            //if (!node_stack.empty()) show(node_stack.top(), 10);
            if (node_stack.size() > N) {
                ++found;
                poly.push_back({ {} });
                bool first = true;
                for (const auto& cell : node_stack.top()->cells) {
                    if (first) {
                        first = false;
                        continue;
                    }
                    if (cell->is_border) {
                        poly.back().emplace_back(cell->y, cell->x);
                    }
                }
            }
        }
        dump(found);
    }

    void show(NodePtr node, int delay = 0) const {
#ifdef HAVE_OPENCV_HIGHGUI
        int grid_size = 50;
        int H = N, W = 2 * N - 1;
        cv::Mat_<cv::Vec3b> img(H * grid_size, W * grid_size, cv::Vec3b(255, 255, 255));
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                cv::Rect roi(j * grid_size, i * grid_size, grid_size, grid_size);
                cv::rectangle(img, roi, cv::Scalar(0, 0, 0));
            }
        }
        for (auto cell : node->cells) {
            int i = cell->y, j = cell->x + N - 1;
            cv::Rect roi(j * grid_size + 1, i * grid_size + 1, grid_size - 2, grid_size - 2);
            if (cell->is_border) {
                cv::rectangle(img, roi, cv::Scalar(0, 255, 0), cv::FILLED);
            }
            cv::putText(img, std::to_string(cell->num), cv::Point(roi.x + grid_size / 2, roi.y + grid_size / 2), cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(0, 0, 0));
        }
        cv::imshow("img", img);
        cv::waitKey(delay);
#endif
    }
};

void output(std::ostream& out, const std::vector<pii>& path) {
    int npoints = path.size() / 8 * 8;
    out << npoints / 8 << '\n';
    for (int i = 0; i < npoints; i++) {
        out << path[i].first + 1 << ' ' << path[i].second + 1 << '\n';
    }
}

int main() {

#ifdef _MSC_VER
    std::ifstream ifs("C:\\dev\\heuristic\\tasks\\RCO2017qualA\\testcase\\in\\1.in");
    std::ofstream ofs("C:\\dev\\heuristic\\tasks\\RCO2017qualA\\testcase\\out\\1.out");
    std::istream& in = ifs;
    std::ostream& out = ofs;
#else
    std::istream& in = std::cin;
    std::ostream& out = std::cout;
#endif

    init(in);

    PolyominoEnumeratorNaive polyenum(K);
    polyenum.run();
    auto polys = polyenum.poly;

    auto calc_score = [&](const std::vector<std::pair<int, int>>& poly, int y, int x) {
        if (!A[y][x]) return 0;
        int prod = 1;
        for (auto [dy, dx] : poly) {
            int ny = y + dy, nx = x + dx;
            if (ny < 0 || ny >= N || nx < 0 || nx >= N || !A[ny][nx]) return 0;
            prod *= A[ny][nx];
        }
        return prod;
    };

    std::vector<std::pair<int, int>> ans;
    ll total_score = 0;
    while (timer.elapsed_ms() < 9800) {
        int best_score = 0;
        int best_y = -1, best_x = -1, best_pid = -1;
        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                if (!A[y][x]) continue;
                for (int i = 0; i < polys.size(); i++) {
                    int score = calc_score(polys[i], y, x);
                    if (best_score < score) {
                        best_score = score;
                        best_y = y;
                        best_x = x;
                        best_pid = i;
                    }
                }
            }
        }
        if (best_pid == -1) break;
        for (auto [dy, dx] : polys[best_pid]) {
            int y = best_y + dy, x = best_x + dx;
            A[y][x] = 0;
            ans.emplace_back(y + 1, x + 1);
        }
        total_score += best_score;
        dump(timer.elapsed_ms(), total_score, best_score, best_y, best_x, polys[best_pid]);
    }

    out << ans.size() / 8 << '\n';
    for (auto [y, x] : ans) {
        out << y << ' ' << x << '\n';
    }

    return 0;
}