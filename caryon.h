#ifndef CARYON_H
#define CARYON_H

#include <cstring>
#include <algorithm>
#include <complex>
#include <functional>
#include <iostream>
#include <numeric>
#include <queue>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <fstream>
#include <mutex>
#include <thread>
#include <filesystem>
#include <ctime>
#include <cstdint>

#define makein(low, high) for (caseIndex = low, invocationCount = 0; caseIndex <= high; caseIndex++, invocationCount = 0)

#define CarYon 4

// 全局变量注释
// datasetName: 用于生成测试数据时的基准数据集名称。
// outputDirectory: 指定输出文件夹路径。
// caseIndex: 当前要写入的测试用例编号（全局上下文）。
// invocationCount: 同一用例内部被写入的次数计数器。
// maxRuntimeMs: 判定 TLE 的全局最大运行时间限制（毫秒）。
// runtimeMs: 记录当前运行时间（用于调试比较）。
// directoryCreated: 标记输出目录是否已创建过。
// runtimeLogStream: 运行时间日志流（调试用）。
// caryon_io_mutex: IO 操作互斥锁，保证多线程写入安全。
std::string datasetName;
std::string outputDirectory;
int caseIndex;
int invocationCount;
int maxRuntimeMs = 1000;
long double runtimeMs;
bool directoryCreated = false;
std::stringstream runtimeLogStream;
std::mutex caryon_io_mutex;

// CaryonConstants: 常用常量
namespace CaryonConstants {
	const long double PI = 3.141592653589793238462643383279502884197169399;
	const long double E = 2.7182818284590452353602874713527;
	const long double PHI = 1.61803398874989484820458683436563811772030917980576;
	const long double SQRT2 = 1.4142135623730950488016887242096980785696718753769;
	const long double SQRT3 = 1.73205080756887729352744634150587236694280525381038;
	const long double SQRT5 = 2.23606797749978969640917366873127623544061835961153;
	const long double LOG2 = 0.693147180559945309417232121458176568075500134360255254120680009;
	const long double LOG10 = 2.3025850929940456840179914546843642076011014886287729760333279095078;
	const std::string ALPHABET_SMALL = "abcdefghijklmnopqrstuvwxyz";
	const std::string ALPHABET_CAPITAL = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	const std::string NUMBERS = "0123456789";
} // namespace CaryonConstants

// CaryonRandom: 随机数据工具集合
namespace CaryonRandom {
	// randomInitialized: 标记 RNG 是否已初始化
	// useStdMT: 如果为 true 则使用 std::mt19937_64，否则使用自定义“legacy”实现
	bool randomInitialized = false;
	bool useStdMT = true;
	int legacyIndex;
	long long legacyMT[624];
	inline std::mt19937_64 rng64;

	/**
	 * seedStd
	 * - 描述: 使用 std::mt19937_64 对随机数生成器进行种子初始化。
	 * - 参数:
	 *    seedValue: 用于初始化标准 MT 引擎的无符号长整型种子。
	 */
	inline void seedStd(unsigned long long seedValue) {
		rng64.seed(seedValue);
		useStdMT = true;
		randomInitialized = true;
	}

	/**
	 * seedLegacy
	 * - 描述: 初始化 legacy（自定义）MT 实现的种子序列。
	 * - 参数:
	 *    initialSeed: 作为 legacyMT[0] 的初始种子（32 位截断）。
	 */
	void seedLegacy(int initialSeed) {
		legacyIndex = 0;
		randomInitialized = 1;
		legacyMT[0] = initialSeed;
		for (int i = 1; i < 624; i++) {
			int t = 1812433253 * (legacyMT[i - 1] ^ (legacyMT[i - 1] >> 30)) + i;
			legacyMT[i] = t & 0xffffffff;
		}
	}

	// legacyGenerate: 生成 legacy MT 的下一批伪随机数（内部使用）
	inline void legacyGenerate() {
		for (int i = 0; i < 624; i++) {
			long long y = (legacyMT[i] & 0x80000000) + (legacyMT[(i + 1) % 624] & 0x7fffffff);
			legacyMT[i] = legacyMT[(i + 397) % 624] ^ (y >> 1);
			if (y % 2 == 1)
				legacyMT[i] ^= 2147483647;
		}
	}

	/**
	 * randInt
	 * - 描述: 返回一个非负的 32 位随机整数。
	 * - 返回: int 随机值范围 [0, 2^31-1)
	 */
	inline int randInt() {
		if (useStdMT) {
			return (int)(rng64() & 0x7fffffff);
		}
		if (!randomInitialized)
			seedLegacy((int)time(nullptr));
		if (legacyIndex == 0)
			legacyGenerate();
		int y = legacyMT[legacyIndex];
		y = y ^ (y >> 11);
		y = y ^ ((y << 7) & 1636928640);
		y = y ^ ((y << 15) & 1022730752);
		y = y ^ (y >> 18);
		legacyIndex = (legacyIndex + 1) % 624;
		return y;
	}

	/**
	 * randLong
	 * - 描述: 返回一个随机的 64 位有符号整数（取自 mt19937_64）。
	 */
	inline long long randLong() {
		if (useStdMT) {
			return (long long)rng64();
		}
		return ((long long)(randInt()) << 31) + randInt();
	}

	/**
	 * randInt
	 * - 描述: 返回区间 [low, high] 之间的随机整数（含端点）。
	 * - 参数:
	 *    low: 下界（包含）
	 *    high: 上界（包含）
	 */
	inline int randInt(int low, int high) {
		if (low > high)
			low = high;
		if (low == high)
			return low;
		if (useStdMT) {
			uint64_t r = rng64();
			return (int)(r % (uint64_t)(high - low + 1) + low);
		}
		else
			return randInt() % (high - low + 1) + low;
	}

	/**
	 * randLong
	 * - 描述: 返回区间 [low, high] 之间的随机长整型（含端点）。
	 * - 参数:
	 *    low: 下界（包含）
	 *    high: 上界（包含）
	 */
	inline long long randLong(long long low, long long high) {
		if (low > high)
			low = high;
		if (low == high)
			return low;
		if (useStdMT) {
			uint64_t r = rng64();
			return (long long)(r % (uint64_t)(high - low + 1) + low);
		}
		else
			return randLong() % (high - low + 1) + low;
	}

	/**
	 * shuffleRange
	 * - 描述: 使用当前 RNG 对 [first, last) 范围内的元素进行 Fisher-Yates 洗牌。
	 * - 模板参数:
	 *    It: 迭代器类型
	 * - 参数:
	 *    first: 区间起始迭代器（包含）
	 *    last: 区间结束迭代器（不包含）
	 */
	template < typename It >
	inline void shuffleRange(It first, It last) {
		if (useStdMT) {
			std::shuffle(first, last, rng64);
		}
		else {
			auto n = std::distance(first, last);
			for (decltype(n) i = n - 1; i > 0; --i) {
				auto j = randInt(0, (int)i);
				using std::iter_swap;
				iter_swap(first + i, first + j);
			}
		}
	}

	/**
	 * randPermutation
	 * - 描述: 生成 1..length 的随机排列并写入 outArray（长度为 length）。
	 * - 参数:
	 *    outArray: 指向输出数组的指针（需至少 capacity >= length）
	 *    length: 要生成的排列长度
	 */
	inline void randPermutation(int *outArray, int length) {
		for (int i = 0; i < length; i++) {
			outArray[i] = i + 1;
		}
		for (int i = length - 1; i > 0; --i) {
			int j = randInt(0, i);
			std::swap(outArray[i], outArray[j]);
		}
	}

	/**
	 * randDouble
	 * - 描述: 返回 [0.0, 1.0) 区间的双精度随机数。
	 */
	inline double randDouble() {
		if (useStdMT) {
			return std::uniform_real_distribution<double>(0.0, 1.0)(rng64);
		}
		return (double)(randInt()) / (double)0x7fffffff;
	}

	/**
	 * seedAllRNG
	 * - 描述: 使用相同的种子同时初始化标准 MT 与 legacy MT。
	 * - 参数:
	 *    seedValue: 要使用的种子
	 */
	inline void seedAllRNG(unsigned long long seedValue) {
		seedStd(seedValue);
		legacyIndex = 0;
		legacyMT[0] = (long long)(seedValue & 0xffffffff);
		for (int i = 1; i < 624; ++i) {
			int t = 1812433253 * (legacyMT[i - 1] ^ (legacyMT[i - 1] >> 30)) + i;
			legacyMT[i] = t & 0xffffffff;
		}
		randomInitialized = true;
	}

	/**
	 * randVector
	 * - 描述: 生成包含 count 个元素的随机整数向量，元素在 [minValue, maxValue]。
	 * - 参数:
	 *    count: 向量长度
	 *    minValue: 元素下界（包含）
	 *    maxValue: 元素上界（包含）
	 */
	inline std::vector<int> randVector(int count, int minValue, int maxValue) {
		std::vector<int> v;
		v.reserve(count);
		for (int i = 0; i < count; ++i) v.push_back(randInt(minValue, maxValue));
		return v;
	}

	/**
	 * sampleUniqueRange
	 * - 描述: 从 [minValue, maxValue] 中随机抽取 count 个互不相同的数（无序返回）。
	 * - 参数:
	 *    minValue: 范围下界（包含）
	 *    maxValue: 范围上界（包含）
	 *    count: 要抽取的唯一元素数量
	 */
	inline std::vector<int> sampleUniqueRange(int minValue, int maxValue, int count) {
		int range = maxValue - minValue + 1;
		if (count <= 0) return {};
		if (count >= range) {
			std::vector<int> all(range);
			for (int i = 0; i < range; ++i) all[i] = minValue + i;
			return all;
		}
		std::vector<int> pool(range);
		for (int i = 0; i < range; ++i) pool[i] = minValue + i;
		for (int i = 0; i < count; ++i) {
			int j = randInt(i, range - 1);
			std::swap(pool[i], pool[j]);
		}
		return std::vector<int>(pool.begin(), pool.begin() + count);
	}

	/**
	 * randString
	 * - 描述: 生成指定长度的随机字符串。
	 * - 参数:
	 *    length: 字符串长度
	 *    charset: 字符集（默认小写字母）
	 * - 返回: 随机字符串
	 */
	inline std::string randString(int length, const std::string &charset = CaryonConstants::ALPHABET_SMALL) {
		std::string result;
		result.reserve(length);
		for (int i = 0; i < length; ++i) {
			result += charset[CaryonRandom::randInt(0, charset.size() - 1)];
		}
		return result;
	}

	/**
	 * randBinaryString
	 * - 描述: 生成指定长度的随机二进制字符串（仅含 '0' 和 '1'）。
	 * - 参数:
	 *    length: 字符串长度
	 * - 返回: 二进制字符串
	 */
	inline std::string randBinaryString(int length) {
		return randString(length, "01");
	}

	/**
	 * randPalindrome
	 * - 描述: 生成指定长度的随机回文串。
	 * - 参数:
	 *    length: 回文串长度
	 *    charset: 字符集
	 * - 返回: 回文串
	 */
	inline std::string randPalindrome(int length, const std::string &charset = CaryonConstants::ALPHABET_SMALL) {
		std::string result(length, ' ');
		int half = length / 2;
		for (int i = 0; i < half; ++i) {
			char c = charset[CaryonRandom::randInt(0, charset.size() - 1)];
			result[i] = c;
			result[length - 1 - i] = c;
		}
		if (length % 2 == 1) {
			result[half] = charset[CaryonRandom::randInt(0, charset.size() - 1)];
		}
		return result;
	}

	/**
	 * shuffleString
	 * - 描述: 随机打乱字符串。
	 * - 参数:
	 *    s: 输入字符串
	 * - 返回: 打乱后的字符串
	 */
	inline std::string shuffleString(const std::string &s) {
		std::vector<char> chars(s.begin(), s.end());
		CaryonRandom::shuffleRange(chars.begin(), chars.end());
		return std::string(chars.begin(), chars.end());
	}
} // namespace CaryonRandom

// CaryonGraph: 基本图结构与生成器
namespace CaryonGraph {
	/**
	 * Edge<T>
	 * - 描述: 图的边表示，包含目标顶点 v 与权重 w。
	 * - 成员:
	 *    v: 目标顶点编号
	 *    w: 边的权重（模板类型）
	 */
	template <typename T>
	struct Edge {
		int v;
		T w;
		bool operator<(const Edge &rw) const {
			return w > rw.w;
		}
	};

	/**
	 * Adjacency<T>
	 * - 描述: 邻接表项，包含若干边。
	 */
	template <typename T>
	struct Adjacency {
		std::vector< Edge<T> > edges;
	};

	/**
	 * Graph<T>
	 * - 描述: 简单的无向/有向图邻接表实现（顶点从 1 开始）。
	 * - 成员:
	 *    n: 当前最大顶点编号
	 *    m: 边数
	 *    adj: 邻接表（按索引存储 Adjacency）
	 */
	template <typename T>
	struct Graph {
		int n = 0, m = 0;
		std::vector< Adjacency<T> > adj;
		Graph() {}

		/**
		 * ensureSize
		 * - 描述: 确保 adj 的大小至少能包含 targetIndex（按顶点编号索引）。
		 * - 参数:
		 *    targetIndex: 需要确保的最大索引
		 */
		void ensureSize(int targetIndex) {
			Adjacency<T> updatemp;
			updatemp.edges.clear();
			while (adj.size() <= targetIndex) adj.push_back(updatemp);
		}

		/**
		 * addEdge
		 * - 描述: 向图中添加一条从 from 到 to 的边，权重为 weight。
		 * - 参数:
		 *    from: 边的起始顶点编号
		 *    to: 边的目标顶点编号
		 *    weight: 边的权重
		 */
		void addEdge(int from, int to, T weight) {
			n = std::max(n, from);
			n = std::max(n, to);
			ensureSize(std::max(from, to));
			m++;
			Edge<T> tmp;
			tmp.v = to;
			tmp.w = weight;
			adj[from].edges.push_back(tmp);
		}

		/**
		 * isConnected
		 * - 描述: 简单 BFS 连通性检测（从随机顶点开始）。
		 * - 返回: 若图的所有顶点都可达则返回 true。
		 */
		bool isConnected() {
			const int vstsize = n + 1;
			int vstn = 0;
			std::vector<bool> visited(vstsize, false);
			std::queue<int> q;
			int start = CaryonRandom::randInt(1, n);
			visited[start] = true;
			vstn = 1;
			q.push(start);
			while (!q.empty()) {
				int cur = q.front();
				q.pop();
				for (int i = 0; i < adj[cur].edges.size(); i++) {
					if (!visited[adj[cur].edges[i].v]) {
						visited[adj[cur].edges[i].v] = true;
						vstn++;
						q.push(adj[cur].edges[i].v);
					}
				}
			}
			return (vstn == n);
		}
	};

	template <typename T>
	std::ostream &operator<<(std::ostream &os, const Graph<T> &c) {
		os << c.n << ' ' << c.m << '\n';
		for (int i = 1; i <= c.n; i++) {
			for (int j = 0; j < c.adj[i].edges.size(); j++)
				os << i << ' ' << c.adj[i].edges[j].v << ' ' << c.adj[i].edges[j].w << '\n';
		}
		return os;
	}

	template <typename T>
	Graph<T> operator+(Graph<T> a, Graph<T> b) {
		Graph<T> ret;
		ret = a;
		ret.m = a.m + b.m;
		ret.ensureSize(ret.n);
		for (int i = 1; i <= b.n; i++) {
			for (int j = 0; j < b.adj[i].edges.size(); j++) {
				ret.adj[i].edges.push_back(b.adj[i].edges[j]);
			}
		}
		return ret;
	}

	/**
	 * randGraph
	 * - 描述: 随机生成具有 nodeCount 个顶点和 edgeCount 条边的图。
	 * - 参数:
	 *    nodeCount: 顶点数量
	 *    edgeCount: 边数量
	 *    minWeight: 权重最小值
	 *    maxWeight: 权重最大值
	 *    randValueFunc: 从 [minWeight, maxWeight] 返回一个随机权重的函数指针
	 */
	template <typename T>
	Graph<T> randGraph(int nodeCount, int edgeCount, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		Graph<T> ret;
		ret.n = nodeCount;
		ret.m = edgeCount;
		ret.ensureSize(nodeCount);
		for (int i = 1; i <= edgeCount; i++) {
			Edge<T> tmp;
			tmp.v = CaryonRandom::randInt(1, nodeCount);
			tmp.w = randValueFunc(minWeight, maxWeight);
			ret.adj[CaryonRandom::randInt(1, nodeCount)].edges.push_back(tmp);
		}
		return ret;
	}

	/**
	 * randDAG
	 * - 描述: 随机生成有向无环图（DAG）。
	 * - 参数:
	 *    nodeCount: 顶点数量
	 *    edgeCount: 边数量
	 *    minWeight: 权重最小值
	 *    maxWeight: 权重最大值
	 *    randValueFunc: 随机权重生成函数
	 */
	template <typename T>
	Graph<T> randDAG(int nodeCount, int edgeCount, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		Graph<T> ret;
		ret.ensureSize(nodeCount);
		ret.n = nodeCount;
		ret.m = edgeCount;
		for (int i = 1; i <= edgeCount; i++) {
			Edge< T > tmp;
			int utmp = CaryonRandom::randInt(1, CaryonRandom::randInt(1, nodeCount - 1));
			tmp.v = CaryonRandom::randInt(utmp + 1, nodeCount);
			tmp.w = randValueFunc(minWeight, maxWeight);
			ret.adj[utmp].edges.push_back(tmp);
		}
		return ret;
	}

	struct CaryonNode {
		int id;
		int sonCount;
	};

	/**
	 * randTree
	 * - 描述: 生成一个随机树（有向结构或父子关系类似），使用最大子节点数控制形状。
	 * - 参数:
	 *    nodeCount: 节点数
	 *    maxChildren: 每个节点最大子节点数（度约束）
	 *    minWeight: 权重最小值
	 *    maxWeight: 权重最大值
	 *    randValueFunc: 随机权重生成函数
	 */
	template < typename T >
	Graph<T> randTree(int nodeCount, int maxChildren, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		Graph<T> ret;
		ret.n = nodeCount;
		ret.m = nodeCount - 1;
		ret.ensureSize(nodeCount);
		std::vector<CaryonNode> t;
		std::vector<int> nodes;
		for (int i = 1; i <= nodeCount; i++) nodes.push_back(i);
		CaryonRandom::shuffleRange(nodes.begin(), nodes.end());
		CaryonNode updatemp;
		updatemp.id = nodes[0];
		updatemp.sonCount = 0;
		t.push_back(updatemp);
		for (int j = 2; j <= nodeCount; j++) {
			int i = nodes[j - 1];
			std::swap(t[CaryonRandom::randInt(0, t.size() - 1)], t[t.size() - 1]);
			t[t.size() - 1].sonCount++;
			if (t[t.size() - 1].sonCount == maxChildren) t.pop_back();
			Edge<T> tmp;
			tmp.v = t[t.size() - 1].id;
			tmp.w = randValueFunc(minWeight, maxWeight);
			ret.adj[i].edges.push_back(tmp);
			updatemp.id = i;
			updatemp.sonCount = 0;
			t.push_back(updatemp);
		}
		return ret;
	}

	/**
	 * completeGraph
	 * - 描述: 生成完全图（每对顶点之间都有边）。
	 * - 参数:
	 *    nodeCount: 顶点数
	 *    minWeight: 最小权重
	 *    maxWeight: 最大权重
	 *    randValueFunc: 权重生成函数
	 * - 返回: 完全图
	 */
	template<typename T>
	CaryonGraph::Graph<T> completeGraph(int nodeCount, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		CaryonGraph::Graph<T> graph;
		graph.n = nodeCount;
		graph.ensureSize(nodeCount);
		for (int i = 1; i <= nodeCount; ++i) {
			for (int j = i + 1; j <= nodeCount; ++j) {
				graph.addEdge(i, j, randValueFunc(minWeight, maxWeight));
				graph.addEdge(j, i, randValueFunc(minWeight, maxWeight));
			}
		}
		return graph;
	}

	/**
	 * bipartiteGraph
	 * - 描述: 生成随机二分图。
	 * - 参数:
	 *    leftSize: 左部顶点数
	 *    rightSize: 右部顶点数
	 *    edgeCount: 边数
	 *    minWeight: 最小权重
	 *    maxWeight: 最大权重
	 *    randValueFunc: 权重生成函数
	 * - 返回: 二分图
	 */
	template<typename T>
	CaryonGraph::Graph<T> bipartiteGraph(int leftSize, int rightSize, int edgeCount, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		CaryonGraph::Graph<T> graph;
		graph.n = leftSize + rightSize;
		graph.ensureSize(graph.n);
		std::vector<std::pair<int, int>> edges;
		for (int i = 1; i <= leftSize; ++i) {
			for (int j = leftSize + 1; j <= leftSize + rightSize; ++j) {
				edges.emplace_back(i, j);
			}
		}
		CaryonRandom::shuffleRange(edges.begin(), edges.end());
		if (edgeCount > edges.size()) edgeCount = edges.size();
		for (int i = 0; i < edgeCount; ++i) {
			auto [u, v] = edges[i];
			graph.addEdge(u, v, randValueFunc(minWeight, maxWeight));
		}
		return graph;
	}

	/**
	 * gridGraph
	 * - 描述: 生成网格图（二维网格，相邻格子有边）。
	 * - 参数:
	 *    rows: 行数
	 *    cols: 列数
	 *    minWeight: 最小权重
	 *    maxWeight: 最大权重
	 *    randValueFunc: 权重生成函数
	 * - 返回: 网格图
	 */
	template<typename T>
	CaryonGraph::Graph<T> gridGraph(int rows, int cols, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		CaryonGraph::Graph<T> graph;
		graph.n = rows * cols;
		graph.ensureSize(graph.n);
		auto id = [&](int r, int c) { return r * cols + c + 1; };
		for (int r = 0; r < rows; ++r) {
			for (int c = 0; c < cols; ++c) {
				if (r + 1 < rows) {
					graph.addEdge(id(r, c), id(r + 1, c), randValueFunc(minWeight, maxWeight));
					graph.addEdge(id(r + 1, c), id(r, c), randValueFunc(minWeight, maxWeight));
				}
				if (c + 1 < cols) {
					graph.addEdge(id(r, c), id(r, c + 1), randValueFunc(minWeight, maxWeight));
					graph.addEdge(id(r, c + 1), id(r, c), randValueFunc(minWeight, maxWeight));
				}
			}
		}
		return graph;
	}

	/**
	 * connectedGraph
	 * - 描述: 生成保证连通的随机图（先生成树，再随机加边）。
	 * - 参数:
	 *    nodeCount: 顶点数
	 *    edgeCount: 总边数（至少 nodeCount-1）
	 *    minWeight: 最小权重
	 *    maxWeight: 最大权重
	 *    randValueFunc: 权重生成函数
	 * - 返回: 连通图
	 */
	template<typename T>
	CaryonGraph::Graph<T> connectedGraph(int nodeCount, int edgeCount, T minWeight, T maxWeight, T(*randValueFunc)(T, T)) {
		if (edgeCount < nodeCount - 1) edgeCount = nodeCount - 1;
		auto tree = CaryonGraph::randTree(nodeCount, nodeCount - 1, minWeight, maxWeight, randValueFunc);
		int remainingEdges = edgeCount - (nodeCount - 1);
		for (int i = 0; i < remainingEdges; ++i) {
			int u = CaryonRandom::randInt(1, nodeCount);
			int v = CaryonRandom::randInt(1, nodeCount);
			if (u == v) continue;
			tree.addEdge(u, v, randValueFunc(minWeight, maxWeight));
		}
		return tree;
	}

	// 生成链状树
	template<typename T>
	Graph<T> randChain(int nodeCount, T minWeight, T maxWeight, T(*randFunc)(T, T)) {
		Graph<T> ret;
		ret.n = nodeCount;
		ret.m = nodeCount - 1;
		ret.ensureSize(nodeCount);
		for (int i = 2; i <= nodeCount; ++i) {
			ret.addEdge(i - 1, i, randFunc(minWeight, maxWeight));
		}
		return ret;
	}

	// 生成菊花图
	template<typename T>
	Graph<T> randStar(int nodeCount, T minWeight, T maxWeight, T(*randFunc)(T, T)) {
		Graph<T> ret;
		ret.n = nodeCount;
		ret.m = nodeCount - 1;
		ret.ensureSize(nodeCount);
		for (int i = 2; i <= nodeCount; ++i) {
			ret.addEdge(1, i, randFunc(minWeight, maxWeight));
		}
		return ret;
	}

	/**
	 * isTree
	 * - 描述: 验证图是否为树（连通且边数 = 顶点数 - 1）。
	 * - 参数:
	 *    graph: 要验证的图对象
	 * - 返回: 如果是树则返回 true
	 */
	template<typename T>
	bool isTree(const CaryonGraph::Graph<T> &graph) {
		return graph.isConnected() && graph.m == graph.n - 1;
	}

	/**
	 * isDAG
	 * - 描述: 使用拓扑排序验证图是否为有向无环图。
	 * - 参数:
	 *    graph: 要验证的图对象
	 * - 返回: 如果是 DAG 则返回 true
	 */
	template<typename T>
	bool isDAG(const CaryonGraph::Graph<T> &graph) {
		std::vector<int> indegree(graph.n + 1, 0);
		for (int i = 1; i <= graph.n; ++i) {
			for (const auto &e : graph.adj[i].edges) {
				indegree[e.v]++;
			}
		}
		std::queue<int> q;
		for (int i = 1; i <= graph.n; ++i) {
			if (indegree[i] == 0) q.push(i);
		}
		int cnt = 0;
		while (!q.empty()) {
			int u = q.front(); q.pop();
			cnt++;
			for (const auto &e : graph.adj[u].edges) {
				if (--indegree[e.v] == 0) {
					q.push(e.v);
				}
			}
		}
		return cnt == graph.n;
	}

	/**
	 * hasDuplicateEdges
	 * - 描述: 检测图中是否存在重复边（相同起点和终点）。
	 * - 参数:
	 *    graph: 要检测的图对象
	 * - 返回: 如果存在重复边则返回 true
	 */
	template<typename T>
	bool hasDuplicateEdges(const CaryonGraph::Graph<T> &graph) {
		for (int i = 1; i <= graph.n; ++i) {
			std::unordered_set<int> targets;
			for (const auto &e : graph.adj[i].edges) {
				if (targets.count(e.v)) return true;
				targets.insert(e.v);
			}
		}
		return false;
	}
} // namespace CaryonGraph

// CaryonIO: 数据集/用例生成与运行相关的辅助工具
namespace CaryonIO {
	/**
	 * writeGraphCase
	 * - 描述: 将图写入指定的测试用例文件。
	 * - 参数:
	 *    graph: 要写入的图对象
	 */
	template <typename T>
	inline void writeGraphCase(CaryonGraph::Graph<T> graph) {
		std::stringstream cci;
		std::string tnmp;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << caseIndex;
			std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
			freopen(name.c_str(), "w", stdout);
			freopen(name.c_str(), "r", stdin);
		}
		graph.print(graph);
		++invocationCount;
	}

	/**
	 * writeCase
	 * - 描述: 将任意可流输出的值写入当前测试用例。
	 * - 参数:
	 *    value: 要写入的值（支持输出流输出）
	 */
	template <typename T>
	inline void writeCase(const T &value) {
		std::string c, tnmp;
		std::stringstream ss, cci;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << caseIndex;
			std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
			freopen(name.c_str(), "w", stdout);
			freopen(name.c_str(), "r", stdin);
		}
		std::cout << value;
		++invocationCount;
	}

	/**
	 * writeSpace
	 * - 描述: 将一个空格符写入当前测试用例。
	 */
	inline void writeSpace() {
		std::stringstream cci;
		std::string tnmp;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << caseIndex;
			std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
			freopen(name.c_str(), "w", stdout);
			freopen(name.c_str(), "r", stdin);
		}
		std::cout << " ";
		++invocationCount;
	}

	/**
	 * writeEndl
	 * - 描述: 将一个换行符写入当前测试用例。
	 */
	inline void writeEndl() {
		std::stringstream cci;
		std::string tnmp;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << caseIndex;
			std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
			freopen(name.c_str(), "w", stdout);
			freopen(name.c_str(), "r", stdin);
		}
		std::cout << "\n";
		++invocationCount;
	}

	/**
	 * executeStd
	 * - 描述: 在命令行上运行指定编号的测试用例对应的可执行程序，并将输入输出重定向到文件。
	 * - 参数:
	 *    caseNumber: 要执行的用例编号
	 */
	inline void executeStd(int caseNumber) {
		freopen("CON.exe", "w", stdout);
		freopen("CON.exe", "r", stdin);
		std::stringstream aa;
		std::string aaa;
		aa << caseNumber;
		aa >> aaa;
		std::string name = "data-" + datasetName + "/" + datasetName + aaa + ".in";
		freopen(name.c_str(), "r", stdin);
		std::string outname = "data-" + datasetName + "/" + datasetName + aaa + ".out";
		freopen(outname.c_str(), "w", stdout);
		system("std.exe");
	}

	/**
	 * executeRangeStd
	 * - 描述: 对一系列 case 从 startIndex 到 endIndex 执行 executeStd。
	 * - 参数:
	 *    startIndex: 起始用例编号
	 *    endIndex: 结束用例编号
	 */
	inline void executeRangeStd(int startIndex, int endIndex) {
		for (int i = startIndex; i <= endIndex; ++i) {
			executeStd(i);
			freopen("CON.exe", "r", stdin);
			freopen("CON.exe", "w", stdout);
			std::stringstream _a;
			std::string _i;
			_a << i;
			_a >> _i;
		}
	}

	inline void closeStreams() {
		freopen("CON.exe", "w", stdout);
		freopen("CON.exe", "r", stdin);
	}

	class CaseFileWriter {
	private:
		std::ofstream ofs;
		int caseIndexLocal;
		std::string baseName;
	public:
		CaseFileWriter() : caseIndexLocal(0) {}
		void setBaseName(const std::string &name) { baseName = name; }
		void setOutputDir(const std::string &dir) { outputDirectory = dir; }
		bool openCase(int idx) {
			caseIndexLocal = idx;
			try {
				if (outputDirectory.empty()) {
					std::string dir = "data-" + datasetName;
					std::filesystem::create_directories(dir);
					outputDirectory = dir;
				}
				std::string fname = outputDirectory + "/" + (baseName.empty() ? datasetName : baseName) + std::to_string(idx) + ".in";
				ofs.open(fname, std::ios::out | std::ios::trunc);
				return ofs.is_open();
			} catch (...) {
				return false;
			}
		}
		template <typename T>
		void write(const T &v) {
			std::lock_guard<std::mutex> lk(caryon_io_mutex);
			ofs << v;
		}
		void writeln(const std::string &s) {
			std::lock_guard<std::mutex> lk(caryon_io_mutex);
			ofs << s << '\n';
		}
		void close() {
			std::lock_guard<std::mutex> lk(caryon_io_mutex);
			if (ofs.is_open()) ofs.close();
		}
	};

	/**
	 * setDatasetName
	 * - 描述: 设置全局数据集名称（用于生成文件夹与文件名）。
	 * - 参数:
	 *    name: 数据集名称字符串
	 */
	inline void setDatasetName(const std::string &name) { datasetName = name; }

	/**
	 * setOutputDirectory
	 * - 描述: 设置输出目录（并尝试创建该目录）。
	 * - 参数:
	 *    dir: 输出目录路径
	 */
	inline void setOutputDirectory(const std::string &dir) { outputDirectory = dir; std::filesystem::create_directories(outputDirectory); }

	/**
	 * setMaxRuntimeGlobal
	 * - 描述: 设置判定 TLE 使用的全局最大运行时间（毫秒）。
	 * - 参数:
	 *    t: 毫秒数
	 */
	inline void setMaxRuntimeGlobal(int t) { maxRuntimeMs = t; }

	/**
	 * batchGenerateCases
	 * - 描述: 多线程批量生成测试用例文件，每个线程负责一段编号区间。
	 * - 参数:
	 *    cases: 总共要生成的用例数量（1...cases）
	 *    generator: 功能对象，签名为 (int caseIndex, CaseFileWriter &writer)
	 *    threads: 并发线程数（默认 1）
	 *    baseName: 基础文件名（默认使用 datasetName）
	 */
	inline void batchGenerateCases(int cases, std::function<void(int, CaseFileWriter &)> generator, int threads = 1, const std::string &baseName = "") {
		if (cases <= 0) return;
		if (threads <= 0) threads = 1;
		std::filesystem::path dir = outputDirectory.empty() ? std::filesystem::path("data-" + datasetName) : std::filesystem::path(outputDirectory);
		std::filesystem::create_directories(dir);
		auto worker = [&](int from, int to) {
			for (int i = from; i <= to; ++i) {
				CaseFileWriter writer;
				writer.setBaseName(baseName.empty() ? datasetName : baseName);
				writer.setOutputDir(dir.string());
				if (!writer.openCase(i)) continue;
				generator(i, writer);
				writer.close();
			}
			};
		int per = cases / threads;
		int rem = cases % threads;
		std::vector<std::thread> ths;
		int cur = 1;
		for (int t = 0; t < threads; ++t) {
			int start = cur;
			int cnt = per + (t < rem ? 1 : 0);
			int end = start + cnt - 1;
			if (cnt <= 0) break;
			ths.emplace_back(worker, start, end);
			cur = end + 1;
		}
		for (auto &t : ths) if (t.joinable()) t.join();
	}

	inline std::string getCaseInputPath(int idx, const std::string &baseName = "") {
		std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
		std::string name = dir + "/" + (baseName.empty() ? datasetName : baseName) + std::to_string(idx) + ".in";
		return name;
	}

	inline int runExecutableOnCase(int idx, const std::string &exeName = "test.exe", const std::string &baseName = "") {
		std::string inPath = getCaseInputPath(idx, baseName);
		std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
		std::string outPath = dir + "/" + (baseName.empty() ? datasetName : baseName) + std::to_string(idx) + ".out";
		std::string command = "\"" + exeName + "\" < \"" + inPath + "\" > \"" + outPath + "\"";
		return system(command.c_str());
	}

	inline void runExecutablesParallel(int start, int end, const std::string &exeName = "test.exe", int threads = 1, const std::string &baseName = "") {
		if (start > end) return;
		if (threads <= 0) threads = 1;
		int total = end - start + 1;
		int per = total / threads;
		int rem = total % threads;
		std::vector<std::thread> ths;
		int cur = start;
		auto worker = [&](int s, int e) {
			for (int i = s; i <= e; ++i) {
				runExecutableOnCase(i, exeName, baseName);
			}
			};
		for (int t = 0; t < threads; ++t) {
			int cnt = per + (t < rem ? 1 : 0);
			if (cnt <= 0) break;
			int s = cur;
			int e = cur + cnt - 1;
			ths.emplace_back(worker, s, e);
			cur = e + 1;
		}
		for (auto &th : ths) if (th.joinable()) th.join();
	}

	inline void generateAndRunCases(int cases, std::function<void(int, CaseFileWriter &)> generator, int genThreads = 1, int runThreads = 1, const std::string &exeName = "test.exe", const std::string &baseName = "") {
		batchGenerateCases(cases, generator, genThreads, baseName);
		runExecutablesParallel(1, cases, exeName, runThreads, baseName);
	}

	class FileCaseWriter {
	private:
		std::ofstream ofs;
		int caseIndex;
		std::string baseName;
	public:
		FileCaseWriter() : caseIndex(0) {}
		void setBaseName(const std::string &name) { baseName = name; }
		void setOutputDir(const std::string &dir) { outputDirectory = dir; }
		bool openCase(int idx) {
			caseIndex = idx;
			try {
				if (outputDirectory.empty()) {
					std::string dir = "data-" + datasetName;
					std::filesystem::create_directories(dir);
					outputDirectory = dir;
				}
				std::string fname = outputDirectory + "/" + (baseName.empty() ? datasetName : baseName) + std::to_string(idx) + ".in";
				ofs.open(fname, std::ios::out | std::ios::trunc);
				return ofs.is_open();
			} catch (...) {
				return false;
			}
		}
		template <typename T>
		void write(const T &v) {
			std::lock_guard<std::mutex> lk(caryon_io_mutex);
			ofs << v;
		}
		void writeln(const std::string &s) {
			std::lock_guard<std::mutex> lk(caryon_io_mutex);
			ofs << s << '\n';
		}
		void close() {
			std::lock_guard<std::mutex> lk(caryon_io_mutex);
			if (ofs.is_open()) ofs.close();
		}
	};

	/**
	 * setDataName
	 * - 描述: 设置全局数据集名称（用于生成文件夹与文件名）。
	 * - 参数:
	 *    name: 数据集名称字符串
	 */
	inline void setDataName(const std::string &name) { datasetName = name; }

	/**
	 * setOutputDir
	 * - 描述: 设置输出目录（并尝试创建该目录）。
	 * - 参数:
	 *    dir: 输出目录路径
	 */
	inline void setOutputDir(const std::string &dir) { outputDirectory = dir; std::filesystem::create_directories(outputDirectory); }

	/**
	 * setMaxTimeGlobal
	 * - 描述: 设置判定 TLE 使用的全局最大运行时间（毫秒）。
	 * - 参数:
	 *    t: 毫秒数
	 */
	inline void setMaxTimeGlobal(int t) { maxRuntimeMs = t; }

	/**
	 * batchGenerate
	 * - 描述: 多线程批量生成测试用例文件，每个线程负责一段编号区间。
	 * - 参数:
	 *    cases: 总共要生成的用例数量（1..cases）
	 *    generator: 功能对象，签名为 (int caseIndex, FileCaseWriter &writer)
	 *    threads: 并发线程数（默认 1）
	 *    baseName: 基础文件名（默认使用 datasetName）
	 */
	inline void batchGenerate(int cases, std::function<void(int, FileCaseWriter &)> generator, int threads = 1, const std::string &baseName = "") {
		if (cases <= 0) return;
		if (threads <= 0) threads = 1;
		std::filesystem::path dir = outputDirectory.empty() ? std::filesystem::path("data-" + datasetName) : std::filesystem::path(outputDirectory);
		std::filesystem::create_directories(dir);
		auto worker = [&](int from, int to) {
			for (int i = from; i <= to; ++i) {
				FileCaseWriter writer;
				writer.setBaseName(baseName.empty() ? datasetName : baseName);
				writer.setOutputDir(dir.string());
				if (!writer.openCase(i)) continue;
				generator(i, writer);
				writer.close();
			}
			};
		int per = cases / threads;
		int rem = cases % threads;
		std::vector<std::thread> ths;
		int cur = 1;
		for (int t = 0; t < threads; ++t) {
			int start = cur;
			int cnt = per + (t < rem ? 1 : 0);
			int end = start + cnt - 1;
			if (cnt <= 0) break;
			ths.emplace_back(worker, start, end);
			cur = end + 1;
		}
		for (auto &t : ths) if (t.joinable()) t.join();
	}

	inline int runTestOnCase(int idx, const std::string &exeName = "test.exe", const std::string &baseName = "") {
		std::string inPath = getCaseInputPath(idx, baseName);
		std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
		std::string outPath = dir + "/" + (baseName.empty() ? datasetName : baseName) + std::to_string(idx) + ".out";
		std::string command = "\"" + exeName + "\" < \"" + inPath + "\" > \"" + outPath + "\"";
		return system(command.c_str());
	}

	inline void runTestsParallel(int start, int end, const std::string &exeName = "test.exe", int threads = 1, const std::string &baseName = "") {
		if (start > end) return;
		if (threads <= 0) threads = 1;
		int total = end - start + 1;
		int per = total / threads;
		int rem = total % threads;
		std::vector<std::thread> ths;
		int cur = start;
		auto worker = [&](int s, int e) {
			for (int i = s; i <= e; ++i) {
				runTestOnCase(i, exeName, baseName);
			}
			};
		for (int t = 0; t < threads; ++t) {
			int cnt = per + (t < rem ? 1 : 0);
			if (cnt <= 0) break;
			int s = cur;
			int e = cur + cnt - 1;
			ths.emplace_back(worker, s, e);
			cur = e + 1;
		}
		for (auto &th : ths) if (th.joinable()) th.join();
	}

	inline void generateAndRun(int cases, std::function<void(int, FileCaseWriter &)> generator, int genThreads = 1, int runThreads = 1, const std::string &exeName = "test.exe", const std::string &baseName = "") {
		batchGenerate(cases, generator, genThreads, baseName);
		runTestsParallel(1, cases, exeName, runThreads, baseName);
	}
} // namespace CaryonIO

// CaryonDebug: 用于生成调试输出、比较输出结果等
namespace CaryonDebug {
	std::stringstream __re;

	/**
	 * makeDebugFile
	 * - 描述: 生成指定编号区间的调试答案文件，并记录运行时间与返回码到 __re / runtimeLogStream。
	 * - 参数:
	 *    startIndex: 起始编号
	 *    endIndex: 结束编号
	 */
	void makeDebugFile(int startIndex, int endIndex) {
		std::string debugDir = "debug-" + datasetName;
		std::filesystem::create_directories(debugDir);
		runtimeLogStream.str("");
		runtimeLogStream.clear();
		__re.str("");
		__re.clear();
		std::string executable = "myprogram.exe";
		if (!std::filesystem::exists(executable)) {
			std::cerr << "Error: Executable file not found" << std::endl;
			return;
		}
		for (int i = startIndex; i <= endIndex; ++i) {
			std::string indexStr = std::to_string(i);
			std::string debugAnsPath = debugDir + "/" + datasetName + indexStr + ".ans";
			std::string inputPath = "data-" + datasetName + "/" + datasetName + indexStr + ".in";
			if (!std::filesystem::exists(inputPath)) {
				std::cerr << "Input file not found" << std::endl;
				continue;
			}
			std::string command = executable + " < \"" + inputPath + "\" > \"" + debugAnsPath + "\"";
			long double clock1 = clock();
			int returnCode = system(command.c_str());
			long double clock2 = clock();
			__re << returnCode << std::endl;
			runtimeMs = (clock2 - clock1) * 1000.0 / CLOCKS_PER_SEC;
			runtimeLogStream << runtimeMs << std::endl;
		}
	}

	/**
	 * compareFile
	 * - 描述: 将生成的 .out 与 .ans 文件进行比较并输出每个用例的评测结果。
	 * - 参数:
	 *    startIndex: 起始编号
	 *    endIndex: 结束编号
	 */
	void compareFile(int startIndex, int endIndex) {
		int acCount = 0;
		std::ofstream logfile("Debug.log", std::ios::trunc);
		if (!logfile.is_open()) {
			std::cerr << "Failed to open Debug.log" << std::endl;
			return;
		}
		logfile << "=== Debug Started ===\n";
		logfile << "Dataset: " << datasetName << "\n";
		logfile << "Test cases: " << startIndex << " to " << endIndex << "\n";
		logfile << "Time limit: " << maxRuntimeMs << " ms\n";
		runtimeLogStream.seekg(0);
		__re.seekg(0);
		for (int i = startIndex; i <= endIndex; i++) {
			std::string indexStr = std::to_string(i);
			std::string ansPath = "debug-" + datasetName + "/" + datasetName + indexStr + ".ans";
			std::string outPath = "data-" + datasetName + "/" + datasetName + indexStr + ".out";
			if (!std::filesystem::exists(ansPath)) {
				logfile << "Answer file not found.\n";
				continue;
			}
			if (!std::filesystem::exists(outPath)) {
				logfile << "Output file not found.\n";
				continue;
			}
			std::ifstream ansFile(ansPath);
			std::ifstream outFile(outPath);
			if (!ansFile.is_open() || !outFile.is_open()) {
				logfile << "Failed to open files.\n";
				continue;
			}
			bool filesEqual = true;
			std::string ansLine, outLine;
			int lineNum = 1;
			while (std::getline(ansFile, ansLine) && std::getline(outFile, outLine)) {
				ansLine.erase(ansLine.find_last_not_of(" \t\r\n") + 1);
				outLine.erase(outLine.find_last_not_of(" \t\r\n") + 1);
				if (ansLine != outLine) {
					filesEqual = false;
					break;
				}
				lineNum++;
			}
			if (filesEqual && (std::getline(ansFile, ansLine) || std::getline(outFile, outLine))) {
				filesEqual = false;
			}
			ansFile.close();
			outFile.close();
			long double runtime = 0;
			int returnCode = 0;
			if (!(runtimeLogStream >> runtime)) runtime = 0;
			if (!(__re >> returnCode)) returnCode = -1;
			std::string result;
			if (returnCode != 0) {
				result = "RE";
			}
			else if (runtime > maxRuntimeMs) {
				result = "TLE";
			}
			else if (filesEqual) {
				result = "AC";
				acCount++;
			}
			else {
				result = "WA";
			}
			logfile << "TestCase " << i << ", result: " << result << "\n";
		}
		logfile << "=== Debug Ended ===\n";
		double score = (acCount * 100.0) / (endIndex - startIndex + 1);
		logfile << "Score: " << std::fixed << std::setprecision(2) << score << "% ("
			<< acCount << "/" << (endIndex - startIndex + 1) << ")\n";
		logfile.close();
	}

	/**
	 * debug
	 * - 描述: 便捷接口，先生成 debug 文件再比较（并暂停）。
	 * - 参数:
	 *    startIndex: 起始编号
	 *    endIndex: 结束编号
	 */
	void debug(int startIndex, int endIndex) {
		makeDebugFile(startIndex, endIndex);
		compareFile(startIndex, endIndex);
	}
} // namespace CaryonDebug

#endif // CARYON_H
