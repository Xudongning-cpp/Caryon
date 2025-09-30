#ifndef CARYON_H
#define CARYON_H

/*
 * Caryon Header
 * -------------
 * 主要功能：
 * - 随机数生成工具
 * - 图结构及生成工具
 * - 数据集生成/写入/运行/调试工具
 */

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
#include <atomic>
#include <sstream>
#include <iomanip>
#include <chrono>

#define makein(low, high) for (caseIndex = low, invocationCount = 0; caseIndex <= high; caseIndex++, invocationCount = 0)
#define CarYon 4

std::string datasetName;
std::string outputDirectory;
std::atomic<int> caseIndex{ 0 };
std::atomic<int> invocationCount{ 0 };
std::atomic<int> maxRuntimeMs{ 1000 };
long double runtimeMs;
std::atomic<bool> directoryCreated{ false };
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
	bool randomInitialized = false;
	bool useStdMT = true;
	int legacyIndex = 0;
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
		randomInitialized = true;
		legacyMT[0] = initialSeed;
		for (int i = 1; i < 624; i++) {
			int t = 1812433253 * (legacyMT[i - 1] ^ (legacyMT[i - 1] >> 30)) + i;
			legacyMT[i] = t & 0xffffffff;
		}
	}

	/*
	* legacyGenerate
	* - 生成 legacy MT 的下一批伪随机数（内部使用）
	*/
	inline void legacyGenerate() {
		for (int i = 0; i < 624; i++) {
			long long y = (legacyMT[i] & 0x80000000LL) + (legacyMT[(i + 1) % 624] & 0x7fffffffLL);
			legacyMT[i] = legacyMT[(i + 397) % 624] ^ (y >> 1);
			if (y % 2 == 1)
				legacyMT[i] ^= 2147483647LL;
		}
	}

	/**
	 * randInt
	 * - 描述: 返回一个非负的 32 位随机整数。
	 * - 返回: int 随机值范围 [0, 2^31-1)
	 */
	inline int randInt() {
		if (useStdMT) {
			return static_cast<int>(rng64() & 0x7fffffff);
		}
		if (!randomInitialized)
			seedLegacy(static_cast<int>(time(nullptr)));
		if (legacyIndex == 0)
			legacyGenerate();
		int y = static_cast<int>(legacyMT[legacyIndex]);
		y = y ^ (y >> 11);
		y = y ^ ((y << 7) & 1636928640);
		y = y ^ ((y << 15) & 1022730752);
		y = y ^ (y >> 18);
		legacyIndex = (legacyIndex + 1) % 624;
		return y & 0x7fffffff;
	}

	/**
	 * randLong
	 * - 描述: 返回一个随机的 64 位有符号整数（取自 mt19937_64）。
	 */
	inline long long randLong() {
		if (useStdMT) {
			return static_cast<long long>(rng64());
		}
		return ((static_cast<long long>(randInt()) << 31) + randInt());
	}

	/**
	 * randInt
	 * - 描述: 返回区间 [low, high] 之间的随机整数（含端点）。
	 * - 参数:
	 *    low: 下界（包含）
	 *    high: 上界（包含）
	 */
	inline int randInt(int low, int high) {
		if (low > high) throw std::invalid_argument("randInt(low, high): low > high");
		if (low == high) return low;
		if (useStdMT) {
			uint64_t r = rng64();
			return static_cast<int>(r % (uint64_t)(high - low + 1) + low);
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
		if (low > high) throw std::invalid_argument("randLong(low, high): low > high");
		if (low == high) return low;
		if (useStdMT) {
			uint64_t r = rng64();
			return static_cast<long long>(r % (uint64_t)(high - low + 1) + low);
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
	template <typename It>
	inline void shuffleRange(It first, It last) {
		if (useStdMT) {
			std::shuffle(first, last, rng64);
		}
		else {
			auto n = std::distance(first, last);
			for (decltype(n) i = n - 1; i > 0; --i) {
				auto j = randInt(0, static_cast<int>(i));
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
		for (int i = 0; i < length; i++) outArray[i] = i + 1;
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
		return static_cast<double>(randInt()) / static_cast<double>(0x7fffffff);
	}

	/**
	 * seedAllRNG
	 * - 描述: 使用相同的种子同时初始化标准 MT 与 legacy MT。
	 * - 参数:
	 *    seedValue: 要使用的种子
	 */
	inline void seedAllRNG(unsigned long long seedValue) {
		seedStd(seedValue);
		seedLegacy(static_cast<int>(seedValue & 0xffffffff));
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
		if (count > range) throw std::invalid_argument("sampleUniqueRange: count > range");
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
			result += charset[randInt(0, static_cast<int>(charset.size()) - 1)];
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
			char c = charset[randInt(0, static_cast<int>(charset.size()) - 1)];
			result[i] = c;
			result[length - 1 - i] = c;
		}
		if (length % 2 == 1) {
			result[half] = charset[randInt(0, static_cast<int>(charset.size()) - 1)];
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
		shuffleRange(chars.begin(), chars.end());
		return std::string(chars.begin(), chars.end());
	}
} // namespace CaryonRandom

// CaryonGraph: 图结构与生成器
namespace CaryonGraph {
	/**
	 * Edge<T>
	 * - 描述: 图的边表示，包含目标顶点 v 与权重 w。
	 * - 成员:
	 *    v: 目标顶点编号
	 *    w: 边的权重（模板类型）
	 */
	template <typename Type>
	struct Edge {
		int v;
		Type w;
		bool operator<(const Edge &rw) const {
			return w > rw.w;
		}
	};

	/**
	 * Adjacency<T>
	 * - 描述: 邻接表项，包含若干边。
	 */
	template <typename Type>
	struct Adjacency {
		std::vector< Edge<Type> > edges;
	};

	/**
	 * Graph<T>
	 * - 描述: 简单的无向/有向图邻接表实现（顶点从 1 开始）。
	 * - 成员:
	 *    n: 当前最大顶点编号
	 *    m: 边数
	 *    adj: 邻接表（按索引存储 Adjacency）
	 */
	template <typename Type>
	struct Graph {
		int n = 0, m = 0;
		std::vector< Adjacency<Type> > adj;
		Graph() {}

		/**
		 * ensureSize
		 * - 描述: 确保 adj 的大小至少能包含 targetIndex（按顶点编号索引）。
		 * - 参数:
		 *    targetIndex: 需要确保的最大索引
		 */
		void ensureSize(int targetIndex) {
			Adjacency<Type> updatemp;
			updatemp.edges.clear();
			while ((int)adj.size() <= targetIndex) adj.push_back(updatemp);
		}

		/**
		 * addEdge
		 * - 描述: 向图中添加一条从 from 到 to 的边，权重为 weight。
		 * - 参数:
		 *    from: 边的起始顶点编号
		 *    to: 边的目标顶点编号
		 *    weight: 边的权重
		 */
		void addEdge(int from, int to, Type weight) {
			if (from <= 0 || to <= 0) return;
			if (from == to) return;
			ensureSize(std::max(from, to));
			n = std::max(n, from);
			n = std::max(n, to);
			m++;
			Edge<Type> tmp{ to, weight };
			adj[from].edges.push_back(tmp);
		}

		/**
		 * isConnected
		 * - 描述: 简单 BFS 连通性检测（从随机顶点开始）。
		 * - 返回: 若图的所有顶点都可达则返回 true。
		 */
		bool isConnected() {
			if (n == 0) return true;
			const int vstsize = n + 1;
			int vstn = 0;
			std::vector<bool> visited(vstsize, false);
			std::queue<int> q;
			int start = 1;
			for (int i = 1; i <= n; ++i) {
				if (!adj[i].edges.empty()) {
					start = i;
					break;
				}
			}
			visited[start] = true;
			vstn = 1;
			q.push(start);
			while (!q.empty()) {
				int cur = q.front();
				q.pop();
				for (const auto &e : adj[cur].edges) {
					if (!visited[e.v]) {
						visited[e.v] = true;
						vstn++;
						q.push(e.v);
					}
				}
			}
			for (int i = 1; i <= n; ++i) {
				if (!adj[i].edges.empty() && !visited[i]) return false;
			}
			return true;
		}
	};

	// 图输出流
	template <typename Type>
	std::ostream &operator<<(std::ostream &os, const Graph<Type> &c) {
		os << c.n << ' ' << c.m << '\n';
		for (int i = 1; i <= c.n; i++) {
			for (const auto &edge : c.adj[i].edges)
				os << i << ' ' << edge.v << ' ' << edge.w << '\n';
		}
		return os;
	}

	// 图合并
	template <typename Type>
	Graph<Type> operator+(Graph<Type> a, Graph<Type> b) {
		Graph<Type> ret = a;
		ret.ensureSize(std::max(a.n, b.n));
		for (int i = 1; i <= b.n; i++) {
			for (const auto &edge : b.adj[i].edges) {
				ret.addEdge(i, edge.v, edge.w);
			}
		}
		ret.n = std::max(a.n, b.n);
		ret.m = a.m + b.m;
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
	template <typename Type>
	Graph<Type> randGraph(int nodeCount, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		Graph<Type> ret;
		ret.n = nodeCount;
		ret.ensureSize(nodeCount);
		std::unordered_set<uint64_t> edge_set;
		int edges_added = 0;
		auto encode = [](int u, int v) { return ((uint64_t)u << 32) | v; };
		while (edges_added < edgeCount) {
			int u = CaryonRandom::randInt(1, nodeCount);
			int v = CaryonRandom::randInt(1, nodeCount);
			if (u == v) continue;
			if (edge_set.count(encode(u, v))) continue;
			edge_set.insert(encode(u, v));
			ret.addEdge(u, v, randValueFunc(minWeight, maxWeight));
			edges_added++;
		}
		ret.m = edges_added;
		return ret;
	}

	/**
	 * randDAG
	 * - 描述: 随机生成有向无环图。
	 * - 参数:
	 *    nodeCount: 顶点数量
	 *    edgeCount: 边数量
	 *    minWeight: 权重最小值
	 *    maxWeight: 权重最大值
	 *    randValueFunc: 随机权重生成函数
	 */
	template <typename Type>
	Graph<Type> randDAG(int nodeCount, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		Graph<Type> ret;
		if (nodeCount < 2) return ret;
		ret.n = nodeCount;
		ret.ensureSize(nodeCount);
		std::unordered_set<uint64_t> edge_set;
		int edges_added = 0;
		auto encode = [](int u, int v) { return ((uint64_t)u << 32) | v; };
		while (edges_added < edgeCount) {
			int u = CaryonRandom::randInt(1, nodeCount - 1);
			int v = CaryonRandom::randInt(u + 1, nodeCount);
			if (u >= v) continue;
			if (edge_set.count(encode(u, v))) continue;
			edge_set.insert(encode(u, v));
			ret.addEdge(u, v, randValueFunc(minWeight, maxWeight));
			edges_added++;
		}
		ret.m = edges_added;
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
	template <typename Type>
	Graph<Type> randTree(int nodeCount, int maxChildren, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		Graph<Type> ret;
		if (nodeCount == 0) return ret;
		ret.n = nodeCount;
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
			std::swap(t[CaryonRandom::randInt(0, static_cast<int>(t.size()) - 1)], t[t.size() - 1]);
			t.back().sonCount++;
			if (t.back().sonCount == maxChildren) t.pop_back();
			ret.addEdge(i, t.back().id, randValueFunc(minWeight, maxWeight));
			updatemp.id = i;
			updatemp.sonCount = 0;
			t.push_back(updatemp);
		}
		ret.m = nodeCount - 1;
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
	template<typename Type>
	Graph<Type> completeGraph(int nodeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		Graph<Type> graph;
		graph.n = nodeCount;
		graph.ensureSize(nodeCount);
		for (int i = 1; i <= nodeCount; ++i) {
			for (int j = i + 1; j <= nodeCount; ++j) {
				graph.addEdge(i, j, randValueFunc(minWeight, maxWeight));
				graph.addEdge(j, i, randValueFunc(minWeight, maxWeight));
			}
		}
		graph.m = nodeCount * (nodeCount - 1);
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
	template<typename Type>
	Graph<Type> bipartiteGraph(int leftSize, int rightSize, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		Graph<Type> graph;
		graph.n = leftSize + rightSize;
		graph.ensureSize(graph.n);
		std::vector<std::pair<int, int>> edges;
		for (int i = 1; i <= leftSize; ++i) {
			for (int j = leftSize + 1; j <= leftSize + rightSize; ++j) {
				edges.emplace_back(i, j);
			}
		}
		CaryonRandom::shuffleRange(edges.begin(), edges.end());
		edgeCount = std::min(edgeCount, static_cast<int>(edges.size()));
		for (int i = 0; i < edgeCount; ++i) {
			auto [u, v] = edges[i];
			graph.addEdge(u, v, randValueFunc(minWeight, maxWeight));
		}
		graph.m = edgeCount;
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
	template<typename Type>
	Graph<Type> gridGraph(int rows, int cols, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		Graph<Type> graph;
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
	template<typename Type>
	Graph<Type> connectedGraph(int nodeCount, int edgeCount, Type minWeight, Type maxWeight, Type(*randValueFunc)(Type, Type)) {
		if (edgeCount < nodeCount - 1) edgeCount = nodeCount - 1;
		auto tree = randTree(nodeCount, nodeCount - 1, minWeight, maxWeight, randValueFunc);
		std::unordered_set<uint64_t> edge_set;
		auto encode = [](int u, int v) { return ((uint64_t)u << 32) | v; };
		for (int i = 1; i <= tree.n; ++i)
			for (const auto &e : tree.adj[i].edges)
				edge_set.insert(encode(i, e.v));
		int remainingEdges = edgeCount - (nodeCount - 1);
		int added = 0;
		while (added < remainingEdges) {
			int u = CaryonRandom::randInt(1, nodeCount);
			int v = CaryonRandom::randInt(1, nodeCount);
			if (u == v) continue;
			if (edge_set.count(encode(u, v))) continue;
			tree.addEdge(u, v, randValueFunc(minWeight, maxWeight));
			edge_set.insert(encode(u, v));
			added++;
		}
		tree.m = edgeCount;
		return tree;
	}

	/**
	 * randChain
	 * - 生成链状树
	 */
	template<typename Type>
	Graph<Type> randChain(int nodeCount, Type minWeight, Type maxWeight, Type(*randFunc)(Type, Type)) {
		Graph<Type> ret;
		ret.n = nodeCount;
		ret.ensureSize(nodeCount);
		for (int i = 2; i <= nodeCount; ++i) {
			ret.addEdge(i - 1, i, randFunc(minWeight, maxWeight));
		}
		ret.m = nodeCount - 1;
		return ret;
	}

	/**
	 * randStar
	 * - 生成菊花图
	 */
	template<typename Type>
	Graph<Type> randStar(int nodeCount, Type minWeight, Type maxWeight, Type(*randFunc)(Type, Type)) {
		Graph<Type> ret;
		ret.n = nodeCount;
		ret.ensureSize(nodeCount);
		for (int i = 2; i <= nodeCount; ++i) {
			ret.addEdge(1, i, randFunc(minWeight, maxWeight));
		}
		ret.m = nodeCount - 1;
		return ret;
	}

	/**
	 * isTree
	 * - 描述: 验证图是否为树（连通且边数 = 顶点数 - 1）。
	 * - 参数:
	 *    graph: 要验证的图对象
	 * - 返回: 如果是树则返回 true
	 */
	template<typename Type>
	bool isTree(Graph<Type> graph) {
		return graph.isConnected() && (graph.m == graph.n - 1);
	}

	/**
	 * isDAG
	 * - 描述: 使用拓扑排序验证图是否为有向无环图。
	 * - 参数:
	 *    graph: 要验证的图对象
	 * - 返回: 如果是 DAG 则返回 true
	 */
	template<typename Type>
	bool isDAG(const Graph<Type> &graph) {
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
	template<typename Type>
	bool hasDuplicateEdges(const Graph<Type> &graph) {
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
		std::lock_guard<std::mutex> lk(caryon_io_mutex);
		std::stringstream cci;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << std::setw(3) << std::setfill('0') << caseIndex;
			std::string name = "data-" + datasetName + "/" + datasetName + cci.str() + ".in";
			freopen(name.c_str(), "w", stdout);
			freopen(name.c_str(), "r", stdin);
		}
		std::cout << graph;
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
		std::lock_guard<std::mutex> lk(caryon_io_mutex);
		std::stringstream cci;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << std::setw(3) << std::setfill('0') << caseIndex;
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
		std::lock_guard<std::mutex> lk(caryon_io_mutex);
		std::stringstream cci;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << std::setw(3) << std::setfill('0') << caseIndex;
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
		std::lock_guard<std::mutex> lk(caryon_io_mutex);
		std::stringstream cci;
		if (invocationCount == 0) {
			if (!directoryCreated) {
				std::filesystem::create_directories("data-" + datasetName);
				directoryCreated = true;
			}
			cci << std::setw(3) << std::setfill('0') << caseIndex;
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
		std::lock_guard<std::mutex> lk(caryon_io_mutex);
		freopen("CON", "w", stdout);
		freopen("CON", "r", stdin);
		std::stringstream aa;
		aa << std::setw(3) << std::setfill('0') << caseNumber;
		std::string name = "data-" + datasetName + "/" + datasetName + aa.str() + ".in";
		freopen(name.c_str(), "r", stdin);
		std::string outname = "data-" + datasetName + "/" + datasetName + aa.str() + ".out";
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
			freopen("CON", "r", stdin);
			freopen("CON", "w", stdout);
		}
	}

	/**
	 * closeStreams
	 * - 描述: 恢复标准输入输出流。
	 */
	inline void closeStreams() {
		std::lock_guard<std::mutex> lk(caryon_io_mutex);
		freopen("CON", "r", stdin);
		freopen("CON", "w", stdout);
	}

	/**
	 * CaseFileWriter
	 * - 描述: 线程安全的文件写入器，支持按用例编号写入文件。
	 */
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
				std::stringstream aa;
				aa << std::setw(3) << std::setfill('0') << idx;
				std::string fname = outputDirectory + "/" + (baseName.empty() ? datasetName : baseName) +
					aa.str() + ".in";
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
	inline void setOutputDirectory(const std::string &dir) {
		outputDirectory = dir;
		std::filesystem::create_directories(outputDirectory);
	}

	/**
	 * setMaxRuntimeGlobal
	 * - 描述: 设置判定 TLE 使用的全局最大运行时间（毫秒）。
	 * - 参数:
	 *    t: 毫秒数
	 */
	inline void setMaxRuntimeGlobal(int t) { maxRuntimeMs = t; }

	/**
	 * getCaseInputPath
	 * - 描述: 获取指定编号用例的输入文件路径。
	 * - 参数:
	 *    idx: 用例编号
	 *    baseName: 可选基础文件名
	 */
	inline std::string getCaseInputPath(int idx, const std::string &baseName = "") {
		std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
		std::stringstream aa;
		aa << std::setw(3) << std::setfill('0') << idx;
		std::string name = dir + "/" + (baseName.empty() ? datasetName : baseName) +
			aa.str() + ".in";
		return name;
	}

	/**
	 * runExecutableOnCase
	 * - 描述: 在指定用例上运行可执行程序。
	 * - 参数:
	 *    idx: 用例编号
	 *    exeName: 可执行文件名
	 *    baseName: 可选基础文件名
	 * - 返回: 系统调用返回值
	 */
	inline int runExecutableOnCase(int idx, const std::string &exeName = "test.exe", const std::string &baseName = "") {
		std::string inPath = getCaseInputPath(idx, baseName);
		std::string dir = outputDirectory.empty() ? ("data-" + datasetName) : outputDirectory;
		std::stringstream aa;
		aa << std::setw(3) << std::setfill('0') << idx;
		std::string outPath = dir + "/" + (baseName.empty() ? datasetName : baseName) +
			aa.str() + ".out";
		std::string command = "\"" + exeName + "\" < \"" + inPath + "\" > \"" + outPath + "\"";
		return system(command.c_str());
	}
} // namespace CaryonIO

// CaryonDebug: 调试与输出比较
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
		if (startIndex > endIndex) return;
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
			std::stringstream aa;
			aa << std::setw(3) << std::setfill('0') << i;
			std::string debugAnsPath = debugDir + "/" + datasetName + aa.str() + ".ans";
			std::string inputPath = "data-" + datasetName + "/" + datasetName + aa.str() + ".in";
			if (!std::filesystem::exists(inputPath)) {
				std::cerr << "Input file not found" << std::endl;
				continue;
			}
			std::string command = executable + " < \"" + inputPath + "\" > \"" + debugAnsPath + "\"";
			auto clock1 = std::chrono::high_resolution_clock::now();
			int returnCode = system(command.c_str());
			auto clock2 = std::chrono::high_resolution_clock::now();
			__re << returnCode << std::endl;
			runtimeMs = std::chrono::duration<long double, std::milli>(clock2 - clock1).count();
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
		if (startIndex > endIndex) return;
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
			std::stringstream aa;
			aa << std::setw(3) << std::setfill('0') << i;
			std::string ansPath = "debug-" + datasetName + "/" + datasetName + aa.str() + ".ans";
			std::string outPath = "data-" + datasetName + "/" + datasetName + aa.str() + ".out";
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
				auto ansTrim = ansLine.find_last_not_of(" \t\r\n");
				auto outTrim = outLine.find_last_not_of(" \t\r\n");
				if (ansTrim != std::string::npos) ansLine.erase(ansTrim + 1);
				else ansLine.clear();
				if (outTrim != std::string::npos) outLine.erase(outTrim + 1);
				else outLine.clear();
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
	 * - 描述: 便捷接口，先生成 debug 文件再比较。
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
