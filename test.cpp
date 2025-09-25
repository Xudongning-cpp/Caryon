#include "testlib.h"
#include "caryon.h"

// 示例：自定义权重生成函数
int randWeight(int minW, int maxW) {
	return CaryonRandom::randInt(minW, maxW);
}

// 示例：自定义字符串生成函数
std::string customCharset() {
	return "!@#$%^&*()"; // 自定义字符集
}

int main(int argc, char* argv[]) {
	// 初始化 testlib 的随机数生成器
	registerGen(argc, argv, 1);

	// 设置数据集名称和输出目录
	datasetName = "DemoDataset";
	outputDirectory = "data-" + datasetName;

	// 设置最大运行时间
	maxRuntimeMs = 2000;

	// 输出常用常量
	std::cout << "=== 常用常量 ===\n";
	std::cout << "PI: " << CaryonConstants::PI << '\n';
	std::cout << "E: " << CaryonConstants::E << '\n';
	std::cout << "SQRT2: " << CaryonConstants::SQRT2 << '\n';
	std::cout << "SQRT3: " << CaryonConstants::SQRT3 << '\n';
	std::cout << "SQRT5: " << CaryonConstants::SQRT5 << '\n';
	std::cout << "LOG2: " << CaryonConstants::LOG2 << '\n';
	std::cout << "LOG10: " << CaryonConstants::LOG10 << '\n';
	std::cout << "小写字母表: " << CaryonConstants::ALPHABET_SMALL << '\n';
	std::cout << "大写字母表: " << CaryonConstants::ALPHABET_CAPITAL << '\n';
	std::cout << "数字表: " << CaryonConstants::NUMBERS << '\n';

	// 设置随机数种子
	CaryonRandom::seedAllRNG(time(0));

	// 生成随机数示例
	std::cout << "\n=== 随机数生成 ===\n";
	std::cout << "1~100 随机整数: " << CaryonRandom::randInt(1, 100) << '\n';
	std::cout << "0~1 随机浮点数: " << CaryonRandom::randDouble() << '\n';
	std::cout << "1~100 随机浮点数: " << CaryonRandom::randInt(0, 99) + CaryonRandom::randDouble() << '\n';
	std::cout << "1~100 随机长整数: " << CaryonRandom::randLong(1, 100) << '\n';

	// 生成随机字符串示例
	std::cout << "\n=== 随机字符串生成 ===\n";
	std::cout << "5位小写字母: " << CaryonRandom::randString(5, CaryonConstants::ALPHABET_SMALL) << '\n';
	std::cout << "5位大写字母: " << CaryonRandom::randString(5, CaryonConstants::ALPHABET_CAPITAL) << '\n';
	std::cout << "5位数字: " << CaryonRandom::randString(5, CaryonConstants::NUMBERS) << '\n';
	std::cout << "5位大小写字母: " << CaryonRandom::randString(5, CaryonConstants::ALPHABET_CAPITAL + CaryonConstants::ALPHABET_SMALL) << '\n';
	std::cout << "5位字母数字: " << CaryonRandom::randString(5, CaryonConstants::ALPHABET_CAPITAL + CaryonConstants::ALPHABET_SMALL + CaryonConstants::NUMBERS) << '\n';
	std::cout << "5位二进制串: " << CaryonRandom::randBinaryString(5) << '\n';
	std::cout << "5位回文串: " << CaryonRandom::randPalindrome(5, CaryonConstants::ALPHABET_CAPITAL + CaryonConstants::ALPHABET_SMALL + CaryonConstants::NUMBERS) << '\n';
	std::cout << "10位默认字符串: " << CaryonRandom::randString(10) << '\n';

	// 打乱操作示例
	std::cout << "\n=== 打乱操作 ===\n";
	std::vector<int> vec = { 1, 2, 3, 4, 5 };
	CaryonRandom::shuffleRange(vec.begin(), vec.end());
	std::cout << "打乱向量: ";
	for (int x : vec) std::cout << x << " ";
	std::cout << '\n';

	std::string s = "12345";
	std::cout << "打乱字符串: " << CaryonRandom::shuffleString(s) << '\n';

	// 生成排列示例
	int arr[5];
	CaryonRandom::randPermutation(arr, 5);
	std::cout << "随机排列: ";
	for (int i = 0; i < 5; ++i) std::cout << arr[i] << " ";
	std::cout << '\n';

	// 生成随机向量示例
	std::vector<int> randVec = CaryonRandom::randVector(5, 1, 100);
	std::cout << "随机向量: ";
	for (int x : randVec) std::cout << x << " ";
	std::cout << '\n';

	// 生成唯一抽样示例
	std::vector<int> uniqueSample = CaryonRandom::sampleUniqueRange(1, 100, 10);
	std::cout << "唯一抽样: ";
	for (int x : uniqueSample) std::cout << x << " ";
	std::cout << '\n';

	// 图生成示例
	std::cout << "\n=== 图生成示例 ===\n";

	// 生成随机图（10个顶点，15条边，权重1~10）
	auto graph = CaryonGraph::randGraph<int>(10, 15, 1, 10, randWeight);
	std::cout << "随机图:\n" << graph << '\n';

	// 生成随机树（10个顶点，最大度3，权重1~5）
	auto tree = CaryonGraph::randTree<int>(10, 3, 1, 5, randWeight);
	std::cout << "随机树:\n" << tree << '\n';

	// 生成DAG（8个顶点，10条边，权重1~5）
	auto dag = CaryonGraph::randDAG<int>(8, 10, 1, 5, randWeight);
	std::cout << "随机DAG:\n" << dag << '\n';

	// 生成完全图（5个顶点，权重1~3）
	auto complete = CaryonGraph::completeGraph<int>(5, 1, 3, randWeight);
	std::cout << "完全图:\n" << complete << '\n';

	// 生成二分图（左3右4，5条边，权重1~5）
	auto bipartite = CaryonGraph::bipartiteGraph<int>(3, 4, 5, 1, 5, randWeight);
	std::cout << "二分图:\n" << bipartite << '\n';

	// 生成网格图（3x3，权重1~2）
	auto grid = CaryonGraph::gridGraph<int>(3, 3, 1, 2, randWeight);
	std::cout << "网格图:\n" << grid << '\n';

	// 生成连通图（6个顶点，8条边，权重1~4）
	auto connected = CaryonGraph::connectedGraph<int>(6, 8, 1, 4, randWeight);
	std::cout << "连通图:\n" << connected << '\n';

	// 生成链状树和菊花图
	auto chain = CaryonGraph::randChain<int>(5, 1, 3, randWeight);
	auto star = CaryonGraph::randStar<int>(5, 1, 3, randWeight);
	std::cout << "链状树:\n" << chain << '\n';
	std::cout << "菊花图:\n" << star << '\n';

	// 图验证示例
	std::cout << "树验证: " << (CaryonGraph::isTree(tree) ? "是树" : "不是树") << '\n';
	std::cout << "DAG验证: " << (CaryonGraph::isDAG(dag) ? "是DAG" : "不是DAG") << '\n';
	std::cout << "重复边检测: " << (CaryonGraph::hasDuplicateEdges(graph) ? "有重复边" : "无重复边") << '\n';

	// 测试用例生成示例
	makein(1, 5) {
		int n = CaryonRandom::randInt(1, 1e5);
		int m = CaryonRandom::randInt(1, 1e5);
		CaryonIO::writeCase(n);
		CaryonIO::writeSpace();
		CaryonIO::writeCase(m);
		CaryonIO::writeSpace();
		CaryonIO::writeEndl();
	}
	CaryonIO::executeRangeStd(1, 5);
	
	// 调试模式：比较输出文件
	CaryonDebug::debug(1, 5);
	return 0;
}
