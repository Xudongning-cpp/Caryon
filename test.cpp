#include "caryon.h"
int main(int argc, char *argv[]) {
	datasetName = "TEST";
	rnd.setSeed(time(0));
	makein(1, 5) {
		int n = rnd.next(1, (int)1e5);
		int m = rnd.next(1, (int)1e2);
		int k = rnd.next(1, (int)1e2);
		int x = rnd.next(1, m);
		int y = rnd.next(x, m);
		int times = rnd.next(0, n / m - 1);
		FileIO::writeCase(n);
		FileIO::writeSpace();
		FileIO::writeCase(m);
		FileIO::writeSpace();
		FileIO::writeCase(k);
		FileIO::writeSpace();
		FileIO::writeCase(x);
		FileIO::writeSpace();
		FileIO::writeCase(y);
		FileIO::writeSpace();
		FileIO::writeCase(times);
	}
	FileIO::executeRangeStd(1, 5);
	Debug::debug(1, 5);
	return 0;
}
