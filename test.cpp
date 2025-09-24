#include "caryon.h"
int main() {
    datasetName = "Florr";
    maxRuntimeMs = 1000;
    CaryonRandom::seedAllRNG(time(0));
    makein(1, 20) {
        int n = CaryonRandom::randInt(1, 1e5);
        int m = CaryonRandom::randInt(1, 1e2);
        int k = CaryonRandom::randInt(1, 1e2);
        int x = CaryonRandom::randInt(1, m);
        int y = CaryonRandom::randInt(x, m);
        int times = CaryonRandom::randInt(0, n / m - 1);
        CaryonIO::writeCase(n);
        CaryonIO::writeSpace();
        CaryonIO::writeCase(m);
        CaryonIO::writeSpace();
        CaryonIO::writeCase(k);
        CaryonIO::writeSpace();
        CaryonIO::writeCase(x);
        CaryonIO::writeSpace();
        CaryonIO::writeCase(y);
        CaryonIO::writeSpace();
        CaryonIO::writeCase(times);
    }
    CaryonIO::executeRangeStd(1, 20);
    CaryonDebug::debug(1, 20);
    return 0;
}
