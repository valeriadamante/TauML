import ROOT
import numpy as np
import psutil
import os
import gc

pid = os.getpid()
py = psutil.Process(pid)

print(py.memory_info())

ROOT.gInterpreter.Declare("""
struct Data {
    Data(size_t n) : x(n), y(n) {}
    std::vector<float> x, y;
};
struct DataLoader {
    DataLoader(size_t _n) : n(_n) {}
    Data LoadNext()
    {
        return Data(n);
    }
private:
    size_t n;
};
""")

data_loader = ROOT.DataLoader(int(1e9))
for n in range(100):
    data = data_loader.LoadNext()
    print('Created', py.memory_info())
    x_np = np.asarray(data.x)
    print('Np allocated', py.memory_info())
    print("average", np.average(x_np))
    #del x_np
    #del x
    #print(gc.collect())
    #print("Gc collected", py.memory_info())
    #ROOT.store.Delete("x")
    #print("Root deleted", py.memory_info())
