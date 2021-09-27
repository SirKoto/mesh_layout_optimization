
#include <map>
#include <string>
#include <vector>


class Args
{
public:
    
    Args(int argc, char** argv);
    
    bool has(const std::string& id) const {
        return mMap.count(id) != 0;
    }
    
    const std::string& get(const std::string& id) const {
        return mMap.at(id);
    }

    const std::string& get(uint32_t i) const{
        return mArgs[i];
    }

    size_t numArgs() const { return mArgs.size(); }
    
    
private:
    std::vector<std::string> mArgs;
    std::map<std::string, std::string> mMap;
    
};
