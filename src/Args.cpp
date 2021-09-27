#include "Args.hpp"


Args::Args(int argc, char** argv)
{
    // the first arg is the executable
    mArgs.reserve(argc);
    
    mArgs.push_back(*argv);
    
    for(uint32_t i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        mArgs.push_back(arg);
        if(arg.size() > 1 &&
            arg[0] == '-') {
            
            std::string id = arg.substr(1);
            uint64_t split = id.find_first_of('=');
            if(split != std::string::npos) {
                std::string val = id.substr(split + 1);
                id = id.substr(0, split);
                
                mMap.emplace(id, val);
            }else
            {
                mMap.emplace(id, "");
            }
        }
        
    }
}
