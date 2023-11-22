#pragma once

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

class Arg {
public:
    int _argc;
    char **_argv;
    std::map<std::string, uint64_t> _arg_idx_map;
    std::map<std::string, bool> _flags_map;
    std::map<std::string, uint64_t> _options_map;
    Arg(int argc, char *argv[], int min_argc, std::vector<std::string> flags, std::vector<std::string> options, std::vector<uint64_t> options_default) : _argc(argc), _argv(argv) {
        if (argc < min_argc) {
            std::cerr << "Not enough arguments" << std::endl;
        }

        std::set<std::string> argv_set;
        for (int i = 0; i < argc; i++) {
            argv_set.insert(argv[i]);
            _arg_idx_map[argv[i]] = i;
        }

        if (argc < min_argc || argv_set.find("--help") != argv_set.end()) {
            std::cerr << "min_argc: " << min_argc << std::endl;
            std::cerr << "flags: " << std::endl;
            for (int i = 0; i < flags.size(); i++) {
                std::cerr << flags[i] << std::endl;
            }
            std::cerr << "options: " << std::endl;
            for (int i = 0; i < options.size(); i++) {
                std::cerr << options[i] << ": default " << options_default[i] << std::endl;
            }
            exit(1);
        }

        for (int i = 0; i < flags.size(); i++) {
            _flags_map[flags[i]] = (argv_set.find(flags[i]) != argv_set.end());
        }

        for (int i = 0; i < options.size(); i++) {
            std::string option = options[i];
            if (argv_set.find(option) != argv_set.end()) {
                int idx              = _arg_idx_map[option];
                _options_map[option] = atoi(argv[idx + 1]);
            } else {
                _options_map[option] = options_default[i];
            }
        }
    }

    int argc() {
        return _argc;
    }

    char *argv(int i) {
        return _argv[i];
    }

    bool has_flag(std::string flag) {
        return _flags_map[flag];
    }

    uint64_t get_option(std::string option) {
        return _options_map[option];
    }
};