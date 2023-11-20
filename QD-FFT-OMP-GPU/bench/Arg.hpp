#pragma once

#include <iostream>
#include <set>
#include <string>

class Arg {
public:
    int argc;
    char **argv;
    std::set<std::string> args;
    Arg(int agrc, char **argv) : argc(agrc), argv(argv) {
        for (int i = 0; i < agrc; i++) {
            args.insert(std::string(argv[i]));
        }
    }
    Arg(int agrc, int min_argc, char **argv) : argc(agrc), argv(argv) {
        if (argc < min_argc) {
            std::cerr << "Not enough argments." << std::endl;
            exit(1);
        }

        for (int i = 0; i < agrc; i++) {
            args.insert(std::string(argv[i]));
        }
    }
    bool has(const std::string &arg) const { return (args.count(arg) > 0); }
};