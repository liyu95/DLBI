#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C" {
#include <getopt.h>
}
#include "util/exception.h"

struct opts {
    char input[255];
	char output[255];
	int size;
};

inline int GetOpts(int argc, char **argv, opts* opts_){

    static struct option longopts[] = {
        { "help",            no_argument,            NULL,              'h' },
        { "input",    	     required_argument,      NULL,              'i' },
		{ "output",    	     required_argument,      NULL,              'o' },
        { NULL,              0,                      NULL,               0  }
    };
	
    if( argc < 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1){
		EX_TRACE("[-i INPUT FILENAME][-o OUTPUT DIR][-s SIZE]\n");
		return -1;
    }
    
    int ch;
    while ((ch = getopt_long(argc, argv, "hi:o:s:", longopts, NULL))!= -1) {
        switch (ch){
        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;

        case 'h':
            EX_TRACE("[-i INPUT FILENAME][-o OUTPUT DIR]\n");
            return 0;

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'o': 
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
		case 's': 
        {
            std::istringstream iss(optarg);
            iss >> opts_->size;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
        case 0: 
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif