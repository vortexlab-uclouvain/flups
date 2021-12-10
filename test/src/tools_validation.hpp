#ifndef _SRC_TOOLS_VALIDATION
#define _SRC_TOOLS_VALIDATION

#include <limits>
#include "flups.h"




double KernelOrder(FLUPS_GreenType kernel){
    switch (kernel) {
        case CHAT_2:
            printf("Warning, you shouldn't ask for the convergence order of a spectral kernel\n");
            return std::numeric_limits<double>::max();
            break;
        case LGF_2:
            return 2.0;
            break;
        case HEJ_2:
            return 2.0;
            break;
        case HEJ_4:
            return 4.0;
            break;
        case HEJ_6: 
            return 6.0;
            break;
        case HEJ_8:
            return 8.0;
            break;
        case HEJ_10:
            return 10.0;
            break;
        case HEJ_0:
            printf("Warning, you shouldn't ask for the convergence order of a spectral kernel\n");
            return std::numeric_limits<double>::max();
            break;
    }
    return 0; 
}

#endif //_SRC_TOOLS_VALIDATION